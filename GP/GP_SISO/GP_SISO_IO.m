classdef GP_SISO_IO < handle
    %GP_SISO_IO Gaussian Process model for SISO input-output systems.
    %
    %   This class implements a Gaussian Process (GP) regression framework
    %   for single-input single-output (SISO) systems. It allows training
    %   on multiple trajectories, constructing regression matrices, and
    %   making predictions for new inputs.
    
    properties
        GP          % Gaussian Process object
        y_cell      % Cell containing training output Data (from each trial)
        u_cell      % Cell containing training input Data (from each trial)
        V           % Design matrix (containing u regression vectors)
        z           % Target Vector (containing y target Values)
        N           % Length of input Trajectory
        J           % Number of training trials

        sigmaF2     % Squared kernel parameter sigma_F^2
        sigmaL2_inv % Squared inverse kernel parameter 1/l^2
        alpha       % Constant term of the GP a = [K + sigma^2*I]^-1 * y
        V_row_norm2 % Euklidean Norm of each row of obj.V
        V_transp    % Transposed design Matrix
        L_chol      % Cholesky Factor of (K + sigma_N^2*I)

        sigma_n     % Optional: Noise standard deviation sigma_n
    end
    
    methods
        function obj = GP_SISO_IO(sigma_n)
            %GP_SISO_IO Constructor.
            %
            %   Optional:
            %       sigma_n : Noise standard deviation sigma_n (Skalar).

            if nargin < 1 || isempty(sigma_n)
                obj.sigma_n = [];
            else
                obj.sigma_n = sigma_n;
            end
        end
        
        function train_GP_model(obj, y_cell, u_cell)
            %train_GP_model Train the Gaussian Process model using provided data.
            %
            %   Inputs:
            %       y_cell : Cell array of output trajectories, one per trial
            %       u_cell : Cell array of input trajectories, one per trial

            % Save training runs
            obj.y_cell = y_cell;
            obj.u_cell = u_cell;

            % Get length and number of training trials
            obj.N = length(obj.u_cell{1});
            obj.J = length(obj.u_cell);

            % Design Matrix V
            obj.constr_design_matrix();

            % Target Vector z
            obj.constr_target_vector();

            % Train GP (only one lengthscale)
            if ~isempty(obj.sigma_n)
                % Train GP and just optimize the kernel parameters (exclude
                % sigma_n)
                obj.GP = fitrgp( ...
                    obj.V, obj.z, ...
                    'BasisFunction','none', ...
                    'KernelFunction','squaredexponential', ...
                    'Standardize', false, ...
                    'SigmaLowerBound', 1.e-6, ...
                    'Sigma', obj.sigma_n, ...
                    'ConstantSigma', true);
            else
                % Train GP and optimize all Hyperparameters (including
                % sigma_n)
                obj.GP = fitrgp( ...
                    obj.V, obj.z, ...
                    'BasisFunction','none', ...
                    'KernelFunction','squaredexponential', ...
                    'Standardize', false);
            end
        end

        function constr_design_matrix(obj)
            %constr_design_matrix Construct the GP design (regressor) matrix.

            % Empty design matrix
            obj.V = zeros(obj.J*obj.N, obj.N);

            % Design Matrix
            for i = 1:obj.J
                % Toeplitz Matrix
                V_temp = toeplitz([0; obj.u_cell{i}(1:end-1)], zeros(1, obj.N));

                % Put Toeplitz into Design Matrix
                obj.V(((i-1)*obj.N + 1):(i*obj.N), :) = V_temp;
            end
    
        end

        function constr_target_vector(obj)
            %constr_target_vector Construct target vector for GP regression.

            % Empty Target vector
            obj.z = zeros(obj.J*obj.N, 1);

            % Target Vector
            for i = 1:obj.J
                obj.z(((i-1)*obj.N + 1):(i*obj.N), :) = obj.y_cell{i};
            end
        end

        function [y_pred, y_std] = predict_trajectory(obj, u_test)
            %predict_GP Predict system output for a new input trajectory.

            u_test = u_test(:);  % Ensure column vector
            
            % Construct test Matrix
            V_test = obj.constr_test_matrix(u_test);

            % Prediction
            [y_pred, y_std] = predict(obj.GP, V_test);
        end

        function V_test = constr_test_matrix(obj, u_test)
            %constr_test_matrix Construct test regression matrix for prediction or linearization.
            %
            %   Inputs:
            %       u_test : Column vector of input samples (test trajectory)
            %
            %   Outputs:
            %       V_test : Toeplitz-structured regression matrix built from u_test

            % Toeplitz test Matrix
            V_test = toeplitz([0; u_test(1:end-1)], zeros(1, obj.N));
        end

        function P = linearize_at_given_trajectory(obj, u_lin)
            %linearize_at_given_trajectory Compute local linearization (Jacobian) of GP at a trajectory.
            %
            %   Inputs:
            %       u_lin : Input trajectory (column vector) around which the GP is linearized
            %
            %   Outputs:
            %       P : Jacobian (linearization) matrix of size N×N, mapping input variations
            %           Δu to output variations Δy near the trajectory u_lin

            % Empty jacobi matrix
            P = zeros(obj.N, obj.N);

            % Construct linearisation Matrix
            V_lin = obj.constr_test_matrix(u_lin);

            % Iteratively compute gradients
            for i = 1:obj.N
                % Extract regression Vector
                vn = V_lin(i, :);

                % Gradient of the GP w.r.t vn
                dy_dv = obj.gradient_wrt_regression_vector(vn);

                % Update jacobi matrix
                obj.J = i - 1;
                for j = 1:obj.J
                    P(i, j) = dy_dv(i - j);
                end
            end
        end

        function dy_dv = gradient_wrt_regression_vector(obj, vn)
            %gradient_wrt_regression_vector Compute GP output gradient w.r.t. a regression vector.
            %
            %   Inputs:
            %       vn : Current regression vector (column or row vector)
            %
            %   Outputs:
            %       dy_dv : Gradient of GP output with respect to vn (column vector)

            % Kernel Parameters
            k_params = obj.GP.KernelInformation.KernelParameters;
            sigmaL = k_params(1);
            sigmaF = k_params(2);

            % Receive alpha vector a = [K + sigma^2*I]^-1 * y
            alpha_vec = obj.GP.Alpha;

            N_alpha = length(alpha_vec);
            N_vn = length(vn);

            % Derivative of the Kernel function dk(vn,V)/dvn
            dk_dvn = zeros(N_vn, N_alpha);
            x1 = vn(:);
            for i = 1:N_alpha
                x2 = obj.V(i, :);
                x2 = x2(:);
                dk_dvn(:, i) = -1/(sigmaL^2) * (x1 - x2) * obj.squared_exponantial_kernel(x1, x2, sigmaL, sigmaF);
            end

            % Gradient w.r.t regression vector
            dy_dv = dk_dvn * alpha_vec;
        end

        function [P, Var_P] = linearize_at_given_trajectory_fast(obj, u_lin)
            %linearize_at_given_trajectory_fast Compute local linearization (Jacobian) of GP at a trajectory.
            %
            %   Inputs:
            %       u_lin : Input trajectory (column vector) around which the GP is linearized
            %
            %   Outputs:
            %       P : Jacobian (linearization) matrix of size N×N, mapping input variations
            %           Δu to output variations Δy near the trajectory u_lin

            % Empty jacobi matrix
            P = zeros(obj.N, obj.N);
            Var_P = zeros(obj.N, obj.N);

            % Construct linearisation Matrix
            V_lin = obj.constr_test_matrix(u_lin);

            % Kernel Parameters
            k_params = obj.GP.KernelInformation.KernelParameters;
            obj.sigmaL2_inv = 1/(k_params(1)^2);
            obj.sigmaF2 = k_params(2)^2;

            % Receive alpha vector a = [K + sigma^2*I]^-1 * y
            obj.alpha = obj.GP.Alpha;

            % Receive Cholesky Factor of (K + sigma_N^2*I)
            obj.L_chol = obj.GP.Impl.LFactor;

            % Precompute euklidean Norm of each row of obj.V and transposed
            % design matrix (each column is one regression Vector)
            obj.V_row_norm2 = sum(obj.V.^2, 2);
            obj.V_transp = obj.V.'; % .' - no complex transposition

            % Iteratively compute gradients
            for i = 1:obj.N
                % Extract regression Vector
                vn = V_lin(i, :);

                % Gradient and Variance of the GP w.r.t vn
                dy_dv = obj.gradient_wrt_regression_vector_fast(vn);
                Var_dy_dv = obj.gradient_variance_wrt_regression_vector_fast(vn);

                % Local variables for parfor
                rowP     = zeros(1, obj.N);
                rowVar_P = zeros(1, obj.N);

                % Update jacobi matrix
                if i > 1
                    rowP(1:(i-1))     = dy_dv((i-1):-1:1);
                    rowVar_P(1:(i-1)) = Var_dy_dv((i-1):-1:1);
                end

                % Update Row (parfor)
                P(i, :) = rowP;
                Var_P(i, :) = rowVar_P;
            end
        end

        function dy_dv = gradient_wrt_regression_vector_fast(obj, vn)
            %gradient_wrt_regression_vector_fast Compute GP output gradient w.r.t. a regression vector.
            %
            %   Inputs:
            %       vn : Current regression vector (column or row vector)
            %
            %   Outputs:
            %       dy_dv : Gradient of GP output with respect to vn (column vector)

            % Column Vector
            vn = vn(:);

            % Calculate euklidean norm between vn and every regression
            % vector in V
            vn_norm2 = sum(vn.^2);
            Vv = obj.V*vn;          % vector containing V_i'*v
            r2 = vn_norm2 + obj.V_row_norm2 - 2*Vv;

            % Calculate kernel Vector (element-wise)
            k = obj.sigmaF2 * exp(-0.5 * obj.sigmaL2_inv * r2);

            % Calculate gradient
            w = obj.alpha .* k;
            dy_dv = -obj.sigmaL2_inv * (sum(w)*vn - obj.V_transp*w);
        end

        function Var_dy_dv = gradient_variance_wrt_regression_vector_fast(obj, vn)
            %gradient_variance_wrt_regression_vector_fast Compute GP output gradient variance w.r.t. a regression vector.
            %
            %   Inputs:
            %       vn : Current regression vector (column or row vector)
            %
            %   Outputs:
            %       Var_dy_dv : Gradient variance of GP output with respect to vn (column vector)

            % Column Vector
            vn = vn(:);

            % Calculate euklidean norm between vn and every regression
            % vector in V
            vn_norm2 = sum(vn.^2);
            Vv = obj.V*vn;          % vector containing V_i'*v
            r2 = vn_norm2 + obj.V_row_norm2 - 2*Vv;

            % Calculate kernel Vector (element-wise)
            k = obj.sigmaF2 * exp(-0.5 * obj.sigmaL2_inv * r2);

            % Difference Matrix, Diff(i,j) = V(i,j) - vn(j) (every row
            % contains the element-wise differences between V(i, :) and vn)
            Diff = obj.V - vn.';

            % g = dk(V, vn) / dvn(j) => G = [g1, ..., gN]
            G = obj.sigmaL2_inv * (Diff .* k);

            % Solve B = K^-1*G
            B = obj.L_chol \ G;
            B = obj.L_chol.' \ B;

            % Calculate Data_Terms d(i) = G(:, i)^T * B(:, i) => Summation
            % in direction 1 (vertical)
            d = sum(G .* B, 1);

            % Calculate varianve
            Var_dy_dv = obj.sigmaF2 * obj.sigmaL2_inv * ones(obj.N, 1) - d(:);
        end

        function [P, Var_P, Cov_dy_dv_cell] = approx_linearisation_at_given_trajectory(obj, u_lin)
            %approx_linearisation_at_given_trajectory Compute local approximated linearization (Jacobian) of GP at a trajectory.
            %
            %   Inputs:
            %       u_lin : Input trajectory (column vector) around which the GP is linearized
            %
            %   Outputs:
            %       P : Jacobian (linearization) matrix of size N×N, mapping input variations
            %           Δu to output variations Δy near the trajectory u_lin

            % Empty jacobi matrix (and Variance)
            P = zeros(obj.N, obj.N);
            Var_P = zeros(obj.N, obj.N);

            % Cell Array for Covariances
            Cov_dy_dv_cell = cell(obj.N, 1);

            % Construct linearisation Matrix
            V_lin = obj.constr_test_matrix(u_lin);

            % Kernel Parameters
            k_params = obj.GP.KernelInformation.KernelParameters;
            obj.sigmaL2_inv = 1/(k_params(1)^2);
            obj.sigmaF2 = k_params(2)^2;

            % Receive Cholesky Factor of (K + sigma_N^2*I)
            obj.L_chol = obj.GP.Impl.LFactor;

            % Step-size
            h = 1e-1;

            % Iteratively compute gradients
            for i = 1:obj.N
                % Extract regression Vector
                vn = V_lin(i, :);

                % Gradient of the GP w.r.t vn
                [dy_dv, Cov_dy_dv] = obj.approx_gradient_wrt_regression_vector(vn, h);
                var_dy_dv = diag(Cov_dy_dv);

                % Save Covariance Matrix in Cell
                Cov_dy_dv_cell{i} = Cov_dy_dv;

                % Update jacobi matrix
                if i > 1
                    P(i, 1:(i-1)) = dy_dv((i-1):-1:1);
                    Var_P(i, 1:(i-1)) = var_dy_dv((i-1):-1:1);
                end
            end
        end

        function [dy_dv, Cov_dy_dv] = approx_gradient_wrt_regression_vector(obj, vn, h)
            %approx_gradient_wrt_regression_vector Compute approximated GP output gradient w.r.t. a regression vector.
            %
            %   Inputs:
            %       vn : Current regression vector (column or row vector)
            %
            %   Outputs:
            %       dy_dv : Gradient of GP output with respect to vn (column vector)

            % General input Matrix containing the regression vector N times
            X = repmat(vn, obj.N, 1);

            % Input Matrix with positive step-size (one input vector unit vector)
            X_plus = X + h*eye(obj.N);

            % Input Matrix with negative step-size (one input vector unit vector)
            X_minus = X - h*eye(obj.N);

            % Perform gp Prediction for positive and negative input (X* =
            % X_plus_minus)
            X_plus_minus = [X_plus; X_minus];
            f = predict(obj.GP, X_plus_minus);

            % Calculation Matrix
            A = [1/(2*h)*eye(obj.N), -1/(2*h)*eye(obj.N)];

            % Compute approximated Gradient
            dy_dv = A*f;

            % Compute kernel Value (X* = X_plus_minus, K** = K(x*, X*))
            K_pm_pm = obj.squared_exponantial_kernel_cross(X_plus_minus, X_plus_minus);

            % Compute kernel Value (K* = K(X*, X_train))
            K_pm_tr = obj.squared_exponantial_kernel_cross(obj.V, X_plus_minus);

            % Efficient computation of the covariance (unsing cholesky
            % factors)
            B = obj.L_chol \ K_pm_tr;           % B = L^-1 * K*
            K_f = K_pm_pm - B.' * B;        % K_f = K** - B^T*B = K** - K*^T [K + sigma^2*I] K*

            % Covariance Matrix
            Cov_dy_dv = A*K_f*A.';
        end

        function K = squared_exponantial_kernel_cross(obj, X1, X2)
            %squared_exponantial_kernel_cross Compute squared exponential (RBF) kernel matrix for set of input vectors.
            %
            %   Inputs:
            %       X1 : First set of input vectors
            %       X2 : Second set of input vector
            %       sigmaL : Length scale parameter (positive scalar)
            %       sigmaF : Signal standard deviation (positive scalar)
            %
            %   Outputs:
            %       K : Kernel Matrix

            % Squared euclidean distances
            D2 = pdist2(X1, X2, 'squaredeuclidean');

            % Squared exponential Kernel
            K = obj.sigmaF2 * exp(-0.5 * D2 * obj.sigmaL2_inv);
        end
    end

    methods (Static)
        function k = squared_exponantial_kernel(x1, x2, sigmaL, sigmaF)
            %squared_exponantial_kernel Compute squared exponential (RBF) kernel value.
            %
            %   Inputs:
            %       x1 : First input vector (column vector)
            %       x2 : Second input vector (column vector)
            %       sigmaL : Length scale parameter (positive scalar)
            %       sigmaF : Signal standard deviation (positive scalar)
            %
            %   Outputs:
            %       k : Kernel value (scalar)

            % Distance
            r2 = sum((x1 - x2).^2);

            % Kernel value
            k = sigmaF^2 * exp(-0.5 * r2 / (sigmaL^2));
        end
    end
end

