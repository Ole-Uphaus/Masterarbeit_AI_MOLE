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
        sigmaF2     % Squared kernel parameter sigma_F^2
        sigmaL2_inv % Squared inverse kernel parameter 1/l^2
        alpha       % Constant term of the GP a = [K + sigma^2*I]^-1 * y
        V_row_norm2 % Euklidean Norm of each row of obj.V
        V_transp    % Transposed design Matrix

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

            N = length(obj.u_cell{1});
            J = length(obj.u_cell);     % Nuber of training trials

            % Empty design matrix
            obj.V = zeros(J*N, N);

            % Design Matrix
            for i = 1:J
                % Toeplitz Matrix
                V_temp = toeplitz([0; obj.u_cell{i}(1:end-1)], zeros(1, N));

                % Put Toeplitz into Design Matrix
                obj.V(((i-1)*N + 1):(i*N), :) = V_temp;
            end
    
        end

        function constr_target_vector(obj)
            %constr_target_vector Construct target vector for GP regression.

            N = length(obj.y_cell{1});
            J = length(obj.y_cell);     % Nuber of training trials

            % Empty Target vector
            obj.z = zeros(J*N, 1);

            % Target Vector
            for i = 1:J
                obj.z(((i-1)*N + 1):(i*N), :) = obj.y_cell{i};
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

            N = length(obj.u_cell{1});

            % Toeplitz test Matrix
            V_test = toeplitz([0; u_test(1:end-1)], zeros(1, N));
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

            N = length(obj.u_cell{1});

            % Empty jacobi matrix
            P = zeros(N, N);

            % Construct linearisation Matrix
            V_lin = obj.constr_test_matrix(u_lin);

            % Iteratively compute gradients
            for i = 1:N
                % Extract regression Vector
                vn = V_lin(i, :);

                % Gradient of the GP w.r.t vn
                dy_dv = obj.gradient_wrt_regression_vector(vn);

                % Update jacobi matrix
                J = i - 1;
                for j = 1:J
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

        function P = linearize_at_given_trajectory_fast(obj, u_lin)
            %linearize_at_given_trajectory Compute local linearization (Jacobian) of GP at a trajectory.
            %
            %   Inputs:
            %       u_lin : Input trajectory (column vector) around which the GP is linearized
            %
            %   Outputs:
            %       P : Jacobian (linearization) matrix of size N×N, mapping input variations
            %           Δu to output variations Δy near the trajectory u_lin

            N = length(obj.u_cell{1});

            % Empty jacobi matrix
            P = zeros(N, N);

            % Construct linearisation Matrix
            V_lin = obj.constr_test_matrix(u_lin);

            % Kernel Parameters
            k_params = obj.GP.KernelInformation.KernelParameters;
            obj.sigmaL2_inv = 1/(k_params(1)^2);
            obj.sigmaF2 = k_params(2)^2;

            % Receive alpha vector a = [K + sigma^2*I]^-1 * y
            obj.alpha = obj.GP.Alpha;

            % Precompute euklidean Norm of each row of obj.V and transposed
            % design matrix (each column is one regression Vector)
            obj.V_row_norm2 = sum(obj.V.^2, 2);
            obj.V_transp = obj.V.'; % .' - no complex transposition

            % Iteratively compute gradients
            for i = 1:N
                % Extract regression Vector
                vn = V_lin(i, :);

                % Gradient of the GP w.r.t vn
                dy_dv = obj.gradient_wrt_regression_vector_fast(vn);

                % Update jacobi matrix
                if i > 1
                    P(i, 1:(i-1)) = dy_dv((i-1):-1:1);
                end
            end
        end

        function dy_dv = gradient_wrt_regression_vector_fast(obj, vn)
            %gradient_wrt_regression_vector Compute GP output gradient w.r.t. a regression vector.
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

        function P = approx_linearisation_at_given_trajectory(obj, u_lin)
            %approx_linearisation_at_given_trajectory Compute local approximated linearization (Jacobian) of GP at a trajectory.
            %
            %   Inputs:
            %       u_lin : Input trajectory (column vector) around which the GP is linearized
            %
            %   Outputs:
            %       P : Jacobian (linearization) matrix of size N×N, mapping input variations
            %           Δu to output variations Δy near the trajectory u_lin

            N = length(obj.u_cell{1});

            % Empty jacobi matrix
            P = zeros(N, N);

            % Construct linearisation Matrix
            V_lin = obj.constr_test_matrix(u_lin);

            % Step-size
            h = 1e-5;

            % Iteratively compute gradients
            for i = 1:N
                % Extract regression Vector
                vn = V_lin(i, :);

                % Gradient of the GP w.r.t vn
                dy_dv = obj.approx_gradient_wrt_regression_vector(vn, h, N);

                % Update jacobi matrix
                if i > 1
                    P(i, 1:(i-1)) = dy_dv((i-1):-1:1);
                end
            end
        end

        function dy_dv = approx_gradient_wrt_regression_vector(obj, vn, h, N)
            %approx_gradient_wrt_regression_vector Compute approximated GP output gradient w.r.t. a regression vector.
            %
            %   Inputs:
            %       vn : Current regression vector (column or row vector)
            %
            %   Outputs:
            %       dy_dv : Gradient of GP output with respect to vn (column vector)

            % General input Matrix containing the regression vector N times
            X = repmat(vn, N, 1);

            % Input Matrix with positive step-size (one input vector unit vector)
            X_plus = X + h*eye(N);

            % Input Matrix with negative step-size (one input vector unit vector)
            X_minus = X - h*eye(N);

            % Perform gp Prediction for positive and negative input
            [f_plus, ~] = predict(obj.GP, X_plus);
            [f_minus, ~] = predict(obj.GP, X_minus);

            % Concatenate Vectors
            f = [f_plus; f_minus];

            % Calculation Matrix
            A = [1/(2*h)*eye(N), -1/(2*h)*eye(N)];

            % Compute approximated Gradient
            dy_dv = A*f;
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

