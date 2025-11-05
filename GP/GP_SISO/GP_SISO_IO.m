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
    end
    
    methods
        function obj = GP_SISO_IO()
            %GP_SISO_IO Constructor.            
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
            obj.GP = fitrgp( ...
                obj.V, obj.z, ...
                'BasisFunction','none', ...
                'KernelFunction','squaredexponential', ...
                'Standardize', false);
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
            alpha = obj.GP.Alpha;

            N_alpha = length(alpha);
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
            dy_dv = dk_dvn * alpha;
        end

        function k = squared_exponantial_kernel(obj, x1, x2, sigmaL, sigmaF)
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

