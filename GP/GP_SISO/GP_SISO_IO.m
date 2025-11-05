classdef GP_SISO_IO < handle
    %GP_SISO_IO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        GP          % Gaussian Process object
        y_cell      % Cell containing training output Data (from each trial)
        u_cell      % Cell containing training input Data (from each trial)
        V           % Design matrix (containing u regression vectors)
        z           % Target Vector (containing y target Values)
    end
    
    methods
        function obj = GP_SISO_IO()
            %GP_SISO_IO Construct an instance of this class
            %   Detailed explanation goes here
            
        end
        
        function train_GP_model(obj, y_cell, u_cell)
            %GP_SISO_IO Construct an instance of this class
            %   Detailed explanation goes here

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

        function [y_pred, y_std] = predict_trajectory(obj, u_test)
            %GP_SISO_IO Construct an instance of this class
            %   Detailed explanation goes here

            u_test = u_test(:);  % Ensure column vector
            
            % Construct test Matrix
            V_test = obj.constr_test_matrix(u_test);

            % Prediction
            [y_pred, y_std] = predict(obj.GP, V_test);
        end

        function P = linearize_at_given_trajectory(obj, u_lin)
            %GP_SISO_IO Construct an instance of this class
            %   Detailed explanation goes here

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

        function constr_design_matrix(obj)
            %GP_SISO_IO Construct an instance of this class
            %   Detailed explanation goes here

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
            %GP_SISO_IO Construct an instance of this class
            %   Detailed explanation goes here

            N = length(obj.y_cell{1});
            J = length(obj.y_cell);     % Nuber of training trials

            % Empty Target vector
            obj.z = zeros(J*N, 1);

            % Target Vector
            for i = 1:J
                obj.z(((i-1)*N + 1):(i*N), :) = obj.y_cell{i};
            end
        end

        function V_test = constr_test_matrix(obj, u_test)
            %GP_SISO_IO Construct an instance of this class
            %   Detailed explanation goes here

            N = length(obj.u_cell{1});

            % Toeplitz test Matrix
            V_test = toeplitz([0; u_test(1:end-1)], zeros(1, N));
        end

        function dy_dv = gradient_wrt_regression_vector(obj, vn)
            %GP_SISO_IO Construct an instance of this class
            %   Detailed explanation goes here

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
            %GP_SISO_IO Construct an instance of this class
            %   Detailed explanation goes here

            % Distance
            r2 = sum((x1 - x2).^2);

            % Kernel value
            k = sigmaF^2 * exp(-0.5 * r2 / (sigmaL^2));
        end
    end
end

