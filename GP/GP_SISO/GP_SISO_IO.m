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
                'KernelFunction','squaredexponential');

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
                V_temp = toeplitz([0; obj.u_cell{1}(1:end-1)], zeros(1, N));

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
                obj.z(((i-1)*N + 1):(i*N), :) = obj.y_cell{1};
            end
        end

        function V_test = constr_test_matrix(obj, u_test)
            %GP_SISO_IO Construct an instance of this class
            %   Detailed explanation goes here

            N = length(obj.u_cell{1});

            % Toeplitz test Matrix
            V_test = toeplitz([0; u_test(1:end-1)], zeros(1, N));
        end

    end
end

