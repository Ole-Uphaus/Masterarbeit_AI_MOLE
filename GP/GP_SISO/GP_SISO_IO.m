classdef GP_SISO_IO < handle
    %GP_SISO_IO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        GP          % Gaussian Process object
        y_cell      % Cell containing training output Data (from each trial)
        u_cell      % Cell containing training input Data (from each trial)
        V           % Design matrix (containing u regression vectors)
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

            % Design Matrix
            obj.calc_design_matrix();


        end

        function calc_design_matrix(obj)
            %GP_SISO_IO Construct an instance of this class
            %   Detailed explanation goes here

            N = length(obj.u_cell{1});

            % Toeplitz design Matrix
            obj.V = toeplitz(zeros(N, 1), [0, obj.u_cell{1}(1:end-1)']);
            
        end
    end
end

