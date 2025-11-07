classdef SISO_MOLE_IO < handle
    %SISO_MOLE_IO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        GP_SISO
        ILC_SISO
        y_cell
        u_cell
        N_iter
        H_trials
        m_delay
    end
    
    methods
        function obj = SISO_MOLE_IO(r_vec, m_delay, u_init, N_iter, H_trials)
            %SISO_MOLE_IO Construct an instance of this class
            %   Detailed explanation goes here
            
            % Generate dynamic paths
            base_dir = fileparts(mfilename("fullpath"));
            GP_path = fullfile(base_dir, '..', '..', '..', 'GP', 'GP_SISO');
            ILC_path = fullfile(base_dir, '..', '..', '..', 'ILC', 'ILC_SISO');

            % Add temporary paths
            addpath(GP_path);
            addpath(ILC_path);

            % Initialize Classes
            obj.GP_SISO = GP_SISO_IO();
            obj.ILC_SISO = ILC_SISO(r_vec, m_delay, u_init);

            % Storage cells
            obj.y_cell = cell(N_iter+1, 1);
            obj.u_cell = cell(N_iter+1, 1);

            % Assign values
            obj.N_iter = N_iter;
            obj.u_cell{1} = u_init;
            obj.H_trials = H_trials;
            obj.m_delay = m_delay;
        end

        function u_vec_new = update_input(obj, y_vec)
            %SISO_MOLE_IO Construct an instance of this class
            %   Detailed explanation goes here

            % Get current iteration counter
            i_iter = sum(~cellfun(@isempty, obj.u_cell));

            % Save Trajectory (measured or simulated)
            obj.y_cell{i_iter} = y_vec;

            % Prepare training data for the GP
            if i_iter <= obj.H_trials
                u_cell_train = obj.u_cell(1:i_iter);
                y_cell_train = obj.y_cell(1:i_iter);
            else
                u_cell_train = obj.u_cell(i_iter-(obj.H_trials-1):i_iter);
                y_cell_train = obj.y_cell(i_iter-(obj.H_trials-1):i_iter);
            end

            % Train GP model
            obj.GP_SISO.train_GP_model(y_cell_train, u_cell_train);

            % Linearize GP model (fast function)
            P = obj.GP_SISO.linearize_at_given_trajectory_fast(obj.u_cell{i_iter});

            % Change dimensions of P (delete first row and last column -
            % because ILC uses a reduced size framework)
            P = P(obj.m_delay+1:end, 1:end-obj.m_delay);

            % Calculate weighting matrices
            W = eye(size(P));
            S = 0.0001 * (norm(P, 2)^2) * eye(size(P));

            % Perform ILC update
            obj.ILC_SISO.init_Quadr_type(W, S, P);
            u_vec_new = obj.ILC_SISO.Quadr_update(y_vec);

            % Save new input
            obj.u_cell{i_iter+1} = [u_vec_new; 0];
        end

        function save_final_trajectory(obj, y_vec)
            %SISO_MOLE_IO Construct an instance of this class
            %   Detailed explanation goes here

            % Get current iteration counter
            i_iter = sum(~cellfun(@isempty, obj.u_cell));

            % Save Trajectory (measured or simulated)
            obj.y_cell{i_iter} = y_vec;

            % Calculate final Error ILC
            obj.ILC_SISO.calculate_final_error(y_vec);
        end
    end
end

