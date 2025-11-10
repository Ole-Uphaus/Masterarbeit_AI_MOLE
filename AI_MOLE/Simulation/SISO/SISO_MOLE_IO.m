classdef SISO_MOLE_IO < handle
    %SISO_MOLE_IO Implementation of the AI-MOLE framework for SISO systems
    
    properties
        GP_SISO         % Gaussian Process model for the SISO system
        ILC_SISO        % Iterative Learning Control (ILC) instance
        y_cell          % Cell array storing output trajectories (per iteration)
        u_cell          % Cell array storing input trajectories (per iteration)
        N_iter          % Total number of learning iterations
        H_trials        % Number of previous trials used for GP training (history length)
        m_delay         % System delay (in samples) - relative degree
    end
    
    methods
        function obj = SISO_MOLE_IO(r_vec, m_delay, u_init, N_iter, H_trials)
            %SISO_MOLE_IO Constructor for the SISO_MOLE_IO class
            %
            %   Inputs:
            %       r_vec    : Reference trajectory (column vector)
            %       m_delay  : System delay (integer)
            %       u_init   : Initial input trajectory (column vector)
            %       N_iter   : Number of learning iterations
            %       H_trials : Number of previous trials used for GP training
            
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
            %update_input Perform one AI-MOLE iteration (GP + ILC update)
            %
            %   Inputs:
            %       y_vec : Measured or simulated output signal of the current iteration
            %
            %   Outputs:
            %       u_vec_new : Updated input signal for the next iteration

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
            %save_final_trajectory Save the final output trajectory and log final error
            %
            %   Inputs:
            %       y_vec : Measured or simulated output of the final iteration

            % Get current iteration counter
            i_iter = sum(~cellfun(@isempty, obj.u_cell));

            % Save Trajectory (measured or simulated)
            obj.y_cell{i_iter} = y_vec;

            % Calculate final Error ILC
            obj.ILC_SISO.calculate_final_error(y_vec);
        end
    end
end

