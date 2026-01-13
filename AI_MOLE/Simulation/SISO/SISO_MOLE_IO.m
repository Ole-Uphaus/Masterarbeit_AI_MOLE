classdef SISO_MOLE_IO < handle
    %SISO_MOLE_IO Implementation of the AI-MOLE framework for SISO systems
    
    properties
        GP_SISO             % Gaussian Process model for the SISO system
        ILC_SISO            % Iterative Learning Control (ILC) instance
        y_cell              % Cell array storing output trajectories (per iteration)
        u_cell              % Cell array storing input trajectories (per iteration)
        N_iter              % Total number of learning iterations
        H_trials            % Number of previous trials used for GP training (history length)
        m_delay             % System delay (in samples) - relative degree
        i_iter              % Current Iteration counter
        weight_init_method  % Tells, how the weighting Matrices are getting Initialized
        nonlin_damping      % Tells, if nonlin damping should be used for Regularisation
        beta                % Parameter for nonlin damping
        alpha_log           % Vector containing the used damping factors alpha
    end
    
    methods
        function obj = SISO_MOLE_IO(r_vec, u_init, params, sigma_n)
            %SISO_MOLE_IO Constructor for the SISO_MOLE_IO class
            %
            %   Inputs:
            %       r_vec               : Reference trajectory (column vector)
            %       u_init              : Initial input trajectory (column vector)
            %       params              : Struct containing model parameters
            %           -m_delay             : System delay (integer)
            %           -N_iter              : Number of learning iterations
            %           -H_trials            : Number of previous trials used for GP training
            %           -weight_init_method  : Method being used for W, R, S initialisation
            %           -nonlin_damping      : if nonlin damping should be used for Regularisation
            %           -beta                : Parameter for nonlin damping
            %       sigma_n             : Measurement Noise (if known before)
            
            % Generate dynamic paths
            base_dir = fileparts(mfilename("fullpath"));
            GP_path = fullfile(base_dir, '..', '..', '..', 'GP', 'GP_SISO');
            ILC_path = fullfile(base_dir, '..', '..', '..', 'ILC', 'Simulation', 'ILC_SISO');

            % Add temporary paths
            addpath(GP_path);
            addpath(ILC_path);

            % Initialize GP (use fixed sigma_n if given)
            if nargin < 4 || isempty(sigma_n)
                obj.GP_SISO = GP_SISO_IO();
            else
                obj.GP_SISO = GP_SISO_IO(sigma_n);
            end

            % Initialize ILC
            obj.ILC_SISO = ILC_SISO(r_vec, params.m_delay, u_init, params.N_iter);

            % Storage cells
            obj.y_cell = cell(params.N_iter+1, 1);
            obj.u_cell = cell(params.N_iter+1, 1);
            obj.alpha_log = NaN(params.N_iter, 1);

            % Assign values
            obj.N_iter = params.N_iter;
            obj.u_cell{1} = u_init;
            obj.H_trials = params.H_trials;
            obj.m_delay = params.m_delay;
            obj.weight_init_method = params.weight_init_method;
            obj.nonlin_damping = params.nonlin_damping;
            obj.beta = params.beta;
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
            obj.i_iter = sum(~cellfun(@isempty, obj.u_cell));

            % Save Trajectory (measured or simulated)
            obj.y_cell{obj.i_iter} = y_vec;

            % Prepare training data for the GP
            if obj.i_iter <= obj.H_trials
                u_cell_train = obj.u_cell(1:obj.i_iter);
                y_cell_train = obj.y_cell(1:obj.i_iter);
            else
                u_cell_train = obj.u_cell(obj.i_iter-(obj.H_trials-1):obj.i_iter);
                y_cell_train = obj.y_cell(obj.i_iter-(obj.H_trials-1):obj.i_iter);
            end

            % Train GP model
            obj.GP_SISO.train_GP_model(y_cell_train, u_cell_train);

            % Linearize GP model (fast function)
            [P, ~, Cov_dy_du_cell] = obj.GP_SISO.linearize_covariance_at_given_trajectory_fast(obj.u_cell{obj.i_iter});
            % [P, ~, Cov_dy_du_cell, ~, H_dy_du_cell] = obj.GP_SISO.linearize_covariance_hessian_at_given_trajectory_fast(obj.u_cell{obj.i_iter});

            % Change dimensions of P (delete first row and last column -
            % because ILC uses a reduced size framework)
            P_reduced = P(obj.m_delay+1:end, 1:end-obj.m_delay);

            % Calculate weighting matrices (different Methods)
            [W, R, S] = obj.design_weighting_matrices(P_reduced, Cov_dy_du_cell);

            % Initialize norm optimal ILC
            obj.ILC_SISO.init_Quadr_type(W, S, R, P_reduced);

            % Perform ILC update (choose between normal and GP method)
            u_vec_new = obj.ILC_SISO.Quadr_update(y_vec);

            % y_vec_GP = obj.GP_SISO.predict_trajectory_measurement(obj.u_cell{obj.i_iter});
            % u_vec_new = obj.ILC_SISO.Quadr_update_GP_based(y_vec, y_vec_GP);

            % Apply damping based on nonlinearity
            u_vec_new = obj.apply_nonlinearity_damping(u_vec_new, P);

            % Update current input trajectory inside the ILC object, to
            % ensure that the next optimisation step uses the right input
            % trajectory to update
            obj.ILC_SISO.set_current_u_vec(u_vec_new);

            % Save new input
            obj.u_cell{obj.i_iter+1} = [u_vec_new; 0];
        end

        function [W, R, S] = design_weighting_matrices(obj, P, Cov_dy_du_cell)
            %design_weighting_matrices
    
            switch obj.weight_init_method
                case 'Meindl'
                    % Meindls method

                    % Choose W, S, R
                    W = eye(size(P));
                    s = (norm(P, 2)^2);
                    S = s * eye(size(P));
                    r = 0;
                    R = r * eye(size(P));

                    % Print
                    fprintf('Iteration = %d | r = %.4e | s = %.4e \n', obj.i_iter, r, s);

                case 'Stochastic'
                    % Stochastic initialisation using covariance matrices

                    % Choose W, R
                    w = 1;
                    r = 0;
                    W = w * eye(size(P));
                    R = r * eye(size(P));

                    % Calculate S using covariances
                    Cov_size = size(Cov_dy_du_cell{end}, 1);
                    Cov_sum = zeros(Cov_size, Cov_size);

                    for i = 1:numel(Cov_dy_du_cell)
                        % Get current covariance matrix and Hessian
                        Cov_i = Cov_dy_du_cell{i};

                        % Check if covariance is empty (fist matrix should
                        % be always empty)
                        if isempty(Cov_i)
                            continue;
                        end

                        n_i = size(Cov_i, 1);

                        % Add matrix to the sum (upper left corner)
                        Cov_sum(1:n_i, 1:n_i) = Cov_sum(1:n_i, 1:n_i) + Cov_i;
                    end

                    S = w * Cov_sum;

                    % Print
                    fprintf('Iteration = %d | r = %.4e | ||S|| = %.4e \n', obj.i_iter, r, norm(S, 2));

                case 'Heuristic'
                    % Heuristic method (prediction new trajectory)

                    % Heuristic needs to be developed
                    error('Noch keine Heuristik implementiert.')

                case 'Robust'
                    % Ensure robust monotonic convergence

                    % Choose W, S
                    W = eye(size(P));
                    s = obj.GP_SISO.GP.Sigma;
                    S = s * eye(size(P));
        
                    % Calculate maximum error in P
                    Sigma_P = sqrt(max(Var_P, 0));
                    Delta_P_max = 3 * Sigma_P;
                    norm_Delta_P_max = norm(Delta_P_max, 2);
        
                    % Smalles Eigenvalue of A = P' W P + S
                    A = P.' * W * P + S;
                    A = (A + A.')/2;        % Numerically symmetric
                    lambda_min_A = min(eig(A));
        
                    % Numerator of the inequation
                    norm_P = norm(P, 2);
                    num = norm(S, 2) + norm_P * norm(W, 2) * norm_Delta_P_max;
        
                    % Calculate lower bound of r
                    r_lower_bound = num - lambda_min_A;
        
                    % Check if r is positive
                    r = max(0, 1.01 * r_lower_bound);
        
                    % Weighting Matrix
                    R = r * eye(size(P));

                    % Print
                    fprintf('Iteration = %d | r = %.4e | s = %.4e | ||P|| = %.4e | ||delta_P|| = %.4e \n', obj.i_iter, r, s, norm_P, norm_Delta_P_max);

                case 'Manual'
                    % Manual parameterisation

                    % Choose W, S, R
                    w = 1;
                    s = 0.01;
                    r = 0.01;
                    W = w * eye(size(P));
                    S = s * eye(size(P));
                    R = r * eye(size(P));

                    % Print
                    fprintf('Iteration = %d | r = %.4e | s = %.4e \n', obj.i_iter, r, s);

                otherwise
                    error('Unbekannte Initialisierungsmethode ausgewählt.')
            end
        end

        function u_vec_new_damped = apply_nonlinearity_damping(obj, u_vec_new, P)

            % Compute delta_u
            delta_u = [u_vec_new; 0] - obj.u_cell{obj.i_iter};

            switch obj.nonlin_damping
                case 'none'
                    % No damping applid so alpha is set to one
                    alpha = 1;
                
                case 'relative_1'
                    % Use the relative error between linearized Gp and
                    % actual GP to parametrize alpha
        
                    % Prediction using GP
                    [y_pred_GP, Cov_y] = obj.GP_SISO.predict_trajectory_covariance([u_vec_new; 0]);
                    y_std_GP = sqrt(diag(Cov_y));
        
                    % Prediction using Linearisation
                    y_pred_lin = obj.GP_SISO.predict_trajectory_measurement(obj.u_cell{obj.i_iter}) + P*delta_u;
        
                    % Relative prediction difference (use variances too)
                    eta = norm((abs(y_pred_GP - y_pred_lin) + obj.beta*y_std_GP), 2) / (norm((P*delta_u), 2) + eps);
        
                    % Damping factor
                    alpha = 1 / (1 + eta);
                    % alpha = max(1-eta, 0.1);

                    % Print parameters
                    fprintf('    eta = %.3f | alpha = %.3f \n', eta, alpha);
                
                case 'relative_2'
                    % Use the relative error between linearized Gp and
                    % actual GP to parametrize alpha
        
                    % Prediction using GP
                    [y_pred_GP, Cov_y] = obj.GP_SISO.predict_trajectory_covariance([u_vec_new; 0]);
                    y_std_GP = sqrt(diag(Cov_y));
        
                    % Prediction using Linearisation
                    y_pred_lin = obj.GP_SISO.predict_trajectory_measurement(obj.u_cell{obj.i_iter}) + P*delta_u;
        
                    % Relative prediction difference (use variances too)
                    eta = norm((abs(y_pred_GP - y_pred_lin) + obj.beta*y_std_GP), 2) / (norm((abs(P*delta_u) + abs(y_pred_GP - y_pred_lin) + obj.beta*y_std_GP), 2) + eps);
        
                    % Damping factor
                    alpha = 1 - eta;

                    % Print parameters
                    fprintf('    eta = %.3f | alpha = %.3f \n', eta, alpha);

                case 'minimize'
                    % Use the GP to Calculate the expected error (with
                    % uncertainty) for different alpha values (grid search)

                    % Set sigma constraint (3 for 99%)
                    k = 3;

                    % Alpha values
                    alpha_vec = linspace(0.1, 1, 10);
                    alpha_vec = alpha_vec(:);

                    % Calulate current iteration squared error norm
                    ej = obj.ILC_SISO.r_vec - obj.y_cell{obj.i_iter};
                    ej_sq_norm = sum(ej.^2);

                    % Expected value of the quadratic error and Variance
                    expected_sq_error_vec = zeros(length(alpha_vec), 1);
                    expected_sq_error_variance_vec = zeros(length(alpha_vec), 1);
                    feasible = false(length(alpha_vec), 1);

                    for i = 1:length(alpha_vec)
                        % GP prediction depending of alpha
                        u_vec_temp = obj.u_cell{obj.i_iter} + alpha_vec(i)*delta_u;
                        [y_pred_temp, Cov_y_temp] = obj.GP_SISO.predict_trajectory_covariance(u_vec_temp);

                        % Calculate Error mean
                        mu_e_jp1 = obj.ILC_SISO.r_vec - y_pred_temp;

                        % Calculate expected squared error
                        expected_sq_error_vec(i) = (mu_e_jp1.' * mu_e_jp1) + trace(Cov_y_temp);

                        % Calculate variance of squared error
                        expected_sq_error_variance_vec(i) = 2*trace(Cov_y_temp * Cov_y_temp) + 4*(mu_e_jp1.' * Cov_y_temp * mu_e_jp1);

                        % Proof feasibility
                        feasible(i) = ((expected_sq_error_vec(i) + k*sqrt(expected_sq_error_variance_vec(i))) < ej_sq_norm);
                    end

                    % Find optimal alpha that is feasible
                    if any(feasible)
                        % Minimize expected error among feasible alphas
                        [~, idx_min] = min(expected_sq_error_vec(feasible));
                        alpha_candidates = alpha_vec(feasible);
                        alpha = alpha_candidates(idx_min);
                        
                        if any(~feasible)
                            status = "feasible with constraints";
                        else
                            status = "feasible for every alpha";
                        end
                    else
                        % Fallback if no alpha is feasible (choose smallest
                        % alpha for update)
                        alpha = min(alpha_vec);

                        status = "infeasible";
                    end

                    % Print parameters
                    fprintf('    alpha = %.3f (%s)\n', alpha, status);

                otherwise
                    error('Unbekannte nichtlineare Dämpfungsmethode ausgewählt.')
            end

            % Calculate damped input Trajectory
            u_vec_new_damped = obj.u_cell{obj.i_iter} + alpha*delta_u;

            % Log alpha
            obj.alpha_log(obj.i_iter) = alpha;

            % Delete last element to ensure reduced size representation
            u_vec_new_damped = u_vec_new_damped(1:end-1);
        end

        function save_final_trajectory(obj, y_vec)
            %save_final_trajectory Save the final output trajectory and log final error
            %
            %   Inputs:
            %       y_vec : Measured or simulated output of the final iteration

            % Get current iteration counter
            obj.i_iter = sum(~cellfun(@isempty, obj.u_cell));

            % Save Trajectory (measured or simulated)
            obj.y_cell{obj.i_iter} = y_vec;

            % Calculate final Error ILC
            obj.ILC_SISO.calculate_final_error(y_vec);
        end
    end
end

