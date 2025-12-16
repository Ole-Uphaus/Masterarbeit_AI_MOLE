classdef ILC_SISO < handle
    %ILC_LINEAR_SISO Iterative Learning Control vor linear and nonlinear SISO Systems.
    %
    % The class implements different Methods of ILC (PD-Type and quadratic
    % optimal)
    
    properties
        r_vec       % Reference Trajectory
        m           % System delay
        u_vec       % System input
        error_vec   % Tracking error
        RMSE_log    % Error log
        P           % Lifted Dynamics
        
        W           % Wheighting Matrix error
        S           % Wheighting Matrix change of u
        R           % Wheighting Matrix input u

        L           % Optimal Learning Matrix
        Q           % Optimal Q-Operator

        a_Q         % Filter Parameter a
        b_Q         % Filter Parameter b
    end
    
    methods
        function obj = ILC_SISO(r_vec, m, u_init)
            %ILC_LINEAR_SISO Construct an instance of this class
            %
            %   Inputs:
            %       r_vec - Reference signal (column vector)
            %       m     - System delay (integer)

            obj.r_vec = r_vec;
            obj.m = m;

            % Initialize input vector u
            obj.u_vec = u_init(1:(size(r_vec, 1)-m), 1);
        end

        function init_Quadr_type(obj, W, S, R, P)
            %init_Quadr_type Initialize quadratic optimal ILC parameters
            %
            %   Inputs:
            %       W - Weighting matrix error
            %       S - Weighting matrix change of u
            %       R - Wheighting Matrix input u
            %       P - Lifted system (Plant) matrix

            % Set W, S, R, P
            if nargin >= 2 && ~isempty(W), obj.W = W; end
            if nargin >= 3 && ~isempty(S), obj.S = S; end
            if nargin >= 4 && ~isempty(R), obj.R = R; end
            if nargin >= 5 && ~isempty(P), obj.P = P; end

            % Optimal learning function (only if everything is given)
            if ~isempty(obj.P) && ~isempty(obj.W) && ~isempty(obj.S) && ~isempty(obj.R)
                % Use the notation from norm optimal ILC Papers (apply the
                % Q-Operator on u and the L-Operator on e - seperate)

                % Calculation Matrix
                M = obj.P'*obj.W*obj.P + obj.R + obj.S;

                % Optimal Q-Operator
                obj.Q = M \ (obj.P'*obj.W*obj.P + obj.S);

                % Optimal L-Operator
                obj.L = M \ (obj.P'*obj.W);
            else
                obj.L = [];
                obj.Q = [];
            end
        end

        function u_vec_new = Quadr_update(obj, y_vec)
            %Quadr_update Perform one Quadratic Optimal ILC iteration
            %
            %   Inputs:
            %       y_vec - Measured output signal of the current trial
            %
            %   Outputs:
            %       u_vec_new - Updated input signal for the next trial

            % Calculate error e(k+1)
            obj.error_vec = obj.r_vec((1+obj.m):end) - y_vec((1+obj.m):end);

            % Log error
            obj.Log_error()

            % Update input u(k) and apply filter
            u_vec_new = obj.apply_Q_lowpass(obj.Q*obj.u_vec + obj.L*obj.error_vec);
            obj.u_vec = u_vec_new;
        end

        function Log_error(obj)
            %Log_error Compute and store RMSE of the current tracking error.
            
            % Calculate RMSE
            RMSE = sqrt(mean(obj.error_vec.^2));
            obj.RMSE_log(end+1) = RMSE;
        end

        function calculate_final_error(obj, y_vec)
            %calculate_final_error Calculate and log the last error during
            %Training
            %
            %   Inputs:
            %       y_vec - Measured output signal of the current trial

            % Calculate error e(k+1)
            obj.error_vec = obj.r_vec((1+obj.m):end) - y_vec((1+obj.m):end);

            % Log error
            obj.Log_error()
        end

        function init_Q_lowpass(obj, fc, order, Ts)
            %init_Q_lowpass Initialize an optional zero-phase low-pass Q-filter.
            %
            %   Inputs:
            %       fc    : Cutoff frequency
            %       order : Butterworth filter order
            %       Ts    : Sampling time in seconds

            fs = 1/Ts;

            % Butterworth lowpass
            [b,a] = butter(order, fc/(fs/2));

            % Set filter parameters
            obj.a_Q = a;
            obj.b_Q = b;
        end

        function x_filt = apply_Q_lowpass(obj, x)
            %apply_Q_lowpass Apply the optional Q-filter to a signal.
            %
            %   Inputs:
            %       x : Column vector to be filtered (e.g., input update)
            %
            %   Outputs:
            %       x_filt : Filtered signal if Q is initialized; otherwise x unchanged
            
            % Apply Q filter if initialized
            if ~isempty(obj.a_Q) && ~isempty(obj.b_Q)
                x_filt = filtfilt(obj.b_Q, obj.a_Q, x);
            else
                x_filt = x;
            end
        end

        function set_current_u_vec(obj, u_vec_new)
            %set_current_u_vec Manually set the current input vector. This
            %is used in AI-MOLE to apply a damping factor to the optimized
            %input trajectory.
            %
            %   Inputs:
            %       u_vec_new : new input vector (damped)

            obj.u_vec = u_vec_new;
        end
    end
end

