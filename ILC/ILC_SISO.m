classdef ILC_SISO < handle
    %ILC_LINEAR_SISO Iterative Learning Control vor linear SISO Systems.
    %
    % The class implements different Methods of ILC (PD-Type and quadratic
    % optimal)
    
    properties
        r_vec       % Reference Trajectory
        m           % System delay
        u_vec       % System input
        kp          % P Gain
        kd          % D Gain
        error_vec   % Tracking error
        RMSE_log    % Error log
        P           % Lifted Dynamics
        W           % Wheighting Matrix error
        S           % Wheighting Matrix change of u
        L           % Learning Matrix
    end
    
    methods
        function obj = ILC_SISO(r_vec, m)
            %ILC_LINEAR_SISO Construct an instance of this class
            %
            %   Inputs:
            %       r_vec - Reference signal (column vector)
            %       m     - System delay (integer)

            obj.r_vec = r_vec;
            obj.m = m;
            % Initiallize input vector u
            obj.u_vec = zeros(size(r_vec, 1)-m, 1);
        end
        
        function init_PD_type(obj, kp, kd)
            %init_PD_type Initialize parameters for PD-type ILC
            %
            %   Inputs:
            %       kp - Proportional gain
            %       kd - Derivative gain

            obj.kp = kp;
            obj.kd = kd;
        end

        function u_vec_new = PD_update(obj, y_vec)
            %PD_update Perform one PD-type ILC iteration
            %
            %   Inputs:
            %       y_vec - Measured output signal of the current trial
            %
            %   Outputs:
            %       u_vec_new - Updated input signal for the next trial

            % Calculate error e(k+1), e(k)
            obj.error_vec = obj.r_vec((1+obj.m):end) - y_vec((1+obj.m):end);
            d_error_vec = [obj.error_vec(1); diff(obj.error_vec)];

            % Log error
            obj.Log_error()

            % Update input u(k)
            u_vec_new = obj.u_vec + obj.kp*obj.error_vec + obj.kd*d_error_vec;
            obj.u_vec = u_vec_new;
        end

        function init_Quadr_type(obj, W, S, P)
            %init_Quadr_type Initialize quadratic optimal ILC parameters
            %
            %   Inputs:
            %       W - Weighting matrix error
            %       S - Weighting matrix change of u
            %       P - Lifted system (Plant) matrix

            % Set W and S
            if nargin >= 2 && ~isempty(W), obj.W = W; end
            if nargin >= 3 && ~isempty(S), obj.S = S; end
        
            % Only set P if given
            if nargin >= 4 && ~isempty(P)
                obj.P = P;
            end

            % Optimal learning function (only if everything is given)
            if ~isempty(obj.P) && ~isempty(obj.W) && ~isempty(obj.S)
                obj.L = (P'*W*P + S) \ (P'*W);
            else
                obj.L = [];
            end
        end

        function u_vec_new = Quadr_update(obj, y_vec, P)
            %Quadr_update Perform one Quadratic Optimal ILC iteration
            %
            %   Inputs:
            %       y_vec - Measured output signal of the current trial
            %       P     - Lifted system (Plant) matrix
            %
            %   Outputs:
            %       u_vec_new - Updated input signal for the next trial

            % If P is not defined
            if isempty(obj.P)
                % Calculate L (current L -> Lc) using the given P
                if nargin >= 3 && ~isempty(P)
                    Lc = (P'*obj.W*P + obj.S) \ (P'*obj.W);
                end
            else
                Lc = obj.L;
            end

            % Calculate error e(k+1)
            obj.error_vec = obj.r_vec((1+obj.m):end) - y_vec((1+obj.m):end);

            % Log error
            obj.Log_error()

            % Update unput u(k)
            u_vec_new = obj.u_vec + Lc*obj.error_vec;
            obj.u_vec = u_vec_new;
        end

        function Log_error(obj)
            %init_Quadr_type Logs RMSE over iterations
            
            % Calculate RMSE
            RMSE = sqrt(mean(obj.error_vec.^2));
            obj.RMSE_log(end+1) = RMSE;
        end 
    end
end

