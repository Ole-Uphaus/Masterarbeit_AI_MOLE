classdef ILC_linear_SISO < handle
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
        function obj = ILC_linear_SISO(r_vec, m)
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

            % Update unput u(k)
            u_vec_new = obj.u_vec + obj.kp*obj.error_vec + obj.kd*d_error_vec;
            obj.u_vec = u_vec_new;
        end

        function init_Quadr_type(obj, P, W, S)
            %init_Quadr_type Initialize quadratic optimal ILC parameters
            %
            %   Inputs:
            %       P - Lifted system (Plant) matrix
            %       W - Weighting matrix error
            %       S - Weighting matrix change of u

            obj.P = P;
            obj.W = W;
            obj.S = S;

            % Optimal learning function
            obj.L = inv(P'*W*P + S)*P'*W;
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

            % Update unput u(k)
            u_vec_new = obj.u_vec + obj.L*obj.error_vec;
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

