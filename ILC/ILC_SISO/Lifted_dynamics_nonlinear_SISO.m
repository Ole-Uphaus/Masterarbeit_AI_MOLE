function P = Lifted_dynamics_nonlinear_SISO(linear_discrete_system, N, m, x_sim)
%LIFTED_DYNAMICS_NONLINEAR_SISO Construct lifted system Matrix for
%nonlinear SISO System. The system is linearized along a given Trajectory.
%   Detailed explanation goes here
%
%   Inputs:
%       linear_discrete_system - Function handle that returns the
%                               linearized discrete-time system matrices
%                               [A_k, B_k, C_k, D_k] at a given state x_k.
%                               Signature: [A, B, C, D] = linear_discrete_system(x)
%
%       N   - Time horizon (total number of samples)
%       m   - System delay (number of steps before first output appears)
%       x_sim - Simulated state trajectory (N × n_x)
%
%   Output:
%       P   - Lifted system matrix ((N-m) × (N-m))
%             Time-varying, lower-triangular matrix constructed from
%             local Markov parameters along the nonlinear trajectory.

% Precompute time varying matrices A_k, B_k, C_k, D_k
Ad_seq = cell(N, 1);
Bd_seq = cell(N, 1);
Cd_seq = cell(N, 1);
Dd_seq = cell(N, 1);
for k = 1:N
    [Ad_seq{k}, Bd_seq{k}, Cd_seq{k}, Dd_seq{k}] = linear_discrete_system(x_sim(k, :));
end

% Compute lifted Matrix
nx = size(Ad_seq{1}, 1);
P = zeros(N-m, N-m);
for k = 1:(N-m)
    % Terms from past inputs
    A_pow = eye(nx);
    for i = (k):-1:1
        % Backwards
        P(k, i) = Cd_seq{k+m} * A_pow * Bd_seq{i};
        A_pow = A_pow * Ad_seq{i};
    end
end
end

