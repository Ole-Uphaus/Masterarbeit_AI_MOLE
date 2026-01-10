function P = Lifted_dynamics_linear_SISO(Ad, Bd, Cd, N, m)
%LIFTED_DYNAMICS_SISO Construct lifted system matrix for a linear SISO system
%   Inputs:
%       Ad  - Discrete-time state matrix  (n_x × n_x)
%       Bd  - Discrete-time input matrix  (n_x × 1)
%       Cd  - Discrete-time output matrix (1 × n_x)
%       N   - Time horizon (total number of samples)
%       m   - System delay (number of steps before first output appears)
%
%   Outputs:
%       P   - Lifted system matrix ( (N-m) × (N-m) )
%             Lower-triangular Toeplitz matrix constructed from Markov parameters

nx = size(Ad, 1);

% Markov-Parameter p (respect delay m)
p = zeros(N-m, 1);
A_pow = eye(nx);

for k = 1:(N-m)
    p(k) = Cd * A_pow * Bd;
    A_pow = A_pow * Ad;
end

% Toeplitz-Matrix P
P = toeplitz(p, [p(1), zeros(1, N-1-m)]); % define first column and first row
end

