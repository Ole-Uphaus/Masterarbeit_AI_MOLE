function P = Lifted_dynamics_nonlinear_SISO(linear_discrete_system, N, m, x_sim)
%LIFTED_DYNAMICS_NONLINEAR_SISO Summary of this function goes here
%   Detailed explanation goes here

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
        P(k, i) = Cd_seq{k+1} * A_pow * Bd_seq{i};
        A_pow = A_pow * Ad_seq{i};
    end
end

end

