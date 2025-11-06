function [q, qd, qdd] = Random_C2_trajectory_1D(num_points, t_vec, sigma)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Parameters
t0 = t_vec(1);
t1 = t_vec(end);

% Time-points equally distributed (except start and end)
time_pts = linspace(t0, t1, num_points + 2);

% Random waypoints (normal distributed)
way_pts = sigma * randn(1, num_points + 2);

% Set initial and final condition tu zero
way_pts(1) = 0;
way_pts(end) = 0;

% Generate trajectory
[q, qd, qdd] = quinticpolytraj(way_pts, time_pts, t_vec);
q = q.';
qd = qd.';
qdd = qdd.';
end

