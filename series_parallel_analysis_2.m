% Rulkov Series vs Parallel Connectivity Analysis (Extended)
% ---------------------------------------------------------------
% Compares burst propagation in extended Rulkov networks:
% 1. Series: A → B → C → D
% 2. Parallel-merged: A → B, A → C, B → D, C → D

clear; clc; close all

%% Simulation Parameters
T = 5000;
alpha = 6.0;
mu = 0.001;
sigma = 0.0;
beta = 0.0;
threshold = -2.0;
gc = -1;  % Fixed inhibitory strength
f = @(x, y) (x <= 0) .* (alpha ./ (1 - x) + y) + ...
           (x > 0 & x < alpha + y) .* (alpha + y) + ...
           (x >= alpha + y) .* (-1);
H = @(v) double(v > 0);

%% Define Network Topologies
% Extended Series: A → B → C → D
A_series = [0 1 0 0;
            0 0 1 0;
            0 0 0 1;
            0 0 0 0];

% Parallel-merged: A → B, A → C, B → D, C → D
A_parallel = [0 1 1 0;
              0 0 0 1;
              0 0 0 1;
              0 0 0 0];

networks = {A_series, A_parallel};
labels = {'Series: A → B → C → D', 'Parallel-Merged: A → {B,C} → D'};

%% Simulate both networks
for net_idx = 1:2
    A = networks{net_idx};
    N = size(A, 1);

    x = -1.5 + 0.01*rand(N, T);
    y = -2 + 0.01*rand(N, T);

    % Stimulate neuron A (neuron 1)
    x(:,1) = -1.5;  % Kick it into bursting

    for t = 1:T-1
        I_inhib = zeros(N, 1);
        for i = 1:N
            for j = 1:N
                if A(i,j) == 1
                    I_inhib(i) = I_inhib(i) + (x(j,t) - threshold) * H(x(j,t) - threshold);
                end
            end
        end
        x(:,t+1) = f(x(:,t), y(:,t) + beta);
        y(:,t+1) = y(:,t) - mu * ((x(:,t) + 1 - sigma) - gc * I_inhib);
    end

    %% Plot membrane traces
    figure('Name', labels{net_idx});
    for i = 1:N
        subplot(N, 1, i);
        plot(1:T, x(i,:), 'k');
        ylim([-3 3]);
        ylabel(sprintf('Neuron %d', i));
    end

    % %% Raster plot
    % subplot(N+1,1,N+1);
    % imagesc(x > 1);
    % colormap(gray);
    % ylabel('Neuron'); xlabel('Time');
    % title('Raster Plot');
    sgtitle(labels{net_idx});
end
