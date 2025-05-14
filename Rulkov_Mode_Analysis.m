% Rulkov Network Mode Analysis (based on working version)
% --------------------------------------------------------
% Computes eigenmodes of the 4x4 lattice and visualizes which mode becomes
% unstable first as g_c increases. Also overlays dominant mode on simulated bursts.

clear; clc;

%% Parameters
N = 16;
T = 10000;
alpha = 6.0;
mu = 0.001;
sigma = 0.0;
beta = 0.0;
threshold = -2.0;
gc_vals = [-0.05, -0.15, -0.20, -0.40, -0.49, -0.55];

%% Rulkov map dynamics
f = @(x, y) (x <= 0) .* (alpha ./ (1 - x) + y) + ...
           (x > 0 & x < alpha + y) .* (alpha + y) + ...
           (x >= alpha + y) .* (-1);
H = @(v) double(v > 0);

%% Generate 4x4 lattice
A = generate_lattice_4x4();
[vecs, D] = eig(A); % Eigenvectors and eigenvalues of A
lambdas = diag(D);

%% Fixed point estimation
x_star = -1.5;
y_star = -2;
a = (alpha / (1 - x_star)^2);
b = 1;

%% Compute stability matrices and dominant mode
fprintf('gc\t\tMax|eig|\t\tUnstable Mode Index\n');
for g = gc_vals
    unstable_modes = [];
    for k = 1:N
        M = [a, b; -mu * (1 - g * lambdas(k)), 1];
        eigsM = eig(M);
        maxmod = max(abs(eigsM));
        if maxmod > 1
            unstable_modes(end+1) = k; %#ok<AGROW>
        end
    end
    if isempty(unstable_modes)
        fprintf('%.2f\t\tStable\n', g);
    else
        fprintf('%.2f\t\t>1\t\tMode(s): %s\n', g, mat2str(unstable_modes));
    end
end

%% Simulate and overlay dominant mode pattern
for idx = 1:length(gc_vals)
    gc = gc_vals(idx);

    x = -1.5 + 0.01*rand(N, T);
    y = -2 + 0.01*rand(N, T);

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

    % Determine dominant mode
    mode_max = 0;
    max_mag = 0;
    for k = 1:N
        M = [a, b; -mu * (1 - gc * lambdas(k)), 1];
        if max(abs(eig(M))) > max_mag
            max_mag = max(abs(eig(M)));
            mode_max = k;
        end
    end

    pattern = vecs(:, mode_max);

    % Plot simulation vs predicted dominant mode
    figure('Name', sprintf('g_c = %.2f | Dominant Mode %d', gc, mode_max));
    for i = 1:N
        subplot(4, 4, i);
        plot(1:T, x(i,:), 'k'); hold on;
        if pattern(i) > 0.1
            fill([1 T T 1], [3 3 -3 -3], [0.8 0.9 1], 'EdgeColor', 'none');
            plot(1:T, x(i,:), 'k');
        elseif pattern(i) < -0.1
            fill([1 T T 1], [3 3 -3 -3], [1 0.9 0.9], 'EdgeColor', 'none');
            plot(1:T, x(i,:), 'k');
        end
        ylim([-3 3]);
        title(sprintf('Neuron %d', i));
    end
    sgtitle(sprintf('g_c = %.2f | Predicted Pattern from Mode %d', gc, mode_max));
end