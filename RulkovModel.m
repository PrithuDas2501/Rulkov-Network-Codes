% Rulkov Map Simulation in MATLAB

% Parameters
alpha = 6.0;     % Controls spiking (bursting for alpha > 4)
mu = 0.001;      % Slow variable step size (must be small)
sigma = -0;    % External input
beta = 0.0;      % External input to fast variable (can be used to excite)

% Simulation time
N = 1000;       % Number of iterations

% Preallocate arrays
x = zeros(1, N); % Fast variable (membrane potential)
y = zeros(1, N); % Slow variable

% Initial conditions
x(1) = -1 - 1e-3;
y(1) = -4 - 1e-3;

% Function for fast subsystem
f = @(x, y) (x <= 0) .* (alpha / (1 - x) + y) + ...
           (x > 0 & x < alpha + y) .* (alpha + y) + ...
           (x >= alpha + y) .* (-1);

% Main loop
for n = 1:N-1
    x(n+1) = f(x(n), y(n) + beta);
    y(n+1) = y(n) - mu * (x(n) - sigma + 1);
end

% Plotting
figure(1); clf;

subplot(2,1,1);
plot(1:N, x, 'b');
xlabel('Iteration (n)');
ylabel('x_n (Membrane potential)');
title('Rulkov Map Neuron - Time Series');

subplot(2,1,2);
plot(1:N, y, 'r');
xlabel('Iteration (n)');
ylabel('y_n');
title('Rulkov Map Neuron - Time Series');
