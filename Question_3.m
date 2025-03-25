
% Question 3 

% Parameters
n = 7; % Number of states
sigma = 1; % Standard deviation of the error term
seed = 2025; % Seed for simulation
T = 100; % Number of periods to simulate

function [y, P] = rouwenhorst(n, gamma1, sigma)
    s = sigma * (sqrt((n - 1) / (1 - gamma1^2)));
    mu = 0.5 / (1 - gamma1);
    y = linspace(mu - s, mu + s, n);
    
    if n == 2
        P = [(1 + gamma1) / 2, (1 - gamma1) / 2;
             (1 - gamma1) / 2, (1 + gamma1) / 2];
    else
        [~, Q] = rouwenhorst(n - 1, gamma1, sigma);
        P = ((1 + gamma1) / 2) * [Q, zeros(n-1, 1); zeros(1, n)];
        P = P + ((1 - gamma1) / 2) * [zeros(n-1, 1), Q; zeros(1, n)];
        P = P + ((1 - gamma1) / 2) * [zeros(1, n); Q, zeros(n-1, 1)];
        P = P + ((1 + gamma1) / 2) * [zeros(1, n); zeros(n-1, 1), Q];
       P(2:end-1, :) = P(2:end-1, :) / 2;
    end
end

% (b) γ1 = 0.85
gamma1_b = 0.85;
[y_b, P_b] = rouwenhorst(n, gamma1_b, sigma);
disp('Transition Matrix (γ1 = 0.85):');
disp(P_b);
disp('State Vector (γ1 = 0.85):');
disp(y_b);

% Simulation function
function y = simulate(grid, pmat, T)
    n = size(grid, 1);
    state0 = randsample(n, 1);
    cmat = cumsum(pmat, 2);
    y = zeros(1, T); 
    for i = 1:T
        y(i) = grid(state0); 
        state1 = find(rand <= cmat(state0, :));
        if ~isempty(state1)
            state0 = state1(1);
        else
        end
    end
end

rng(2025); % Set seed
y_value = simulate(y_b, P_b, T);
figure;
plot(1:T, y_value, 'b', 'LineWidth', 1.5);
title('Simulated Markov Chains for γ1 = 0.85');
xlabel('Periods');
ylabel('States');
legend('Location', 'best');
grid on;

% (d) Repeat (c) for γ1 = 0.75, 0.85, 0.95, 0.99
gamma1_values = [0.75, 0.85, 0.95, 0.99];
figure;
hold on;

for i = 1:length(gamma1_values)
    gamma1_d = gamma1_values(i);
    [y_d, P_d] = rouwenhorst(n, gamma1_d, sigma);
    rng(seed); % Reset seed for each simulation
    state0 = randi(n); % Use same initial state for each gamma
    y_values_d = simulate(y_d, P_d, T);
    plot(1:T, y_values_d, 'LineWidth', 1.5);
end

title('Simulated Markov Chains for Different γ1');
xlabel('Periods');
ylabel('States');
legend('γ1 = 0.75', 'γ1 = 0.85', 'γ1 = 0.95', 'γ1 = 0.99', 'Location', 'best');
grid on;
hold off;