%% NOTE:
% This script analyzes an aggregate general equilibrium model with capital for autonomous long haul trucks.

clear
clc

%% Parameters
D = 0.1; % disutility
Ay = 1; % productivity of good Y
As = 1; % productivity of short haul
Al = 1; % productivity of long haul
alpha = 0.5; % the share of Y
P = 1; % the price of good Q (normalize to 1)
beta_range = linspace(0,1,50); % the share of autonomous truck
Pk_range = linspace(0.5,1.5,50); % the price of capital for autonomous truck


%% short run labor market
Ly = 1;
Ls = 0.75;
Ll = 0.25;

% Results
Kl_s = zeros(length(beta_range), length(Pk_range));
Wy_s = zeros(length(beta_range), length(Pk_range));
Ws_s = zeros(length(beta_range), length(Pk_range));
Wl_s = zeros(length(beta_range), length(Pk_range));
Q_s = zeros(length(beta_range), length(Pk_range));

Uy_s = zeros(length(beta_range), length(Pk_range));
Us_s = zeros(length(beta_range), length(Pk_range));
Ul_s = zeros(length(beta_range), length(Pk_range));


% Short-run solution
for i = 1:length(beta_range)
    for j = 1:length(Pk_range)
        Kl_s(i,j) = (((1-alpha) * beta_range(i) * (Ay * Ly) ^ alpha * (Al * Ll ^ (1 - beta_range(i))) ^ (1 - alpha)) / Pk_range(j)) ^ (1 / (1 - (1 - alpha) * beta_range(i))) ; % Autonomous Trucks
        kl_s(i,j) = Kl_s(i,j) / Ly; % Autonomous truck per capita
        Q_s(i,j) = (Ay * Ly) ^ alpha * (Al * Kl_s(i,j) ^ beta_range(i) * Ll ^ (1 - beta_range(i))) ^ (1 - alpha);
        Wy_s(i,j) = alpha * Q_s(i,j) / Ly;
        Ws_s(i,j) = (1 - alpha) * Q_s(i,j) / Ls;
        Wl_s(i,j) = (1 - alpha) * (1 - beta_range(i)) * Q_s(i,j) / Ll;
        
        Uy_s(i,j) = log((Wy_s(i,j) + Pk_range(j) * kl_s(i,j)) / P);
        Us_s(i,j) = log(Ws_s(i,j) / P);
        Ul_s(i,j) = log(Wl_s(i,j) / P) - D;
    end
end

%% long run labor market
L_bar = 2;

% Results
Wy_l = zeros(length(beta_range), length(Pk_range));
Ws_l = zeros(length(beta_range), length(Pk_range));
Wl_l = zeros(length(beta_range), length(Pk_range));
Ly_l = zeros(length(beta_range), length(Pk_range));
Ls_l = zeros(length(beta_range), length(Pk_range));
Ll_l = zeros(length(beta_range), length(Pk_range));
Kl_l = zeros(length(beta_range), length(Pk_range));

Uy_l = zeros(length(beta_range), length(Pk_range));
Us_l = zeros(length(beta_range), length(Pk_range));
Ul_l = zeros(length(beta_range), length(Pk_range));

function F = longrun(x_log, params)
    x = exp(x_log);  % Log-to-level transformation
    % Unpack variables
    W = x(1:3); 
    % W(1): Wage in producing Y
    % W(2): Wage in local transportation
    % W(3): Wage in interstate transportation

    L = x(4:6); % Labor allocations
    % L(1): Labor allocated to the good sector
    % L(2): Labor allocated to the local transportation sector 
    % L(3): Labor allocated to the interstate transportation sector

    K = x(7); % Autonomous Trucks

    % Unpack parameters
    Ay = params.Ay;
    As = params.As;
    Al = params.Al;
    alpha = params.alpha;
    beta = params.beta;
    Pk = params.Pk;
    P = params.P;
    L_bar = params.L_bar;
    D = params.D;

    % Equations
    F(1) = alpha * (Ay * L(1)) ^ alpha * (As * L(2)) ^ (1 - alpha) - W(1) * L(1);
    F(2) = (1 - alpha) * (Ay * L(1)) ^ alpha * (As * L(2)) ^ (1 - alpha) - W(2) * L(2);
    F(3) = (1 - alpha) * (1 - beta) * (Ay * L(1)) ^ alpha * (Al * K ^ beta * L(3) ^ (1 - beta)) ^ (1 - alpha) - W(3) * L(3);
    F(4) = (1 - alpha) * beta * (Ay * L(1)) ^ alpha * (Al * K ^ beta * L(3) ^ (1 - beta)) ^ (1 - alpha) - K * Pk;

    % Utility equalization
    F(5) = log((W(1) + Pk * K / L(1)) / P) - log(W(2) / P);
    F(6) = log((W(1) + Pk * K / L(1)) / P) - log(W(3) / P) + D;

    % Labor market clearing
    F(7) = L(1) + L(2) + L(3) - L_bar;
end

for i = 1:length(beta_range)
    for j = 1:length(Pk_range)
        
        % set parameter
        params.Ay = Ay;
        params.As = As;
        params.Al = Al;
        params.alpha = alpha;
        params.beta = beta_range(i);
        params.Pk = Pk_range(j);
        params.P = P;
        params.L_bar = L_bar;
        params.D = D;
        
        % Initial guesses for the variables
        x0 = log([0.5 * ones(3, 1); % Wage
          1 ; 0.75; 0.25; % Labor allocation
          0.1]); % Capital
        % Solve the system of equations
        options = optimoptions('fsolve', ...
        'Display', 'iter', ...
        'TolFun', 1e-16, ...
        'TolX', 1e-16, ...
        'MaxIterations', 10000, ...
        'MaxFunctionEvaluations', 50000);
        x_sol = fsolve(@(x) longrun(x, params), x0, options);
        x = exp(x_sol);

        % Store the solution
        Wy_l(i,j) = x(1);
        Ws_l(i,j) = x(2);
        Wl_l(i,j) = x(3);
        Ly_l(i,j) = x(4);
        Ls_l(i,j) = x(5);
        Ll_l(i,j) = x(6);
        Kl_l(i,j) = x(7);
        Q_l(i,j) = (Ay * x(4)) ^ alpha * (As * x(5)) ^ (1 - alpha);

        Uy_l(i,j) = log((x(1) + Pk_range(j) * x(7) / x(4)) / P);
        Us_l(i,j) = log(x(2) / P);
        Ul_l(i,j) = log(x(3) / P) - D;
    end
end

%% Plot Utility vs. Beta (Immobile Labor)
% Set the value of Pk you want to plot
Pk_plot = 1; % <-- Change this value to plot for a different Pk (0.5 - 1.5)

% Find the index where Pk is closest to Pk_plot
[~, pk_idx] = min(abs(Pk_range - Pk_plot));

% Prepare data for plotting
beta = beta_range;
Uy = Uy_s(:, pk_idx);
Us = Us_s(:, pk_idx);
Ul = Ul_s(:, pk_idx);

% Plot
figure;
plot(beta, Uy, 'r-', 'LineWidth', 2); hold on;
plot(beta, Us, 'b--', 'LineWidth', 2);
plot(beta, Ul, 'g-.', 'LineWidth', 2);
hold off;

xlabel('\beta (Share of Autonomous Trucks)');
ylabel('Utility');
title(['Utility vs. Share of Autonomous Trucks (\beta), P_k = ' num2str(Pk_plot) ' (Immobile Labor)']);
legend('Goods Sector (Uy)', 'Short-haul (Us)', 'Long-haul (Ul)', 'Location', 'Best');
grid on;
saveas(gcf, 'Utility_vs_Beta_Immobile.png');

%% Plot Utility vs. Pk (Immobile Labor) for a given beta

% Set the value of beta you want to plot
beta_plot = 0.5; % <-- Change this value to plot for a different beta (0 to 1)

% Find the index where beta is closest to beta_plot
[~, beta_idx] = min(abs(beta_range - beta_plot));

% Prepare data for plotting
Pk = Pk_range;
Uy = Uy_s(beta_idx, :);
Us = Us_s(beta_idx, :);
Ul = Ul_s(beta_idx, :);

% Plot
figure;
plot(Pk, Uy, 'r-', 'LineWidth', 2); hold on;
plot(Pk, Us, 'b--', 'LineWidth', 2);
plot(Pk, Ul, 'g-.', 'LineWidth', 2);
hold off;

xlabel('P_k (Price of Capital)');
ylabel('Utility');
title(['Utility vs. Price of Capital (P_k), \beta = ' num2str(beta_plot) ' (Immobile Labor)']);
legend('Goods Sector (Uy)', 'Short-haul (Us)', 'Long-haul (Ul)', 'Location', 'Best');
grid on;
saveas(gcf, 'Utility_vs_Pk_Immobile.png');

%% Plot Utility vs. Beta (Mobile Labor)
% Set the value of Pk you want to plot
Pk_plot = 1; % <-- Change this value to plot for a different Pk (0.5 - 1.5)

% Find the index where Pk is closest to Pk_plot
[~, pk_idx] = min(abs(Pk_range - Pk_plot));

% Prepare data for plotting
beta = beta_range;
Uy = Uy_l(:, pk_idx);
Us = Us_l(:, pk_idx);
Ul = Ul_l(:, pk_idx);

% Plot
figure;
plot(beta, Uy, 'r-', 'LineWidth', 2); hold on;
plot(beta, Us, 'b--', 'LineWidth', 2);
plot(beta, Ul, 'g-.', 'LineWidth', 2);
hold off;

xlabel('\beta (Share of Autonomous Trucks)');
ylabel('Utility');
title(['Utility vs. Share of Autonomous Trucks (\beta), P_k = ' num2str(Pk_plot) ' (Mobile Labor)']);
legend('Goods Sector (Uy)', 'Short-haul (Us)', 'Long-haul (Ul)', 'Location', 'Best');
grid on;
saveas(gcf, 'Utility_vs_Beta_Mobile.png');

%% Plot Utility vs. Pk (Mobile Labor) for a given beta

% Set the value of beta you want to plot
beta_plot = 0.5; % <-- Change this value to plot for a different beta (0 to 1)

% Find the index where beta is closest to beta_plot
[~, beta_idx] = min(abs(beta_range - beta_plot));

% Prepare data for plotting
Pk = Pk_range;
Uy = Uy_l(beta_idx, :);
Us = Us_l(beta_idx, :);
Ul = Ul_l(beta_idx, :);

% Plot
figure;
plot(Pk, Uy, 'r-', 'LineWidth', 2); hold on;
plot(Pk, Us, 'b--', 'LineWidth', 2);
plot(Pk, Ul, 'g-.', 'LineWidth', 2);
hold off;

xlabel('P_k (Price of Capital)');
ylabel('Utility');
title(['Utility vs. Price of Capital (P_k), \beta = ' num2str(beta_plot) ' (Mobile Labor)']);
legend('Goods Sector (Uy)', 'Short-haul (Us)', 'Long-haul (Ul)', 'Location', 'Best');
grid on;
saveas(gcf, 'Utility_vs_Pk_Mobile.png');

%% Plot Labor Allocations vs. Beta (Mobile Labor)
% Set the value of Pk you want to plot
Pk_plot = 1; % <-- Change this value to plot for a different Pk

% Find the index where Pk is closest to Pk_plot
[~, pk_idx] = min(abs(Pk_range - Pk_plot));

% Prepare data for plotting
beta = beta_range;
Ly = Ly_l(:, pk_idx);
Ls = Ls_l(:, pk_idx);
Ll = Ll_l(:, pk_idx);

% Plot
figure;
plot(beta, Ly, 'r-', 'LineWidth', 2); hold on;
plot(beta, Ls, 'b--', 'LineWidth', 2);
plot(beta, Ll, 'g-.', 'LineWidth', 2);
hold off;

xlabel('\beta (Share of Autonomous Trucks)');
ylabel('Labor Allocation');
title(['Labor Allocation vs. Share of Autonomous Trucks (\beta), P_k = ' num2str(Pk_plot) ' (Mobile Labor)']);
legend('Goods Sector (Ly)', 'Short-haul (Ls)', 'Long-haul (Ll)', 'Location', 'Best');
grid on;
saveas(gcf, 'Labor_vs_Beta_Mobile.png');

%% Plot Labor Allocations vs. Pk (Mobile Labor) for a given beta

% Set the value of beta you want to plot
beta_plot = 0.5; % <-- Change this value to plot for a different beta

% Find the index where beta is closest to beta_plot
[~, beta_idx] = min(abs(beta_range - beta_plot));

% Prepare data for plotting
Pk = Pk_range;
Ly = Ly_l(beta_idx, :);
Ls = Ls_l(beta_idx, :);
Ll = Ll_l(beta_idx, :);

% Plot
figure;
plot(Pk, Ly, 'r-', 'LineWidth', 2); hold on;
plot(Pk, Ls, 'b--', 'LineWidth', 2);
plot(Pk, Ll, 'g-.', 'LineWidth', 2);
hold off;

xlabel('P_k (Price of Capital)');
ylabel('Labor Allocation');
title(['Labor Allocation vs. Price of Capital (P_k), \beta = ' num2str(beta_plot) ' (Mobile Labor)']);
legend('Goods Sector (Ly)', 'Short-haul (Ls)', 'Long-haul (Ll)', 'Location', 'Best');
grid on;
saveas(gcf, 'Labor_vs_Pk_Mobile.png');
