%% NOTE:
% This script analyzes an aggregate general equilibrium model with capital for autonomous long haul trucks.

clear
clc

%% Parameters
D = 0.9; % disutility
Ay = 3; % productivity of good Y
As = 1; % productivity of short haul
Al = 1; % productivity of long haul
alpha = 0.5; % the share of Y
gamma = 0.96; % discount factor
delta = 0.1; % depreciation rate
L_bar = 1; % labor supply
beta_range = linspace(0.1,0.9,50); % the share of autonomous truck


%% Solutions
Pk = (1 / gamma) + delta - 1;
Wy = zeros(1,length(beta_range));
Ws = zeros(1,length(beta_range));
Wl = zeros(1,length(beta_range));
Ly = zeros(1,length(beta_range));
Ls = zeros(1,length(beta_range));
Ll = zeros(1,length(beta_range));
K = zeros(1,length(beta_range));
C = zeros(1,length(beta_range));
Q = zeros(1,length(beta_range));
U = zeros(1,length(beta_range));
EXITFLAG = zeros(1,length(beta_range));
FVAL = zeros(1,length(beta_range));


function F = AT(x_log, params, W1)
  x = exp(x_log);  % Log-to-level transformation
  % Unpack variables
  W = [W1; x(1:2)]; 
  % W(1): Wage in producing Y (fixed)
  % W(2): Wage in local transportation
  % W(3): Wage in interstate transportation

  L = x(3:5); % Labor allocations
  % L(1): Labor allocated to the good sector
  % L(2): Labor allocated to the local transportation sector 
  % L(3): Labor allocated to the interstate transportation sector

  K = x(6); % Autonomous Trucks

  C = x(7); % Consumption

  % Unpack parameters
  D = params.D;
  Ay = params.Ay;
  As = params.As;
  Al = params.Al;
  alpha = params.alpha;
  delta = params.delta;
  L_bar = params.L_bar;
  Pk = params.Pk;
  beta = params.beta;

  % Equations (labor market clearance removed)
  F(1) = W(3) - C * D;
  F(2) = C - ((Ay * L(1)) ^ alpha * (min(As * L(2), Al * K ^ beta * L(3) ^ (1 - beta))) ^ (1 - alpha) - delta * K);
  F(3) = W(3) * L(3) - (1 - alpha) * (1 - beta) * (Ay * L(1)) ^ alpha * (min(As * L(2), Al * K ^ beta * L(3) ^ (1 - beta))) ^ (1 - alpha);
  F(4) = W(2) * L(2) - (1 - alpha) * (Ay * L(1)) ^ alpha * (min(As * L(2), Al * K ^ beta * L(3) ^ (1 - beta))) ^ (1 - alpha);
  F(5) = W(1) * L(1) - alpha * (Ay * L(1)) ^ alpha * (min(As * L(2), Al * K ^ beta * L(3) ^ (1 - beta))) ^ (1 - alpha);
  F(6) = W(1) - W(2);
  F(7) = Pk * K - (1 - alpha) * beta * (Ay * L(1)) ^ alpha * (min(As * L(2), Al * K ^ beta * L(3) ^ (1 - beta))) ^ (1 - alpha);
end


for i = 1:length(beta_range)
  % set parameter
  params.D = D;
  params.Ay = Ay;
  params.As = As;
  params.Al = Al;
  params.alpha = alpha;
  params.delta = delta;
  params.L_bar = L_bar;
  params.Pk = Pk;
  params.beta = beta_range(i);

  % Initial guess for wage (W1)
  W1 = 1;
  tol = 1e-6;
  max_iter = 50;
  iter = 0;
  excess = 1;
  while abs(excess) > tol && iter < max_iter
    % Initial guesses for the remaining variables
    x0 = log([1; 1; 1; 1; 1; 1; 1]); % [W2, W3, L1, L2, L3, K, C]
    options = optimoptions('fsolve', ...
      'Display', 'off', ...
      'TolFun', 1e-10, ...
      'TolX', 1e-10, ...
      'MaxIterations', 1000, ...
      'MaxFunctionEvaluations', 5000);
    [x_sol, fval, exitflag] = fsolve(@(x) AT(x, params, W1), x0, options);
    x = exp(x_sol);
    W2 = x(1); W3 = x(2);
    L1 = x(3); L2 = x(4); L3 = x(5);
    Kval = x(6); Cval = x(7);
    excess = (L1 + L2 + L3) - L_bar;
    % Adjust W1 based on excess demand/supply
    W1 = W1 * (1 + 0.5 * excess / L_bar);
    iter = iter + 1;
  end
  Wy(i) = W1;
  Ws(i) = W2;
  Wl(i) = W3;
  Ly(i) = L1;
  Ls(i) = L2;
  Ll(i) = L3;
  K(i) = Kval;
  C(i) = Cval;
  Q(i) = (Ay * L1) ^ alpha * (min(As * L2, Al * Kval ^ beta_range(i) * L3 ^ (1 - beta_range(i)))) ^ (1 - alpha);
  U(i) = ((1 - gamma) / gamma) * (log(C(i)) - Ll(i) * D);
  EXITFLAG(i) = exitflag;
  FVAL(i) = norm(fval);
end



%% Plot Consumption vs. Beta
% Prepare data for plotting
beta = beta_range;

% Plot
figure;
plot(beta, C, 'r-', 'LineWidth', 2); hold on;
hold off;

xlabel('\beta (Share of Autonomous Trucks)');
ylabel('Consumption');
title('Consumption vs. Share of Autonomous Trucks (\beta)');
legend('Consumption', 'Location', 'Best');
grid on;
saveas(gcf, 'C_vs_Beta.png');

%% Plot Production vs. Beta
% Prepare data for plotting
beta = beta_range;

% Plot
figure;
plot(beta, Q, 'r-', 'LineWidth', 2); hold on;
hold off;

xlabel('\beta (Share of Autonomous Trucks)');
ylabel('Production');
title('Production vs. Share of Autonomous Trucks (\beta)');
legend('Production', 'Location', 'Best');
grid on;
saveas(gcf, 'Q_vs_Beta.png');


%% Plot Labor allocation vs. Beta
% Prepare data for plotting
beta = beta_range;

% Plot
figure;
plot(beta, Ly, 'r-', 'LineWidth', 2); hold on;
plot(beta, Ls, 'b--', 'LineWidth', 2);
plot(beta, Ll, 'g-.', 'LineWidth', 2);
hold off;

xlabel('\beta (Share of Autonomous Trucks)');
ylabel('Labor allocations');
title('Labor allocations vs. Share of Autonomous Trucks (\beta)');
legend('Goods Sector (Ly)', 'Short-haul (Ls)', 'Long-haul (Ll)', 'Location', 'Best');
grid on;
saveas(gcf, 'Labor_vs_Beta.png');

%% Plot Wage vs. Beta
% Prepare data for plotting
beta = beta_range;

% Plot
figure;
plot(beta, Wy, 'r-', 'LineWidth', 2); hold on;
plot(beta, Ws, 'b--', 'LineWidth', 2);
plot(beta, Wl, 'g-.', 'LineWidth', 2);
hold off;

xlabel('\beta (Share of Autonomous Trucks)');
ylabel('Wage');
title('Wage vs. Share of Autonomous Trucks (\beta)');
legend('Goods Sector (Wy)', 'Short-haul (Ws)', 'Long-haul (Wl)', 'Location', 'Best');
grid on;
saveas(gcf, 'Wage_vs_Beta.png');

%% Plot Utility vs. Beta
% Prepare data for plotting
beta = beta_range;

% Plot
figure;
plot(beta, U, 'r-', 'LineWidth', 2); hold on;
hold off;

xlabel('\beta (Share of Autonomous Trucks)');
ylabel('Lifetime Utility');
title('Lifetime Utility vs. Share of Autonomous Trucks (\beta)');
legend('Lifetime Utility');
grid on;
saveas(gcf, 'U_vs_Beta.png');

EXITFLAG
FVAL