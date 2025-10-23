%% NOTE:
% This script analyzes an aggregate general equilibrium model with capital for autonomous long haul trucks.

clear
clc
close all

%% Parameters
gamma = 0.96; % discount factor
D = 1; % disutility
Ay = 7; % productivity of good Y
As = 1; % productivity of short haul
Al = 1; % productivity of long haul
alpha = 0.5; % the share of Y
beta_range = linspace(0.1,0.9,100); % the share of autonomous truck
delta = 0.05; % depreciation rate
L_bar = 1; % labor supply
Pk = 1 / gamma + delta - 1;


%% Solutions
Wy = zeros(1,length(beta_range));
Ws = zeros(1,length(beta_range));
Wl = zeros(1,length(beta_range));
PT = zeros(1,length(beta_range));
Ly = zeros(1,length(beta_range));
Ls = zeros(1,length(beta_range));
Ll = zeros(1,length(beta_range));
K = zeros(1,length(beta_range));
C = zeros(1,length(beta_range));
Q = zeros(1,length(beta_range));
U = zeros(1,length(beta_range));
EXITFLAG = zeros(1,length(beta_range));
FVAL = zeros(1,length(beta_range));

function F = AT(x_log, params)
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

    C = x(8); % Consumption

    PT = x(9); % Price index in transportation sector

    % Unpack parameters
    D = params.D;
    Ay = params.Ay;
    As = params.As;
    Al = params.Al;
    alpha = params.alpha;
    beta = params.beta;
    delta = params.delta;
    L_bar = params.L_bar;
    Pk = params.Pk;

    % Equations
    F(1) = W(3) - W(1) - C * D;
    F(2) = W(1) - W(2);
    F(3) = As * L(2) - Al * K ^ beta * L(3) ^ (1 - beta);
    Q_here = (Ay * L(1)) ^ alpha * (As * L(2)) ^ (1 - alpha); 
    F(4) = C + delta * K - Q_here;
    F(5) = W(1) - alpha * Q_here / L(1);
    F(6) = W(2) - PT * As;
    F(7) = Pk * K - beta * PT * Al * K ^ beta * L(3) ^ (1 - beta);
    T_here = As * L(2);
    F(8) = PT - (1 - alpha) * Q_here / T_here;
    F(9) = L(1) + L(2) + L(3) - L_bar;
end

% Initial guesses for the variables
x0 = log([1;1;1;   1/3;1/3;1/3;   1;   1;   1]);
%           W's     Ly  Ls  Ll    K    C    PT

x_guess = x0;
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


        % Solve the system of equations
        options = optimoptions('fsolve', ...
        'Display','off', ...
        'FunctionTolerance',1e-10, ...
        'StepTolerance',1e-10, ...
        'MaxIterations',20000, ...
        'MaxFunctionEvaluations',100000);
        [x_sol, fval, exitflag] = fsolve(@(x) AT(x, params), x_guess, options);
        x_guess = x_sol;  
        x = exp(x_sol);

        % Store the solution
        Wy(i) = x(1);
        Ws(i) = x(2);
        Wl(i) = x(3);
        Ly(i) = x(4);
        Ls(i) = x(5);
        Ll(i) = x(6);
        K(i) = x(7);
        C(i) = x(8);
        PT(i) = x(9);
        Q(i) = (Ay * x(4)) ^ alpha * (As * x(5)) ^ (1 - alpha);
        U(i) = (log(C(i)) - Ll(i) * D) / (1 - gamma);

        % exitflag
        EXITFLAG(i) = exitflag
        FVAL(i) = norm(fval)
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

