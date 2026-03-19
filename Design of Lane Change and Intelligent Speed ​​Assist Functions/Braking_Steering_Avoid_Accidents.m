%% Brake + Steering Avoidance Analysis
% Based on:
%   s_steer = v_x * sqrt(2*d / a_y,max)
%   s_stop  = v_x^2 / (2*|a_x|)
%   sqrt(a_x^2 + a_y^2) <= mu*g
%
% Advanced engineering-style implementation:
% - command vs realized acceleration separation
% - friction-circle enforcement
% - infeasible-case handling
% - summary table + diagnostic plots
%
% Author intent:
% Motion Planning & Control style MATLAB analysis

clear; clc; close all;

%% Scenario Parameters
P.vx0   = 25.0;                 % [m/s] initial longitudinal speed
P.axCmd = -5.0;                 % [m/s^2] commanded emergency braking
P.muVec = [0.15, 0.50, 0.85];   % friction coefficients
P.g     = 9.81;                 % [m/s^2]
P.d     = 3.5;                  % [m] required lateral displacement (ASSUMPTION)
P.dt    = 0.01;                 % [s] simulation step
P.epsAy = 1e-8;                 % threshold for lateral reserve

%% Run all cases
nCases = numel(P.muVec);
R = repmat(initResultStruct(), nCases, 1);

for k = 1:nCases
    R(k) = simulateCase(P, P.muVec(k));
end

%% Print engineering summary
summaryTbl = buildSummaryTable(R);
disp('=== Engineering Summary ===');
disp(summaryTbl);

%% Plot 1: Steering-only avoidance distance during emergency braking
figure('Name','Steering Avoidance Distance vs Time','Color','w');
hold on; grid on; box on;

finiteMax = 0;
for k = 1:nCases
    if any(isfinite(R(k).sSteer))
        finiteMax = max(finiteMax, max(R(k).sSteer(isfinite(R(k).sSteer))));
    end
end
if finiteMax == 0
    finiteMax = 1;
end

for k = 1:nCases
    if all(~isfinite(R(k).sSteer))
        % No lateral reserve -> steering avoidance distance is effectively infinite
        plot(R(k).t, nan(size(R(k).t)), 'LineWidth', 2.0, ...
            'DisplayName', sprintf('\\mu = %.2f (yanal rezerv yok)', R(k).mu));
        
        % Visual annotation near end of case horizon
        text(R(k).t(end)*0.10, finiteMax*(0.92 - 0.08*k), ...
            sprintf('\\mu=%.2f: |a_x| > \\mu g, a_y rezervi \\approx 0', R(k).mu), ...
            'FontSize', 10);
    else
        plot(R(k).t, R(k).sSteer, 'LineWidth', 2.2, ...
            'DisplayName', sprintf('\\mu = %.2f', R(k).mu));
    end
end

xlabel('Zaman [s]');
ylabel('s_{steer}(t) [m]');
title('Acil Frenleme Sırasında Sadece Dümenleme ile Kaçınma Mesafesi');
legend('Location','best');
ylim([0, 1.15*finiteMax]);

%% Plot 2: Longitudinal speed profile
figure('Name','Longitudinal Speed vs Time','Color','w');
hold on; grid on; box on;

for k = 1:nCases
    plot(R(k).t, R(k).vx, 'LineWidth', 2.0, ...
        'DisplayName', sprintf('\\mu = %.2f, a_{x,real}=%.3f m/s^2', R(k).mu, R(k).axReal));
end

xlabel('Zaman [s]');
ylabel('v_x(t) [m/s]');
title('Fiziksel Olarak Gerçekleşen Frenleme Altında Boylamsal Hız');
legend('Location','best');

%% Plot 3: Remaining stopping distance during braking
figure('Name','Remaining Stopping Distance vs Time','Color','w');
hold on; grid on; box on;

for k = 1:nCases
    plot(R(k).t, R(k).sStop, 'LineWidth', 2.0, ...
        'DisplayName', sprintf('\\mu = %.2f', R(k).mu));
end

xlabel('Zaman [s]');
ylabel('s_{stop}(t) [m]');
title('Acil Frenleme Sırasında Kalan Durma Mesafesi');
legend('Location','best');

%% Optional combined diagnostic plot
figure('Name','Combined Diagnostics','Color','w');
tiledlayout(2,1, 'Padding','compact', 'TileSpacing','compact');

nexttile; hold on; grid on; box on;
for k = 1:nCases
    if any(isfinite(R(k).sSteer))
        plot(R(k).t, R(k).sSteer, 'LineWidth', 2.0, ...
            'DisplayName', sprintf('\\mu = %.2f', R(k).mu));
    end
end
xlabel('Zaman [s]');
ylabel('s_{steer}(t) [m]');
title('Steering Avoidance Distance');
legend('Location','best');

nexttile; hold on; grid on; box on;
for k = 1:nCases
    plot(R(k).t, R(k).sStop, 'LineWidth', 2.0, ...
        'DisplayName', sprintf('\\mu = %.2f', R(k).mu));
end
xlabel('Zaman [s]');
ylabel('s_{stop}(t) [m]');
title('Remaining Stopping Distance');
legend('Location','best');

%% ---------------- Local Functions ----------------
function S = initResultStruct()
S.mu         = NaN;
S.mu_g       = NaN;
S.axCmd      = NaN;
S.axReal     = NaN;
S.ayMax      = NaN;
S.isFeasible = false;
S.t          = [];
S.vx         = [];
S.sStop      = [];
S.sSteer     = [];
S.tStop      = NaN;
S.sStop0     = NaN;
S.sSteer0    = NaN;
end

function S = simulateCase(P, mu)
% Simulate one friction case under commanded emergency braking.
%
% Engineering logic:
% 1) Commanded braking is axCmd.
% 2) Realized braking cannot exceed mu*g.
% 3) Remaining lateral authority comes from friction circle.
% 4) If lateral authority is zero, steering-only avoidance distance is infinite.

    S = initResultStruct();

    S.mu   = mu;
    S.mu_g = mu * P.g;
    S.axCmd = P.axCmd;

    axCmdMag = abs(P.axCmd);
    axRealMag = min(axCmdMag, S.mu_g);   % physically realizable magnitude
    S.axReal = -axRealMag;

    % Feasibility of commanded emergency braking
    S.isFeasible = (axCmdMag <= S.mu_g + 1e-12);

    % Remaining lateral acceleration reserve from friction circle
    ayMaxSq = max(S.mu_g^2 - axRealMag^2, 0.0);
    S.ayMax = sqrt(ayMaxSq);

    % Time to stop with realized deceleration
    if axRealMag < 1e-12
        S.t = 0;
        S.vx = P.vx0;
        S.sStop = inf;
        S.sSteer = inf;
        S.tStop = inf;
        S.sStop0 = inf;
        S.sSteer0 = inf;
        return;
    end

    S.tStop = P.vx0 / axRealMag;
    S.t = 0:P.dt:S.tStop;
    if abs(S.t(end) - S.tStop) > 1e-12
        S.t = [S.t, S.tStop];
    end

    % Longitudinal speed profile
    S.vx = max(P.vx0 + S.axReal * S.t, 0.0);

    % Remaining stopping distance during braking
    S.sStop = (S.vx.^2) ./ (2*axRealMag);
    S.sStop0 = S.sStop(1);

    % Steering-only avoidance distance during braking
    % s_steer(t) = v_x(t) * sqrt(2*d / a_y,max)
    if S.ayMax > P.epsAy
        S.sSteer = S.vx .* sqrt(2*P.d / S.ayMax);
        S.sSteer0 = S.sSteer(1);
    else
        % No lateral reserve -> effectively impossible by steering alone
        S.sSteer = nan(size(S.t));
        S.sSteer0 = inf;
    end
end

function T = buildSummaryTable(R)
    n = numel(R);

    mu       = zeros(n,1);
    mu_g     = zeros(n,1);
    axCmd    = zeros(n,1);
    axReal   = zeros(n,1);
    ayMax    = zeros(n,1);
    feasible = false(n,1);
    tStop    = zeros(n,1);
    sStop0   = zeros(n,1);
    sSteer0  = zeros(n,1);

    for i = 1:n
        mu(i)       = R(i).mu;
        mu_g(i)     = R(i).mu_g;
        axCmd(i)    = R(i).axCmd;
        axReal(i)   = R(i).axReal;
        ayMax(i)    = R(i).ayMax;
        feasible(i) = R(i).isFeasible;
        tStop(i)    = R(i).tStop;
        sStop0(i)   = R(i).sStop0;
        sSteer0(i)  = R(i).sSteer0;
    end

    T = table(mu, mu_g, axCmd, axReal, ayMax, feasible, tStop, sStop0, sSteer0, ...
        'VariableNames', {'mu','mu_g','ax_cmd','ax_real','ay_max','ax_cmd_feasible','t_stop','s_stop_t0','s_steer_t0'});
end