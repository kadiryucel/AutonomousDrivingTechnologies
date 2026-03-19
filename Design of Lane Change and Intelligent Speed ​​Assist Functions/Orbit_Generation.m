%% Longitudinal Trajectory Generation - Advanced Engineering Layout
% Quintic polynomial candidate generation for longitudinal motion planning
% Based on the provided equations:
%
% s(t) = a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5
%
% Initial state:
% S0 = [s0, sdot0, sddot0]
%
% Terminal candidate state for each (i,j):
% [s1, sdot1, sddot1]_(i,j) = [s_target(Tj)+Delta_s_i, sdot_target(Tj), sddot_target(Tj)]
%
% Target path:
% s_target(t) = s_lv(t) - [D0 + tau * sdot_lv(t)]
%
% Cost:
% Ct = kj*Jt + kt*T + ks*(s1 - sd)^2
% Jt = integral( jerk(t)^2 dt )
%
% Notes:
% - Equation (4) in the image is interpreted as a quintic polynomial in time.
% - sd is taken as the desired target terminal position s_target(T).

clear; clc; close all;

%% =========================
% 1) Scenario Parameters
% ==========================
P = struct();

% Ego vehicle initial state
P.ego.s0 = 0.0;      % [m]
P.ego.v0 = 25.0;     % [m/s]
P.ego.a0 = 0.0;      % [m/s^2]

% Lead vehicle state
P.lv.s0 = 50.0;      % [m]
P.lv.v0 = 20.0;      % [m/s]
P.lv.a0 = 0.0;       % [m/s^2]

% Following policy parameters
P.follow.D0  = 5.0;  % [m] standstill spacing
P.follow.tau = 1.5;  % [s] time gap

% Candidate grid
P.grid.Tset  = 2.0 : 0.5 : 6.0;     % [s] terminal time candidates
P.grid.dSset = -15.0 : 1.0 : 10.0;  % [m] terminal position offsets

% Numerical sampling
P.dt = 0.02;  % [s]

% Feasibility limits (engineering guardrails)
P.lim.vMin = 0.0;    % [m/s]
P.lim.vMax = 40.0;   % [m/s]
P.lim.aMin = -8.0;   % [m/s^2]
P.lim.aMax =  3.0;   % [m/s^2]
P.lim.jMax = 15.0;   % [m/s^3]
P.lim.minGap = 0.0;  % [m] collision-free constraint

%% =========================
% 2) Weight Cases for Analysis
% ==========================
% We analyze the effect of kj, kt, ks using multiple scenarios.
W = [];

W(1).name = 'Baseline';
W(1).kj   = 1.0;
W(1).kt   = 0.2;
W(1).ks   = 1.0;

W(2).name = 'High_kj';
W(2).kj   = 10.0;
W(2).kt   = 0.2;
W(2).ks   = 1.0;

W(3).name = 'High_kt';
W(3).kj   = 1.0;
W(3).kt   = 3.0;
W(3).ks   = 1.0;

W(4).name = 'High_ks';
W(4).kj   = 1.0;
W(4).kt   = 0.2;
W(4).ks   = 10.0;

%% =========================
% 3) Plan for Each Weight Set
% ==========================
Results(1) = planLongitudinalTrajectories(P, W(1));

for c = 2:numel(W)
    Results(c) = planLongitudinalTrajectories(P, W(c));
end

%% =========================
% 4) Summary Table
% ==========================
summaryTbl = buildSummaryTable(Results);
disp(' ');
disp('========== OPTIMAL TRAJECTORY SUMMARY ==========');
disp(summaryTbl);

%% =========================
% 5) Plots
% ==========================
plotCandidateBundles(P, Results);
plotOptimalComparisons(P, Results);

%% =========================
% 6) Console Interpretation
% ==========================
disp(' ');
disp('========== INTERPRETATION ==========');
disp('- High kj  -> smoother / lower jerk trajectories are preferred.');
disp('- High kt  -> shorter horizon T is preferred, usually more aggressive behavior.');
disp('- High ks  -> terminal position stays closer to s_target(T), i.e. smaller Delta_s is preferred.');
disp('- Inspect the summary table and plots together for final engineering interpretation.');

%% =====================================================================
%                          LOCAL FUNCTIONS
% =====================================================================

function R = planLongitudinalTrajectories(P, W)
% Generate candidate trajectories over (Delta_s, T), evaluate them,
% and pick the minimum-cost feasible one.

    cand = [];
    idx = 0;

    s0 = P.ego.s0;
    v0 = P.ego.v0;
    a0 = P.ego.a0;

    for T = P.grid.Tset
        [sDes, vDes, aDes] = targetTerminalState(P, T);

        for dS = P.grid.dSset
            idx = idx + 1;

            s1 = sDes + dS;
            sd = sDes;  % desired target terminal position

            coeff = solveQuinticBoundaryValue(s0, v0, a0, s1, vDes, aDes, T);
            traj  = sampleQuinticTrajectory(coeff, T, P.dt);

            % Lead vehicle over same time horizon
            lv = leadTrajectory(P, traj.t);

            % Cost
            Jt = trapz(traj.t, traj.j.^2);
            Ct = W.kj*Jt + W.kt*T + W.ks*(s1 - sd)^2;

            % Feasibility
            feas = checkFeasibility(traj, lv, P);

            cand(idx).T        = T;
            cand(idx).dS       = dS;
            cand(idx).sDes     = sDes;
            cand(idx).vDes     = vDes;
            cand(idx).aDes     = aDes;
            cand(idx).s1       = s1;
            cand(idx).coeff    = coeff;
            cand(idx).traj     = traj;
            cand(idx).lead     = lv;
            cand(idx).Jt       = Jt;
            cand(idx).Ct       = Ct;
            cand(idx).feasible = feas.isFeasible;
            cand(idx).reason   = feas.reason;
            cand(idx).minGap   = min(lv.s - traj.s);
            cand(idx).maxAbsA  = max(abs(traj.a));
            cand(idx).maxAbsJ  = max(abs(traj.j));
        end
    end

    feasibleIdx = find([cand.feasible]);
    if isempty(feasibleIdx)
        error('No feasible candidate trajectory found. Relax constraints or expand candidate set.');
    end

    feasibleCosts = [cand(feasibleIdx).Ct];
    [~, localBest] = min(feasibleCosts);
    bestIdx = feasibleIdx(localBest);

    R.W          = W;
    R.candidates = cand;
    R.best       = cand(bestIdx);
end

function coeff = solveQuinticBoundaryValue(s0, v0, a0, s1, v1, a1, T)
% Solve quintic polynomial:
% s(t) = c0 + c1 t + c2 t^2 + c3 t^3 + c4 t^4 + c5 t^5
%
% Boundary conditions:
% s(0)=s0, sdot(0)=v0, sddot(0)=a0
% s(T)=s1, sdot(T)=v1, sddot(T)=a1

    c0 = s0;
    c1 = v0;
    c2 = a0 / 2.0;

    A = [ T^3,    T^4,     T^5;
          3*T^2,  4*T^3,   5*T^4;
          6*T,   12*T^2,  20*T^3 ];

    b = [ s1 - (c0 + c1*T + c2*T^2);
          v1 - (c1 + 2*c2*T);
          a1 - (2*c2) ];

    x = A \ b;

    c3 = x(1);
    c4 = x(2);
    c5 = x(3);

    coeff = [c0; c1; c2; c3; c4; c5];
end

function traj = sampleQuinticTrajectory(coeff, T, dt)
% Sample position, velocity, acceleration, jerk

    t = 0:dt:T;
    if abs(t(end)-T) > 1e-12
        t = [t, T];
    end

    c0 = coeff(1); c1 = coeff(2); c2 = coeff(3);
    c3 = coeff(4); c4 = coeff(5); c5 = coeff(6);

    s = c0 + c1*t + c2*t.^2 + c3*t.^3 + c4*t.^4 + c5*t.^5;
    v = c1 + 2*c2*t + 3*c3*t.^2 + 4*c4*t.^3 + 5*c5*t.^4;
    a = 2*c2 + 6*c3*t + 12*c4*t.^2 + 20*c5*t.^3;
    j = 6*c3 + 24*c4*t + 60*c5*t.^2;

    traj.t = t;
    traj.s = s;
    traj.v = v;
    traj.a = a;
    traj.j = j;
end

function lv = leadTrajectory(P, t)
% Lead vehicle motion with constant acceleration

    s = P.lv.s0 + P.lv.v0*t + 0.5*P.lv.a0*t.^2;
    v = P.lv.v0 + P.lv.a0*t;
    a = P.lv.a0 * ones(size(t));

    lv.t = t;
    lv.s = s;
    lv.v = v;
    lv.a = a;
end

function [sTarget, vTarget, aTarget] = targetTerminalState(P, T)
% s_target(t) = s_lv(t) - [D0 + tau * sdot_lv(t)]

    sLv = P.lv.s0 + P.lv.v0*T + 0.5*P.lv.a0*T^2;
    vLv = P.lv.v0 + P.lv.a0*T;
    aLv = P.lv.a0;

    sTarget = sLv - (P.follow.D0 + P.follow.tau * vLv);
    vTarget = vLv - P.follow.tau * aLv;
    aTarget = aLv;
end

function feas = checkFeasibility(traj, lv, P)
% Engineering feasibility checks:
% 1) nonnegative speed
% 2) speed upper bound
% 3) acceleration bounds
% 4) jerk bound
% 5) collision-free gap

    gap = lv.s - traj.s;

    cond1 = all(traj.v >= P.lim.vMin - 1e-9);
    cond2 = all(traj.v <= P.lim.vMax + 1e-9);
    cond3 = all(traj.a >= P.lim.aMin - 1e-9) && all(traj.a <= P.lim.aMax + 1e-9);
    cond4 = all(abs(traj.j) <= P.lim.jMax + 1e-9);
    cond5 = all(gap >= P.lim.minGap - 1e-9);

    feas.isFeasible = cond1 && cond2 && cond3 && cond4 && cond5;

    reasons = {};
    if ~cond1, reasons{end+1} = 'negative_speed'; end
    if ~cond2, reasons{end+1} = 'speed_upper_bound'; end
    if ~cond3, reasons{end+1} = 'acceleration_bound'; end
    if ~cond4, reasons{end+1} = 'jerk_bound'; end
    if ~cond5, reasons{end+1} = 'collision_gap'; end

    if isempty(reasons)
        feas.reason = 'ok';
    else
        feas.reason = strjoin(reasons, ', ');
    end
end

function T = buildSummaryTable(Results)
    n = numel(Results);

    CaseName = strings(n,1);
    kj       = zeros(n,1);
    kt       = zeros(n,1);
    ks       = zeros(n,1);
    Tstar    = zeros(n,1);
    dSstar   = zeros(n,1);
    Cost     = zeros(n,1);
    JerkCost = zeros(n,1);
    sEnd     = zeros(n,1);
    vEnd     = zeros(n,1);
    aEnd     = zeros(n,1);
    minGap   = zeros(n,1);
    maxAbsA  = zeros(n,1);
    maxAbsJ  = zeros(n,1);

    for i = 1:n
        b = Results(i).best;
        CaseName(i) = string(Results(i).W.name);
        kj(i)       = Results(i).W.kj;
        kt(i)       = Results(i).W.kt;
        ks(i)       = Results(i).W.ks;
        Tstar(i)    = b.T;
        dSstar(i)   = b.dS;
        Cost(i)     = b.Ct;
        JerkCost(i) = b.Jt;
        sEnd(i)     = b.traj.s(end);
        vEnd(i)     = b.traj.v(end);
        aEnd(i)     = b.traj.a(end);
        minGap(i)   = b.minGap;
        maxAbsA(i)  = b.maxAbsA;
        maxAbsJ(i)  = b.maxAbsJ;
    end

    T = table(CaseName, kj, kt, ks, Tstar, dSstar, Cost, JerkCost, ...
              sEnd, vEnd, aEnd, minGap, maxAbsA, maxAbsJ);
end

function plotCandidateBundles(P, Results)

    figure('Name','Candidate Bundles - Position','Color','w');
    tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    for i = 1:numel(Results)
        nexttile; hold on; grid on; box on;

        cand = Results(i).candidates;
        best = Results(i).best;

        % Plot feasible candidates lightly
        for k = 1:numel(cand)
            if cand(k).feasible
                plot(cand(k).traj.t, cand(k).traj.s, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
            end
        end

        % Lead vehicle and target for best horizon
        tBest = best.traj.t;
        lv = leadTrajectory(P, tBest);

        sTarget = zeros(size(tBest));
        for n = 1:numel(tBest)
            [sTarget(n),~,~] = targetTerminalState(P, tBest(n));
        end

        plot(tBest, lv.s, 'k--', 'LineWidth', 1.8, 'DisplayName', 'Lead Vehicle');
        plot(tBest, sTarget, 'm-.', 'LineWidth', 1.8, 'DisplayName', 'Target Following Position');
        plot(best.traj.t, best.traj.s, 'b', 'LineWidth', 2.5, 'DisplayName', 'Optimal');

        xlabel('Time [s]');
        ylabel('s(t) [m]');
        title(sprintf('%s | k_j=%.2f, k_t=%.2f, k_s=%.2f', ...
            Results(i).W.name, Results(i).W.kj, Results(i).W.kt, Results(i).W.ks));
        legend('Location','best');
    end
end

function plotOptimalComparisons(P, Results)

    figure('Name','Optimal Trajectory Comparison','Color','w');
    tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

    % -------- speed
    nexttile; hold on; grid on; box on;
    for i = 1:numel(Results)
        plot(Results(i).best.traj.t, Results(i).best.traj.v, 'LineWidth', 2.0, ...
            'DisplayName', Results(i).W.name);
    end
    yline(P.lv.v0, 'k--', 'Lead speed');
    xlabel('Time [s]');
    ylabel('v(t) [m/s]');
    title('Optimal Speed Profiles');
    legend('Location','best');

    % -------- acceleration
    nexttile; hold on; grid on; box on;
    for i = 1:numel(Results)
        plot(Results(i).best.traj.t, Results(i).best.traj.a, 'LineWidth', 2.0, ...
            'DisplayName', Results(i).W.name);
    end
    yline(P.lim.aMin, 'r--', 'a_{min}');
    yline(P.lim.aMax, 'r--', 'a_{max}');
    xlabel('Time [s]');
    ylabel('a(t) [m/s^2]');
    title('Optimal Acceleration Profiles');
    legend('Location','best');

    % -------- jerk
    nexttile; hold on; grid on; box on;
    for i = 1:numel(Results)
        plot(Results(i).best.traj.t, Results(i).best.traj.j, 'LineWidth', 2.0, ...
            'DisplayName', Results(i).W.name);
    end
    yline(P.lim.jMax, 'r--', 'j_{max}');
    yline(-P.lim.jMax, 'r--', '-j_{max}');
    xlabel('Time [s]');
    ylabel('j(t) [m/s^3]');
    title('Optimal Jerk Profiles');
    legend('Location','best');
end
Results = repmat(struct('field1',[],'field2',[]), 1, numElements);