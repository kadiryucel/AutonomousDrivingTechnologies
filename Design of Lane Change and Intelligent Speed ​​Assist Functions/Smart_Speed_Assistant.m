%% Intelligent Speed Assistant (ISA) - Advanced Engineering Layout
% Based on the provided equations:
%
% Route segment propagation:
%   s_{i+1}   = s_i + L_i
%   psi_{i+1} = psi_i + L_i*kappa_i
%   X_{i+1}   = X_i + L_i*cos(psi_i + kappa_i*L_i/2)
%   Y_{i+1}   = Y_i + L_i*sin(psi_i + kappa_i*L_i/2)
%
% Road speed limit due to curvature:
%   v_lim_road = sqrt(a_y_comfort / kappa)
%
% Longitudinal motion:
%   s_dot  = v_x
%   v_dot  = a_x_des
%
% Trigger distance:
%   d_trig = (v_curr^2 - v_lim^2) / (2*a_x_comfort)
%
% Controller logic:
% - Combined speed limit: v_lim = min(v_lim_traffic, v_lim_road)
% - Preview the next lower speed limit zone
% - If distance to next lower zone <= d_trig, command comfortable braking
% - If speed limit increases, accelerate gently back to allowed speed

clear; clc; close all;

%% =========================
% 1) Parameters
% =========================
P = struct();

% Route discretization
P.route.sEnd = 3000.0;    % [m]
P.route.ds   = 1.0;       % [m] segment length L_i

% Comfort parameters
P.comfort.ay = 2.0;       % [m/s^2] comfortable lateral acceleration
P.comfort.axDecel = 2.0;  % [m/s^2] comfortable braking magnitude
P.comfort.axAccel = 1.0;  % [m/s^2] comfortable acceleration magnitude

% Longitudinal controller
P.ctrl.kv = 0.6;          % [1/s] proportional speed recovery gain
P.ctrl.speedTol = 0.10;   % [m/s] tolerance for mode switching

% Simulation
P.sim.dt   = 0.05;        % [s]
P.sim.tEnd = 250.0;       % [s] simulation cap

% Initial condition
P.init.s0  = 0.0;         % [m]
P.init.v0  = 90/3.6;      % [m/s] start at 90 km/h

%% =========================
% 2) Build route profiles
% =========================
Route = buildRouteProfiles(P);
Route = integrateRouteGeometry(Route);

%% =========================
% 3) Run ISA simulation
% =========================
R = simulateISA(P, Route);

%% =========================
% 4) Summary
% =========================
printSummary(R, Route);

%% =========================
% 5) Plots
% =========================
plotISAResults(Route, R, P);

%% ============================================================
%                        LOCAL FUNCTIONS
% ============================================================

function Route = buildRouteProfiles(P)
% Build s-grid, curvature profile, traffic speed limits,
% road speed limits, and combined speed limits.

    s = 0:P.route.ds:P.route.sEnd;
    N = numel(s);

    kappa = zeros(1, N);          % [1/m]
    vTraffic_kmh = zeros(1, N);   % [km/h]

    for i = 1:N
        si = s(i);

        % Curvature profile from the statement
        if si <= 400
            kappa(i) = 0.0;
        elseif si <= 1000
            kappa(i) = 0.005;
        elseif si <= 1800
            kappa(i) = 0.0;
        elseif si <= 2400
            kappa(i) = 0.0025;
        else
            kappa(i) = 0.0;
        end

        % Traffic speed limit profile from the statement
        if si <= 400
            vTraffic_kmh(i) = 90;
        elseif si <= 1000
            vTraffic_kmh(i) = 50;
        elseif si <= 1800
            vTraffic_kmh(i) = 70;
        elseif si <= 2400
            vTraffic_kmh(i) = 50;
        else
            vTraffic_kmh(i) = 90;
        end
    end

    % Convert traffic speed limit to m/s
    vTraffic = vTraffic_kmh / 3.6;

    % Road speed limit due to curvature
    % v_lim_road = sqrt(a_y_comfort / kappa)
    vRoad = inf(1, N);
    idxCurve = abs(kappa) > 1e-12;
    vRoad(idxCurve) = sqrt(P.comfort.ay ./ abs(kappa(idxCurve)));

    % Combined limit
    vLim = min(vTraffic, vRoad);

    Route.s = s;
    Route.kappa = kappa;
    Route.vTraffic = vTraffic;
    Route.vTraffic_kmh = vTraffic_kmh;
    Route.vRoad = vRoad;
    Route.vRoad_kmh = 3.6 * vRoad;
    Route.vLim = vLim;
    Route.vLim_kmh = 3.6 * vLim;
end

function Route = integrateRouteGeometry(Route)
% Integrate route centerline in Cartesian coordinates using:
%   psi_{i+1} = psi_i + L_i*kappa_i
%   X_{i+1}   = X_i + L_i*cos(psi_i + kappa_i*L_i/2)
%   Y_{i+1}   = Y_i + L_i*sin(psi_i + kappa_i*L_i/2)

    s = Route.s;
    ds = s(2) - s(1);
    N = numel(s);

    X = zeros(1, N);
    Y = zeros(1, N);
    psi = zeros(1, N);

    for i = 1:N-1
        ki = Route.kappa(i);

        psi(i+1) = psi(i) + ds*ki;
        X(i+1) = X(i) + ds*cos(psi(i) + ki*ds/2);
        Y(i+1) = Y(i) + ds*sin(psi(i) + ki*ds/2);
    end

    Route.X = X;
    Route.Y = Y;
    Route.psi = psi;
end

function R = simulateISA(P, Route)
% Preview-based Intelligent Speed Assistant.
%
% Logic:
% 1) If current speed exceeds local speed limit -> brake comfortably now
% 2) Else preview the next lower limit zone
% 3) Compute d_trig = (v_curr^2 - v_target^2)/(2*a_x_comfort)
% 4) If distance to lower zone <= d_trig -> brake comfortably
% 5) Else accelerate gently toward current local limit

    dt   = P.sim.dt;
    tEnd = P.sim.tEnd;

    maxSteps = floor(tEnd/dt) + 1;

    tLog      = zeros(1, maxSteps);
    sLog      = zeros(1, maxSteps);
    vLog      = zeros(1, maxSteps);
    aLog      = zeros(1, maxSteps);
    vLocLog   = zeros(1, maxSteps);
    vTarLog   = nan(1, maxSteps);
    dTrigLog  = nan(1, maxSteps);
    dNextLog  = nan(1, maxSteps);
    modeLog   = strings(1, maxSteps);

    s = P.init.s0;
    v = P.init.v0;

    k = 1;
    tLog(k) = 0.0;
    sLog(k) = s;
    vLog(k) = v;

    while s < Route.s(end) && k < maxSteps

        % Local speed limit at current position
        vLocal = getProfileAtS(Route.s, Route.vLim, s);
        vLocLog(k) = vLocal;

        % Determine control action
        [sNextLower, vNextLower, foundLower] = findNextLowerLimitZone(Route, s, P.ctrl.speedTol);

        if v > vLocal + P.ctrl.speedTol
            % Immediate local overspeed correction
            axDes = -P.comfort.axDecel;
            mode = "local_limit_brake";
            dTrig = max((v^2 - vLocal^2)/(2*P.comfort.axDecel), 0);
            dNext = 0;
            vTarget = vLocal;

        elseif foundLower
            dNext = sNextLower - s;
            dTrig = max((v^2 - vNextLower^2)/(2*P.comfort.axDecel), 0);
            vTarget = vNextLower;

            if dNext <= dTrig && v > vNextLower + P.ctrl.speedTol
                axDes = -P.comfort.axDecel;
                mode = "preview_brake";
            else
                if v < vLocal - P.ctrl.speedTol
                    axDes = min(P.ctrl.kv*(vLocal - v), P.comfort.axAccel);
                    mode = "speed_recovery";
                else
                    axDes = 0.0;
                    mode = "hold";
                end
            end

        else
            % No more restrictive limit ahead
            dNext = nan;
            dTrig = nan;
            vTarget = nan;

            if v < vLocal - P.ctrl.speedTol
                axDes = min(P.ctrl.kv*(vLocal - v), P.comfort.axAccel);
                mode = "speed_recovery";
            else
                axDes = 0.0;
                mode = "hold";
            end
        end

        % Integrate longitudinal dynamics
        % s_dot = v
        % v_dot = a_x_des
        vNext = max(v + axDes*dt, 0.0);
        sNext = s + v*dt + 0.5*axDes*dt^2;
        sNext = min(sNext, Route.s(end));

        % Advance
        aLog(k)     = axDes;
        vTarLog(k)  = vTarget;
        dTrigLog(k) = dTrig;
        dNextLog(k) = dNext;
        modeLog(k)  = mode;

        s = sNext;
        v = vNext;

        k = k + 1;
        tLog(k) = tLog(k-1) + dt;
        sLog(k) = s;
        vLog(k) = v;
    end

    % Final sample completion
    vLocLog(k) = getProfileAtS(Route.s, Route.vLim, sLog(k));
    aLog(k) = aLog(max(k-1,1));
    vTarLog(k) = vTarLog(max(k-1,1));
    dTrigLog(k) = dTrigLog(max(k-1,1));
    dNextLog(k) = dNextLog(max(k-1,1));
    modeLog(k) = modeLog(max(k-1,1));

    % Trim logs
    R.t      = tLog(1:k);
    R.s      = sLog(1:k);
    R.v      = vLog(1:k);
    R.a      = aLog(1:k);
    R.vLocal = vLocLog(1:k);
    R.vTarget = vTarLog(1:k);
    R.dTrig  = dTrigLog(1:k);
    R.dNext  = dNextLog(1:k);
    R.mode   = modeLog(1:k);

    % Derived quantities
    R.v_kmh      = 3.6*R.v;
    R.vLocal_kmh = 3.6*R.vLocal;
    R.ayActual   = R.v.^2 .* getProfileAtSVector(Route.s, Route.kappa, R.s);
end

function [sNext, vNext, found] = findNextLowerLimitZone(Route, sCurr, tol)
% Find the first future location where the combined speed limit
% becomes lower than the current local limit.

    idxCurr = getIndexAtS(Route.s, sCurr);
    vCurrLim = Route.vLim(idxCurr);

    found = false;
    sNext = nan;
    vNext = nan;

    for j = idxCurr+1:numel(Route.s)
        if Route.vLim(j) < vCurrLim - tol
            found = true;
            sNext = Route.s(j);
            vNext = Route.vLim(j);
            return;
        end
    end
end

function idx = getIndexAtS(sGrid, sQuery)
    ds = sGrid(2) - sGrid(1);
    idx = floor(sQuery/ds) + 1;
    idx = max(1, min(idx, numel(sGrid)));
end

function val = getProfileAtS(sGrid, profile, sQuery)
    idx = getIndexAtS(sGrid, sQuery);
    val = profile(idx);
end

function vals = getProfileAtSVector(sGrid, profile, sQueryVec)
    vals = zeros(size(sQueryVec));
    for i = 1:numel(sQueryVec)
        vals(i) = getProfileAtS(sGrid, profile, sQueryVec(i));
    end
end

function printSummary(R, Route)

    speedError = R.v - R.vLocal;
    maxOverspeed = max(speedError);
    maxOverspeed = max(maxOverspeed, 0);

    travelTime = R.t(end);
    finalSpeed = R.v(end);
    maxBrake   = min(R.a);
    maxAccel   = max(R.a);
    maxAy      = max(abs(R.ayActual));

    % Transition points in the combined speed limit
    transitions = find(abs(diff(Route.vLim)) > 1e-9) + 1;

    fprintf('\n========== INTELLIGENT SPEED ASSISTANT SUMMARY ==========\n');
    fprintf('Route length [m]                 : %.1f\n', Route.s(end));
    fprintf('Travel time [s]                  : %.2f\n', travelTime);
    fprintf('Final speed [km/h]               : %.2f\n', 3.6*finalSpeed);
    fprintf('Max overspeed wrt local limit[km/h]: %.3f\n', 3.6*maxOverspeed);
    fprintf('Max commanded accel [m/s^2]      : %.3f\n', maxAccel);
    fprintf('Max commanded brake [m/s^2]      : %.3f\n', maxBrake);
    fprintf('Max lateral accel [m/s^2]        : %.3f\n', maxAy);
    fprintf('Comfort lateral accel limit[m/s^2]: %.3f\n', 2.0);
    fprintf('Number of limit transitions      : %d\n', numel(transitions));

    fprintf('\nSpeed limit transitions:\n');
    for i = 1:numel(transitions)
        idx = transitions(i);
        fprintf('  s = %6.1f m | new combined limit = %5.1f km/h\n', ...
            Route.s(idx), Route.vLim_kmh(idx));
    end
end

function plotISAResults(Route, R, P)

    figure('Name','Intelligent Speed Assistant','Color','w');
    tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

    % 1) Route geometry
    nexttile; hold on; grid on; box on;
    plot(Route.X, Route.Y, 'b', 'LineWidth', 2.0);
    xlabel('X [m]');
    ylabel('Y [m]');
    title('Route Geometry in Cartesian Plane');
    axis equal;

    % 2) Curvature profile
    nexttile; hold on; grid on; box on;
    plot(Route.s, Route.kappa, 'k', 'LineWidth', 2.0);
    xlabel('s [m]');
    ylabel('\kappa [1/m]');
    title('Route Curvature Profile');

    % 3) Speed limits along route
    nexttile; hold on; grid on; box on;
    plot(Route.s, Route.vTraffic_kmh, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Traffic Limit');
    plot(Route.s, Route.vRoad_kmh, 'm-.', 'LineWidth', 2.0, 'DisplayName', 'Road Limit');
    plot(Route.s, Route.vLim_kmh, 'b', 'LineWidth', 2.5, 'DisplayName', 'Combined Limit');
    xlabel('s [m]');
    ylabel('Speed [km/h]');
    title('Traffic / Road / Combined Speed Limits');
    legend('Location','best');

    % 4) Vehicle speed vs route position
    nexttile; hold on; grid on; box on;
    plot(Route.s, Route.vLim_kmh, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Combined Limit');
    plot(R.s, R.v_kmh, 'r', 'LineWidth', 2.2, 'DisplayName', 'Vehicle Speed');
    xlabel('s [m]');
    ylabel('Speed [km/h]');
    title('ISA Speed Tracking Along Route');
    legend('Location','best');

    % 5) Speed and acceleration vs time
    nexttile; yyaxis left; hold on; grid on; box on;
    plot(R.t, R.v_kmh, 'b', 'LineWidth', 2.0, 'DisplayName', 'Vehicle Speed');
    ylabel('Speed [km/h]');
    yyaxis right; hold on;
    plot(R.t, R.a, 'r', 'LineWidth', 1.8, 'DisplayName', 'a_{x,des}');
    ylabel('a_{x,des} [m/s^2]');
    xlabel('Time [s]');
    title('Longitudinal Response vs Time');

    % 6) Preview logic diagnostics
    nexttile; hold on; grid on; box on;
    plot(R.s, R.dTrig, 'b', 'LineWidth', 2.0, 'DisplayName', 'd_{trig}');
    plot(R.s, R.dNext, 'r--', 'LineWidth', 2.0, 'DisplayName', 'Distance to Next Lower Limit');
    xlabel('s [m]');
    ylabel('Distance [m]');
    title('Preview Braking Logic');
    legend('Location','best');

    sgtitle(sprintf('ISA | a_y^{comfort}=%.1f m/s^2, a_x^{decel}=%.1f m/s^2', ...
        P.comfort.ay, P.comfort.axDecel));
end