%% Double Lane Change - Reference Generation and Lateral Tracking
% Advanced engineering layout
%
% Model:
%   r_dot   = -(1/tau)*r + (vx/tau)*k_des
%   psi_dot = r
%   X_dot   = vx*cos(psi)
%   Y_dot   = vx*sin(psi)
%
% Reference path:
%   Y_ref = dy1/2*(1+tanh(z1)) - dy2/2*(1+tanh(z2))
%   z1 = alpha/dx1*(X-Xs1) - alpha/2
%   z2 = alpha/dx2*(X-Xs2) - alpha/2
%
% Controller:
%   k_des = kappa_ref + ky*e_y + kpsi*e_psi
%
% Notes:
% - The image likely contains a typo in the kinematic equations:
%   second X_dot should be Y_dot = vx*sin(psi)
% - Table values ds1, ds2 are interpreted as Xs1, Xs2

clear; clc; close all;

%% =========================
% 1) Parameters
% =========================
P = struct();

% Vehicle / model
P.model.vx   = 20.0;    % [m/s]
P.model.tau  = 0.25;    % [s] closed-loop yaw-rate time constant
P.model.kMax = 0.20;    % [1/m] curvature command saturation

% Simulation
P.sim.dt   = 0.01;      % [s]
P.sim.tEnd = 12.0;      % [s]

% Reference path parameters from table
P.path.dy1   = 3.5;     % [m]
P.path.dx1   = 25.0;    % [m]
P.path.Xs1   = 60.0;    % [m] interpreted from ds1
P.path.dy2   = 3.5;     % [m]
P.path.dx2   = 25.0;    % [m]
P.path.Xs2   = 130.0;   % [m] interpreted from ds2
P.path.alpha = 2.4;     % [-]

% Controller selection
% 'ff_fb'       : curvature feedforward + proportional feedback
% 'pure_pursuit': optional pure pursuit implementation
P.ctrl.mode = 'ff_fb';

% Feedback controller gains
P.ctrl.ky   = 0.015;    % [1/m^2]
P.ctrl.kpsi = 0.90;     % [1/m]

% Pure pursuit lookahead
P.ctrl.Ld = 15.0;       % [m]

% Initial condition
x0 = [0; 0; 0; 0];      % [X; Y; psi; r]

%% =========================
% 2) Run simulation
% =========================
R = simulateLaneChange(P, x0);

%% =========================
% 3) Performance summary
% =========================
fprintf('\n========== DLC TRACKING SUMMARY ==========\n');
fprintf('Controller mode        : %s\n', P.ctrl.mode);
fprintf('Final X [m]            : %.3f\n', R.X(end));
fprintf('Final Y [m]            : %.3f\n', R.Y(end));
fprintf('Max |e_y| [m]          : %.4f\n', max(abs(R.eY)));
fprintf('RMS e_y [m]            : %.4f\n', rmsLocal(R.eY));
fprintf('Max |e_psi| [deg]      : %.4f\n', rad2deg(max(abs(R.ePsi))));
fprintf('RMS e_psi [deg]        : %.4f\n', rad2deg(rmsLocal(R.ePsi)));
fprintf('Max |k_des| [1/m]      : %.4f\n', max(abs(R.kDes)));
fprintf('Max |r| [deg/s]        : %.4f\n', rad2deg(max(abs(R.r))));
fprintf('Final lateral error[m] : %.4f\n', R.eY(end));

%% =========================
% 4) Plots
% =========================
plotResults(R, P);

%% ============================================================
%                     LOCAL FUNCTIONS
% ============================================================

function R = simulateLaneChange(P, x0)

    dt   = P.sim.dt;
    tVec = 0:dt:P.sim.tEnd;
    N    = numel(tVec);

    % State vector = [X; Y; psi; r]
    x = zeros(4, N);
    x(:,1) = x0;

    % Logs
    Yref     = zeros(1, N);
    psiref   = zeros(1, N);
    kapparef = zeros(1, N);
    kDes     = zeros(1, N);
    eY       = zeros(1, N);
    ePsi     = zeros(1, N);

    for k = 1:N-1

        X   = x(1,k);
        Y   = x(2,k);
        psi = x(3,k);
        r   = x(4,k);

        % Reference evaluated at current X
        ref = referencePath(X, P.path);

        Yref(k)     = ref.Y;
        psiref(k)   = ref.psi;
        kapparef(k) = ref.kappa;

        % Tracking errors
        eY(k)   = ref.Y - Y;
        ePsi(k) = wrapAngle(ref.psi - psi);

        % Controller
        switch lower(P.ctrl.mode)

            case 'ff_fb'
                kCmd = ref.kappa + P.ctrl.ky * eY(k) + P.ctrl.kpsi * ePsi(k);

            case 'pure_pursuit'
                XLd   = X + P.ctrl.Ld;
                refLd = referencePath(XLd, P.path);
                alpha = wrapAngle(atan2(refLd.Y - Y, XLd - X) - psi);
                kCmd  = 2.0 * sin(alpha) / P.ctrl.Ld;

            otherwise
                error('Unknown controller mode: %s', P.ctrl.mode);
        end

        % Curvature saturation
        kCmd = max(min(kCmd, P.model.kMax), -P.model.kMax);
        kDes(k) = kCmd;

        % Integrate vehicle model with RK4
        x(:,k+1) = rk4Step(x(:,k), kCmd, dt, P.model);
    end

    % Final sample reference/log completion
    refEnd       = referencePath(x(1,end), P.path);
    Yref(end)     = refEnd.Y;
    psiref(end)   = refEnd.psi;
    kapparef(end) = refEnd.kappa;
    eY(end)       = Yref(end) - x(2,end);
    ePsi(end)     = wrapAngle(psiref(end) - x(3,end));
    kDes(end)     = kDes(end-1);

    R.t        = tVec;
    R.X        = x(1,:);
    R.Y        = x(2,:);
    R.psi      = x(3,:);
    R.r        = x(4,:);
    R.Yref     = Yref;
    R.psiref   = psiref;
    R.kapparef = kapparef;
    R.kDes     = kDes;
    R.eY       = eY;
    R.ePsi     = ePsi;
end

function xNext = rk4Step(x, kDes, dt, M)

    f = @(xx) vehicleDynamics(xx, kDes, M);

    k1 = f(x);
    k2 = f(x + 0.5*dt*k1);
    k3 = f(x + 0.5*dt*k2);
    k4 = f(x + dt*k3);

    xNext = x + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4);
end

function dx = vehicleDynamics(x, kDes, M)
    X   = x(1);
    Y   = x(2);
    psi = x(3);
    r   = x(4);

    %#ok<NASGU> % X,Y included for completeness
    dx = zeros(4,1);
    dx(1) = M.vx * cos(psi);                  % X_dot
    dx(2) = M.vx * sin(psi);                  % Y_dot
    dx(3) = r;                               % psi_dot
    dx(4) = -(1/M.tau)*r + (M.vx/M.tau)*kDes;% r_dot
end

function ref = referencePath(X, Pp)
% Generates Y_ref(X), psi_ref(X), kappa_ref(X)

    b1 = Pp.alpha / Pp.dx1;
    b2 = Pp.alpha / Pp.dx2;

    z1 = b1*(X - Pp.Xs1) - Pp.alpha/2;
    z2 = b2*(X - Pp.Xs2) - Pp.alpha/2;

    sech2_z1 = 1 ./ cosh(z1).^2;
    sech2_z2 = 1 ./ cosh(z2).^2;

    % Y_ref
    Y = Pp.dy1/2 * (1 + tanh(z1)) ...
      - Pp.dy2/2 * (1 + tanh(z2));

    % dY/dX
    dYdX = (Pp.dy1/2) * sech2_z1 * b1 ...
         - (Pp.dy2/2) * sech2_z2 * b2;

    % psi_ref = atan(dY/dX)
    psi = atan(dYdX);

    % d2Y/dX2
    d2YdX2 = -Pp.dy1 * b1^2 * sech2_z1 * tanh(z1) ...
             + Pp.dy2 * b2^2 * sech2_z2 * tanh(z2);

    % curvature
    kappa = d2YdX2 / (1 + dYdX^2)^(3/2);

    ref.Y     = Y;
    ref.psi   = psi;
    ref.kappa = kappa;
end

function ang = wrapAngle(ang)
    ang = mod(ang + pi, 2*pi) - pi;
end

function val = rmsLocal(x)
    val = sqrt(mean(x.^2));
end

function plotResults(R, P)

    figure('Name','Double Lane Change Tracking','Color','w');
    tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

    % 1) Path in XY plane
    nexttile; hold on; grid on; box on;
    plot(R.X, R.Yref, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Reference');
    plot(R.X, R.Y, 'b', 'LineWidth', 2.2, 'DisplayName', 'Vehicle');
    xlabel('X [m]');
    ylabel('Y [m]');
    title('Path Tracking in Cartesian Plane');
    legend('Location', 'best');

    % 2) Heading vs X
    nexttile; hold on; grid on; box on;
    plot(R.X, rad2deg(R.psiref), 'k--', 'LineWidth', 2.0, 'DisplayName', '\psi_{ref}');
    plot(R.X, rad2deg(R.psi), 'b', 'LineWidth', 2.2, 'DisplayName', '\psi');
    xlabel('X [m]');
    ylabel('\psi [deg]');
    title('Heading Tracking');
    legend('Location', 'best');

    % 3) Lateral error
    nexttile; hold on; grid on; box on;
    plot(R.X, R.eY, 'r', 'LineWidth', 2.0);
    xlabel('X [m]');
    ylabel('e_y [m]');
    title('Lateral Tracking Error');

    % 4) Heading error
    nexttile; hold on; grid on; box on;
    plot(R.X, rad2deg(R.ePsi), 'm', 'LineWidth', 2.0);
    xlabel('X [m]');
    ylabel('e_\psi [deg]');
    title('Heading Error');

    % 5) Curvature command
    nexttile; hold on; grid on; box on;
    plot(R.X, R.kapparef, 'k--', 'LineWidth', 2.0, 'DisplayName', '\kappa_{ref}');
    plot(R.X, R.kDes, 'b', 'LineWidth', 2.0, 'DisplayName', 'k_{des}');
    yline(P.model.kMax, 'r--', 'k_{max}');
    yline(-P.model.kMax, 'r--', '-k_{max}');
    xlabel('X [m]');
    ylabel('Curvature [1/m]');
    title('Reference vs Commanded Curvature');
    legend('Location', 'best');

    % 6) Yaw rate
    nexttile; hold on; grid on; box on;
    plot(R.X, rad2deg(R.r), 'b', 'LineWidth', 2.0);
    xlabel('X [m]');
    ylabel('r [deg/s]');
    title('Yaw Rate Response');
end