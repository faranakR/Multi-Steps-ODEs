%% Adaptive Time-Stepping Numerical Integration Solver
% This code implements adaptive numerical integration using Adams-Bashforth-Moulton 
% predictor-corrector methods for solving various types of differential equations.
%
% Author: Faranak Rajabi
% Institution: University of California, Santa Barbara
%
% Methods implemented:
%   - First-order AB-AM for initial step
%   - Second-order AB-AM for subsequent steps
%   - Adaptive time-stepping based on local error estimates

% %% Add methods directory to path
addpath(fullfile(pwd, 'methods'));

clear
close all
clc

%%% User Input and Problem Selection
% Available problems:
%   1. Linear ODE system: dy1/dt = -y1, dy2/dt = -10(y2-t^2) + 2t
%   2. Predator-Prey system: Lotka-Volterra equations
%   3. Van der Pol oscillator
%   4. Advection equation with spatial discretization
Q = input('Which problem do you want to solve? Choose between 1, 2, 3, 4: ');

%%% Parameters Initialization
h = 0.1;        % Initial step size
t = zeros(1);   % Time array initialization
Frac = 0.9;     % Safety fraction for step size prediction

%%% Problem-Specific Setup
if Q == 1 
    % Problem 1: Linear ODE System
    t_end = 1;                  % Integration interval
    ETol = input('Enter ETOL: '); % Error tolerance
    Y = [1; 2];                % Initial conditions [y1(0); y2(0)]
    % System: dy1/dt = -y1, dy2/dt = -10(y2-t^2) + 2t
    fun = @(Y,t) [-Y(1); -10*(Y(2)-t^2)+2*t];
    
elseif Q == 2
    % Problem 2: Predator-Prey System (Lotka-Volterra)
    t_end = 100;
    ETol = input('Enter ETOL: ');
    Y = [10; 10];              % Initial populations [prey; predator]
    % System: dy1/dt = 0.25y1 - 0.01y1y2 (prey)
    %         dy2/dt = -y2 + 0.01y1y2 (predator)
    fun = @(Y,t) [0.25*Y(1)-0.01*Y(1)*Y(2); -Y(2)+0.01*Y(1)*Y(2)];
    
elseif Q == 3
    % Problem 3: Van der Pol Oscillator
    t_end = 11;
    ETol = input('Enter ETOL: ');
    Y = [2; 0];                % Initial conditions [position; velocity]
    % System: dy1/dt = y2
    %         dy2/dt = 2(1-y1^2)y2 - y1
    fun = @(Y,t) [Y(2); 2*(1-Y(1)^2)*Y(2)-Y(1)];
    
elseif Q == 4
    % Problem 4: Advection Equation
    t_end = 1;
    ETol = input('Enter ETOL: ');
    toff = [0 0.25 0.5 0.6 0.8 1];  % Output time points
    L = 1;                          % Spatial domain length
    N = 100;                        % Number of spatial grid points
    x = linspace(0,L,N);           % Spatial grid
    Dx = x(2)-x(1);                % Grid spacing
    Y = exp(-10*x)';               % Initial condition
    % Construct advection matrix A for du/dt = Au
    A = diag([0 -1/Dx*ones(1,N-1)]) + diag(1/Dx*ones(1,N-1),-1);
    fun = @(Y,t) A*Y;
end

%%% Integration Loop Setup
Conv = 1;    % Convergence flag for time stepping
End = 0;     % Flag for reaching end of integration
i = 1;       % Time step counter

t(i+1) = t(i) + h(i);  % Initial time step

%%% Main Integration Loop
while 1
    % Check if integration interval is complete
    if t(i+1) > t_end && End == 0
        End = 1;
    end
    
    % Time stepping: First vs. subsequent steps
    if i == 1
        % First step: Use first-order AB-AM
        [Y(:,i+1), h(i+1), Conv] = AB_AM_1(Y(:,1:i), fun, h(i), t, ETol, Frac);
    else
        % Subsequent steps: Use second-order AB-AM
        [Y(:,i+1), h([i,i+1]), Conv] = AB_AM_2(Y(:,1:i), fun, h([i-1,i]), t, ETol, Frac);
    end
    
    % Handle non-convergent steps
    if Conv == 0
        h(i) = h(i+1);
        t(i+1) = t(i) + h(i);
        End = 0;
        continue;
    end
    
    % Exit if integration is complete
    if End == 1
        break
    end
    
    % Advance to next time step
    i = i + 1;
    t(i+1) = t(i) + h(i);
end

%%% Post-Processing for Problem 4
if Q == 4
    % Interpolate solution at specified output times
    Yoff = Y(:,1);
    for i = 2:length(toff)
        a = find(t <= toff(i));
        b = find(t > toff(i));
        if toff(i) == a(end)
            Yoff(:,i) = Y(:,a(end));
            continue
        end
        % Quadratic interpolation using solution and its derivatives
        Yoff(:,i) = Y(:,a(end)) + fun(Y(:,a(end)),t(a(end)))*(toff(i)-t(a(end))) + ...
            ((fun(Y(:,a(end)),t(a(end)))-fun(Y(:,a(end)),t(a(end))))/(2*h(a(end-1))))* ...
            (toff(i)-t(a(end)))^2;
    end
end

%%% Visualization
% Different plotting routines based on problem type
if Q == 1 || Q == 2 || Q == 3
    % Standard problems: Plot solution and step size history
    figure(1)
    plot(t,Y(1,:),t,Y(2,:),'--','linewidth',1.5);
    % ... [rest of plotting code]
else
    % Advection problem: Plot solution at specified times
    figure(1)
    plot(x,Yoff(:,1),x,Yoff(:,2),'^',x,Yoff(:,3),'--',x,Yoff(:,4),'+',...
         x,Yoff(:,5),'*',x,Yoff(:,6),'s','linewidth',1.5);
    % ... [rest of plotting code]
end