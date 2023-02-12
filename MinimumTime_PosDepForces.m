function MinimumTime_PosDepForces
%%  < File Description >
%    Author:        Ruthvik Bommena
%    File Name:     MinimumTime_PosDepForces.m
%    Compiler:      MATLAB R2022b
%    Date:          12 February, 2023
%    Affiliation:   Department of Aerospace Engineering, University of Illinois Urbana-Champaign.
%    Description:   Function to solve minimum-time paths through a three-dimensional 
%                   region of position dependent forces.
%    References:    Ch 3. Applied Optimal Control, 1975, A.E. Bryson. Jr, Yu-Chi Ho

clear all; clc;

% [x; y; z; vx; vy; vz]
x0 = [12000;5000;7500;50;10;10]; % m & m/s
xf = [0;0;0;0;0;0];

% Aircraft Information - Piper Archer PA28A
m = 1120; % kg
Cd = 0.0296; % coefficient of drag
uMax = 5; % m/s^2
area = 20; % m^2 
density = 1.225; % air density at sea level kg/m^3

t0 = 0;
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,...
    'MaxIter',1e3,'TolFun',1e-14,'TolX',1e-16,...
    'UseParallel',true); % fsolve

%% Numerical Solution 
    if nargin < 6
       lam0_guess = rand(7,1);
    end 

rho = .5; rho_f = 1e-5; iter = 1;
    while rho > rho_f
        [lam0,fval] = fsolve(@costFunction,lam0_guess,options,x0,xf,opts_ode,t0,Cd,uMax,density,area,rho,m);
        rho = rho*0.5;
        lam0_guess = lam0;
        iter = iter + 1;
    end
tf = lam0(7);
[t_minT,X_minT] = ode45(@ThreeDimAircraft,[t0 tf],[x0; lam0(1:6)],opts_ode,Cd,uMax,density,area,rho,m);  % ode propagator

% Plots
PlotSolution(X_minT,t_minT,x0,xf); % plots

end


%% Functions
function hf = PlotSolution(X_minT,t_minT,x0,xf)

    for ii = 1:length(t_minT)
        rel_Pos(ii) = norm(X_minT(ii,1:3)-xf(1:3)');
        rel_Vel(ii) = norm(X_minT(ii,4:6)-xf(4:6)');
    end
    
% Plot
figure; grid on; hold on;
plot3(X_minT(:,1),X_minT(:,2),X_minT(:,3));
xlabel('Distance (m)'); ylabel('Distance (m)'); zlabel('Distance (m)'); 
title('Problem 2 - Aircraft Trajectory',['x0 = ',num2str(x0(1:3)'),'  |  ','xf = ',num2str(xf(1:3)')]); 
view(3);

figure; grid on; hold on; 
plot(t_minT/60,rel_Pos);
xlabel('Time (min)'); ylabel('Distance (m)');
title('Problem 2 - Relative Position of Aircraft',['x0 = ',num2str(x0(1:3)'),'  |  ','xf = ',num2str(xf(1:3)')]);

end 


function Xdot = ThreeDimAircraft(t,X,Cd,uMax,density,area,rho,m)

x = X(1:6);
lambda = X(7:12);

S = norm(lambda(4:6));
delta = 0.5*(1+tanh(S/rho));

aircraftAxis = -lambda/norm(lambda); % unit vector

% State Dynamics
xDot = [x(4);x(5);x(6);...
        delta*uMax*aircraftAxis(1) - (0.5*density*(x(4)^2)*Cd*area)/m;...
        delta*uMax*aircraftAxis(2) - (0.5*density*(x(5)^2)*Cd*area)/m;...
        delta*uMax*aircraftAxis(3) - (0.5*density*(x(6)^2)*Cd*area)/m];

% Costate Dynamics
lambdaDot = [0;0;0;...
             0.5*density*(2*x(4))*Cd*area;...
             0.5*density*(2*x(5))*Cd*area;...
             0.5*density*(2*x(6))*Cd*area];

Xdot = [xDot; lambdaDot];

end


function H = hamiltonian(t,X,Cd,uMax,density,area,rho,m)

x   = X(:,1:6);
lambda = X(:,7:12)';

S = vecnorm(lambda(4:6,:));
delta = 0.5*(1+tanh(S/rho));

magLambda = vecnorm(lambda);

    for ii = 1:length(X)
        aircraftAxis = -lambda(:,ii)./magLambda(ii);
        % State Differential Equations
        xDot = [x(4);x(5);x(6);...
                delta(ii)*uMax*aircraftAxis(1) - (0.5*density*(x(ii,4)^2)*Cd*area)/m;...
                delta(ii)*uMax*aircraftAxis(2) - (0.5*density*(x(ii,5)^2)*Cd*area)/m;...
                delta(ii)*uMax*aircraftAxis(3) - (0.5*density*(x(ii,6)^2)*Cd*area)/m];
        % Hamiltonian
        H(ii) = (lambda(:,ii)' * xDot); 
    end

end


function err = costFunction(lam0_guess,x0,xf,opts_ode,t0,Cd,uMax,density,area,rho,m)

lam0 = lam0_guess(1:6); 
tf = lam0_guess(7); 

[t,X] = ode45(@ThreeDimAircraft,[t0 tf],[x0; lam0],opts_ode,Cd,uMax,density,area,rho,m);

H = hamiltonian(t,X,Cd,uMax,density,area,rho,m);
H = H'; 

err =  [X(end,1:6)'; H(end)+1] - [xf; 0]; % transversality condition for minT (H_f = 0)

end