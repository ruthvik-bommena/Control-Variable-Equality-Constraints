function MinimumTime_PosDepVelMag
%%  < File Description >
%    Author:        Ruthvik Bommena
%    File Name:     MinimumTime_PosDepVelMag.m
%    Compiler:      MATLAB R2022b
%    Date:          12 February, 2023
%    Affiliation:   Department of Aerospace Engineering, University of Illinois Urbana-Champaign.
%    Description:   Function to solve minimum-time paths through a three-dimensional 
%                   region of position dependent velocity magnitude.
%    References:    Ch 3. Sec 3.2, Applied Optimal Control, 1975, A.E. Bryson. Jr, Yu-Chi Ho

clear; close all; clc;

% [x; y; z]
x0 = [12000;5000;7500]; % m
xf = [0;0;0];

t0 = 0;

V = [50;10;10]; % aircraft velocity vector m/s (assuming constant velocity)

opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15); % ode
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,...
    'MaxIter',1e3,'TolFun',1e-14,'TolX',1e-16,...
    'UseParallel',true); % fsolve

%% Numerical Solution 
lam0_guess = [0.0203425681285667; 0.0423803502678474; 0.0635705254017710; 932.791509395320];

iter = 1; err = 1; errTol = 1e-10; iterMax = 100;
   while err > errTol && iter < iterMax
        [lam0,~] = fsolve(@costFunction,lam0_guess,options,x0,xf,opts_ode,V,t0);
        err = norm(costFunction(lam0_guess,x0,xf,opts_ode,V,t0));
        lam0_guess = lam0;
        iter = iter+1; 
    end

% Trajectory
tf = lam0(4);
[t_minT,X_minT] = ode45(@ThreeDimAircraft,[t0 tf],[x0; lam0(1:3)],opts_ode,V);  % ode propagator

% Plots
PlotSolution(X_minT,t_minT,x0,xf); % plots

end


%% Functions
function PlotSolution(X_minT,t_minT,x0,xf)

    for ii = 1:length(t_minT)
        rel_Pos(ii) = norm(X_minT(ii,1:3)-xf');
    end
    
    % Plot
    figure; grid on; hold on;
    plot3(X_minT(:,1),X_minT(:,2),X_minT(:,3));
    xlabel('Distance (m)'); ylabel('Distance (m)'); zlabel('Distance (m)'); 
    title('Problem 1 - Aircraft Trajectory',['x0 = ',num2str(x0(1:3)'),'  |  ','xf = ',num2str(xf(1:3)')]); 
    view(3);

    figure; grid on; hold on; 
    plot(t_minT/60,rel_Pos);
    xlabel('Time (min)'); ylabel('Distance (m)');
    title('Problem 1 - Relative Position of Aircraft',['x0 = ',num2str(x0(1:3)'),'  |  ','xf = ',num2str(xf(1:3)')]);

end 


function Xdot = ThreeDimAircraft(~,X,V)

lambda = X(4:6);

aircraftAxis = -lambda/norm(lambda); % unit vector

% State Dynamics
xDot = [V(1)*aircraftAxis(1);...
        V(2)*aircraftAxis(2);...
        V(3)*aircraftAxis(3)];

% Costate Dynamics
lambdaDot = [0;...
             0;...
             0];

Xdot = [xDot; lambdaDot];

end


function H = hamiltonian(~,X,V)

lambda = X(:,4:6)';

magLambda = vecnorm(lambda);

    for ii = 1:length(X)
        
        aircraftAxis = -lambda(:,ii)./magLambda(ii);
        
        % State Differential Equations
       xDot = [V(1)*aircraftAxis(1);...
               V(2)*aircraftAxis(2);...
               V(3)*aircraftAxis(3)];
        
        % Hamiltonian
        H(ii) = (lambda(:,ii)' * xDot);
        
    end

end


function err = costFunction(lam0_guess,x0,xf,opts_ode,V,t0)

lam0 = lam0_guess(1:3); 
tf = lam0_guess(4); 

[t,X] = ode45(@ThreeDimAircraft,[t0 tf],[x0; lam0],opts_ode,V);

H = hamiltonian(t,X,V);
H = H'; 

err =  [X(end,1:3)'; H(end)+1] - [xf; 0]; % transversality condition for minT (H_f = 0)

end