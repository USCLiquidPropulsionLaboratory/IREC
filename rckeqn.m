function [ hb,ub,tb,h,t  ] = rckeqn( Mp,q1,q2 )
% solve the rocket equation with drag to obtain the burnout properties as
% well as the max altitude, assuming perfectly vertical flight

% simulation is currently 1-D, with positive being up

%% Extract Data

% extract data from array
mdot = q1(1);
g = q1(2);
Isp = q1(3);
d = q1(4);
Cd = q1(5);
R_sp = q1(6);
mu = q1(7);
dt = q1(8);

%% Compute Masses and Volumes

% calculate cross-sectional area of rocket
A = 0.25*pi*d^2;

[~,~,~,~,~,~,~,~,~,~,~,~,~,~,M0,Mb] = getMassAndVolume(q2,Mp);

%% Initialize Variables for Numerical Integration
mass = M0;      % get initial wet mass
t= 0;           % specify start time [sec]
u = 0;          % specify initial velocity (starting on launchpad) [m/s]
alt = 0;        % specify initial altitude (launchpad) [m]

%% Numerical Integration for Boost Phase

% continue iterations while there is fuel in the tank
% numerical integration uses differential form of the rocket equation with
% drag force D=0.5*rho*u^2*A*Cd
while mass > Mb
    grav = mu/(R_sp+alt)^2; % calculate local gravity
    rho = rho_alt(alt);     % calculate local atmospheric density
    
    % find change in velocity over timestep due to the effect of thrust,
    % drag, and gravity
    % Drag force is defined as 0.5*rho*A*Cd*u^2, where:
    %   rho = density of fluid
    %     A = cross sectional area (perpendicular to velocity vector)
    %    Cd = drag coefficient (pretty much a fudge factor, needs to be
    %         obtained empirically or using CFD - yay...)
    %     u = magnitude of velocity vector
    
    du = (g*Isp*mdot/mass - 0.5*rho*u^2*A*Cd/mass - grav)*dt;
    mass = mass - mdot*dt;  % get updated mass as propellant is ejected
    
    u = u + du;         % update velocity
    alt = alt + u*dt;   % update altitude
    t = t + dt;         % increment time step
end

hb = alt;   % get burnout altitude
ub = u;     % get burnout velocity
tb = t;     % get burnout time

%% Numerical Integration for Coast Phase

% continue iterations while velocity is positive (until apex)
while u > 0
    grav = mu/(R_sp+alt)^2; % calculate local gravity
    rho = rho_alt(alt);     % calculate local atmospheric density
    
    % find change in velocity over timestep due to the effect of thrust,
    % drag, and gravity
    du = -(0.5*rho*u^2*A*Cd/mass + grav)*dt;
    
    u = u + du;         % update velocity
    alt = alt + u*dt;   % update altitude  
    t = t + dt;         % increment time step
end

h = alt;    % get max altitude


end

