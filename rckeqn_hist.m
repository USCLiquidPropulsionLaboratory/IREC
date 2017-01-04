function [ h,u,T,t ] = rckeqn_hist( Mp,q1,q2 )
% solve the rocket equation with drag; outputs time history of trajectory
% data for plotting purposes

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
n = 1;          % initialize index counter
mass = M0;      % get initial wet mass
u(1) = 0;       % specify initial velocity (starting on launchpad) [m/s]
h(1) = 0;       % specify initial altitude (launchpad) [m]
T(1) = T_alt(h(1)); % Air temperature at altitude (launchpad) [K]
t(1) = 0;       % time at initial altitude (launchpad) [sec]

%% Numerical Integration for Boost Phase

% continue iterations while there is fuel in the tank
% numerical integration uses differential form of the rocket equation with
% drag force D=0.5*rho*u^2*A*Cd
while mass > Mb
    n = n + 1;              % increment counter
    grav = mu/(R_sp+h(n-1))^2; % calculate local gravity
    rho = rho_alt(h(n-1));     % calculate local atmospheric density
    
    % find change in velocity over timestep due to the effect of thrust,
    % drag, and gravity
    % Drag force is defined as 0.5*rho*A*Cd*u^2, where:
    %   rho = density of fluid
    %     A = cross sectional area (perpendicular to velocity vector)
    %    Cd = drag coefficient (pretty much a fudge factor, needs to be
    %         obtained empirically or using CFD - yay...)
    %     u = magnitude of velocity vector
    
    du = (g*Isp*mdot/mass - 0.5*rho*u(n-1)^2*A*Cd/mass - grav)*dt;
    mass = mass - mdot*dt;  % get updated mass as propellant is ejected
   
    % append data to arrays for output
    u(n) = u(n-1) + du;
    h(n) = h(n-1) + u(n)*dt;
    T(n) = T_alt(h(n));
    t(n) = t(n-1) + dt;
end


%% Numerical Integration for Coast Phase

% continue iterations while velocity is positive (until apex)
while u(n) > 0
    n = n + 1;              % increment counter
    grav = mu/(R_sp+h(n-1))^2; % calculate local gravity
    rho = rho_alt(h(n-1));     % calculate local atmospheric density
    
    % find change in velocity over timestep due to the effect of thrust,
    % drag, and gravity
    du = -(0.5*rho*u(n-1)^2*A*Cd/mass + grav)*dt;
    
    % append data to arrays for output
    u(n) = u(n-1) + du;
    h(n) = h(n-1) + u(n)*dt;  
    T(n) = T_alt(h(n));
    t(n) = t(n-1) + dt;
end

% if iterations go (slightly) past apogee, delete last station value
if u(end) < 0
    u = u(1:end-1);
    h = h(1:end-1);
    T = T(1:end-1);
    t = t(1:end-1);
end


end

