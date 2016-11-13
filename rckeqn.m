function [ hb,ub,tb,h,t  ] = rckeqn( Mp,mdot,Ms_0,Ml,g,Isp,alpha,A,Cd )
% solve the rocket equation with drag to obtain the burnout properties as
% well as the max altitude, assuming perfectly vertical flight

% Assume structural mass of rocket is the base structural mass, plus an
% overhead as a fraction of propellant mass (accounts for oxidizer, fuel
% and ox tankage, and fuselage)
Ms = Ms_0 + alpha*Mp;

% calculate wet and dry masses of rocket
M0 = Ms + Mp + Ml;
Mb = Ms + Ml;

%% Initialize Variables for Numerical Integration
mass_prop = Mp; % get initial propellant mass
mass = M0;      % get initial wet mass
dt = 0.01;      % specify timestep [sec]
t= 0;           % specify start time [sec]
u = 0;          % specify initial velocity (starting on launchpad) [m/s]
alt = 0;        % specify initial altitude (launchpad) [m]

%% Numerical Integration for Boost Phase

% continue iterations while there is fuel in the tank
% numerical integration uses differential form of the rocket equation with
% drag force D=0.5*rho*u^2*A*Cd
while mass_prop > 0
    mass = mass - mdot*dt;
    du = (g*Isp*mdot/mass - 0.5*rho_alt(alt)*u^2*A*Cd/mass - g)*dt;
    u = u + du;
    alt = alt + u*dt;
    mass_prop = mass - Mb;
    t = t + dt;
end

hb = alt;   % get burnout altitude
ub = u;     % get burnout velocity
tb = t;     % get burnout time

%% Numerical Integration for Coast Phase

% continue iterations while velocity is positive (until apex)
while u > 0
    du = -(0.5*rho_alt(alt)*u^2*A*Cd/mass + g)*dt;
    u = u + du;
    alt = alt + u*dt;  
    t = t + dt;
end

h = alt;    % get max altitude


end

