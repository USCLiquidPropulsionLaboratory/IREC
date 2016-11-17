function [ hb,ub,tb,h,t  ] = rckeqn( Mp,q )
% solve the rocket equation with drag to obtain the burnout properties as
% well as the max altitude, assuming perfectly vertical flight

% extract data from array
mdot = q(1);
Ms_0 = q(2);
Ml = q(3);
g = q(4);
Isp = q(5);
alpha = q(6);
d = q(7);
Cd = q(8);
of_ratio = q(9);
R_sp = q(10);
mu = q(11);
rho_f = q(12);
rho_ox = q(13);
rho_tank = q(14);
sig_tank = q(15);
FS_tank = q(16);
P_f = q(17);
P_ox = q(18);
d_tank = q(19);

% calculate cross-sectional area of rocket
A = 0.25*pi*d^2;

Mf = Mp/(1+of_ratio);   % fuel mass
Mox = Mp - Mf;          % oxidizer mass

Vf = Mf/rho_f;      % volume of propellant
Vox = Mox/rho_ox;   % volume of oxidizer

% we are using cylindrical tanks with domed ends, so lets find the height
% of the cylindrical portion of the tanks
W_f = (Vf - pi*d_tank^3/6)/(pi*d_tank^2/4);  
W_ox = (Vox - pi*d_tank^3/6)/(pi*d_tank^2/4);

% determine total heigh of tanks
h_tank_f = W_f + d_tank;
h_tank_ox = W_ox + d_tank;

% determine masses of tanks, assuming a safety factor on yield stress
M_tank_f = 0.5*pi*d_tank^2/2*(d_tank/2+W_f)*FS_tank*P_f*rho_tank/sig_tank;
M_tank_ox = 0.5*pi*d_tank^2/2*(d_tank/2+W_ox)*FS_tank*P_ox*rho_tank/sig_tank;

% Assume structural mass of rocket is the base structural mass, plus
% tankage masses, and an overhead for the rest of the fuselage
Ms = Ms_0 + M_tank_f + M_tank_ox + alpha*Mp;
% Ms = Ms_0 + alpha*Mp*of_ratio;

% calculate wet and dry masses of rocket
M0 = Ms + Mp + Ml;
Mb = Ms + Ml;

%% Initialize Variables for Numerical Integration
mass = M0;      % get initial wet mass
dt = 0.01;      % specify timestep [sec]
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
    
    mass = mass - mdot*dt;
    du = (g*Isp*mdot/mass - 0.5*rho*u^2*A*Cd/mass - grav)*dt;
    u = u + du;
    alt = alt + u*dt;
    t = t + dt;
end

hb = alt;   % get burnout altitude
ub = u;     % get burnout velocity
tb = t;     % get burnout time

%% Numerical Integration for Coast Phase

% continue iterations while velocity is positive (until apex)
while u > 0
    grav = mu/(R_sp+alt)^2; % calculate local gravity
    rho = rho_alt(alt);     % calculate local atmospheric density
    
    du = -(0.5*rho*u^2*A*Cd/mass + grav)*dt;
    u = u + du;
    alt = alt + u*dt;  
    t = t + dt;
end

h = alt;    % get max altitude


end

