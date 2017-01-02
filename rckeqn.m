function [ hb,ub,tb,h,t  ] = rckeqn( Mp,q )
% solve the rocket equation with drag to obtain the burnout properties as
% well as the max altitude, assuming perfectly vertical flight

% simulation is currently 1-D, with positive being up

%% Extract Data

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
rho_press = q(14);
rho_tank = q(15);
sig_tank = q(16);
FS_tank = q(17);
P_f = q(18);
P_ox = q(19);
P_press = q(20);
d_tank = q(21);
dt = q(22);

%% Compute Masses and Volumes

% calculate cross-sectional area of rocket
A = 0.25*pi*d^2;

% split propellant mass into fuel mass and oxidizer mass
Mf = Mp/(1+of_ratio);   % fuel mass
Mox = Mp - Mf;          % oxidizer mass

% get volumes
Vf = Mf/rho_f;      % volume of propellant
Vox = Mox/rho_ox;   % volume of oxidizer

% lets calculate the pressurant volume based off the fuel properties
R_univ = 8.3144598;
MM_He = 4.002602e-3;
R_He = R_univ/MM_He;
T_uf = 298;
P_u = 7e5; 
P_comb = 7e6;
gam_He = 1.33;
Zuf = 1;    % assume ideal gas for now
Mpress = P_u*Vf/R_He/T_uf/Zuf/(1 - (P_comb/P_press)^(1/gam_He));
Vpress = Mpress/rho_press;

% we are using cylindrical tanks with domed ends, so lets find the height
% of the cylindrical portion of the tanks
W_f = (Vf - pi*d_tank^3/6)/(pi*d_tank^2/4);  
W_ox = (Vox - pi*d_tank^3/6)/(pi*d_tank^2/4);
W_press = (Vpress - pi*d_tank^3/6)/(pi*d_tank^2/4);

% determine masses of tanks, assuming a safety factor on yield stress
M_tank_f = 0.5*pi*d_tank^2/2*(d_tank/2+W_f)*FS_tank*P_f*rho_tank/sig_tank;
M_tank_ox = 0.5*pi*d_tank^2/2*(d_tank/2+W_ox)*FS_tank*P_ox*rho_tank/sig_tank;
M_tank_press = 0.5*pi*d_tank^2/2*(d_tank/2+W_press)*FS_tank*P_press*...
                rho_press/sig_tank;

% determine mass of airframe and housing (estimated as an overhead to the
% propellant mass
M_frame = alpha*Mp;

% Assume structural mass of rocket is the base structural mass, plus
% tankage masses, and an overhead for the rest of the fuselage
Ms = Ms_0 + M_tank_f + M_tank_ox + M_tank_press + M_frame;

% calculate wet and dry masses of rocket
M0 = Ms + Mp + Mpress + Ml;
Mb = Ms + Ml;

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

