%   Rahul Rughani
%   USC SERC
%   LPL IREC - Altitude Calcs
%   November 12, 2016

clc;
format long;


%% Define Inputs

% Requirements
h =     30000;  %   altitude [ft]
I_max = 4e4;    %   max total impulse [N-sec]
Ml = 4;         %   payload mass [kg]

% Engine Properties
T =     4400;   %   thrust [N]
Isp =   295;    %   Isp [sec]
mdot =  1.5;    %   mass flow rate [kg/s]

% Rocket Properties
Ms_0 =  40;     %   fixed mass [kg]
alpha = 0.3;    %   fuel mass overhead
of_ratio = 1.8; %   Oxidizer to fuel mass ratio
rho_p = 810;    %   propellant density [kg/m^3]
P_p = 500;      %   Propellant pressure [Psi]
P_ox = 2000;    %   ox tank pressure [Psi]

% Drag Properties    
d = 0.15;       %   rocket diameter [m]
Cd = 0.5;       %   drag coefficient 

% Physical Properties
g = 9.8066;     %   gravitational accel at Earth's surface [m/s^2]
Req = 6378.14e3;%   Earth's equatorial radius [m]
Rp = 6356.8e3;  %   Earth's polar radius [m]  
lat = 34.42132; %   Latitude of SpacePort America [deg]
mu = 3.986004418e14;    %   Earth's gravitational parameter [m^2/s^3]  
R_univ = 8.3144598;     %   Universal gas constant [J/mol-K]
MM_ox = 31.9988e-3;     %   Molar mass of oxygen gas (O2) [kg/mol]
T_amb = 298.15; %   Ambient temperature [K]

rho_tank = 8000;%   Density of 304 SS [kg/m^3]
sig_tank = 2.15e8;      %   Yield Stress of 304 SS [Pa]
FS_tank = 1.2;  %   Safety factor on tank yielding


%% Calculations

P_p = P_p * 6894.75729;     % convert to Pascals
P_ox = P_ox * 6894.75729;   % convert to Pascals
rho_ox = P_ox/(R_univ/MM_ox)/T_amb; % get density of oxygen in tank

R_sp = sqrt( ((Req^2*cosd(lat))^2 + (Rp^2*sind(lat))^2) / ...
             ((Req*cosd(lat))^2 + (Rp*sind(lat))^2) );

h = h*0.3048;   % convert to meters

Mp = 10;    % initial guess of propellant mass [kg]

% create function handle for solver input (solver requires a function
% with only one input variable, so need to specify other inputs in a 
% function handle)
q = [mdot,Ms_0,Ml,g,Isp,alpha,d,Cd,of_ratio,R_sp,mu,rho_p,rho_ox,...
    rho_tank,sig_tank,FS_tank,P_p,P_ox];
f = @(x) rckeqn_solve(x,q,h);

% solve for propellant mass required to achieve desired altitude
Mp = fsolve(f,Mp);

I = Mp*g*Isp;   % compute total impulse [N-sec]

% get mass values of rocket from converged propellant mass
Ms = Ms_0 + alpha*Mp;
M0 = Ms + Mp + Ml;
Mb = Ms + Ml;

% get max acceleration
a_max = T/Mb;       % [m/s^2]
g_max = a_max/g;    % [g's]

% get rocket mass ratio
R = M0/Mb;

% compute burnout and apex properties with converged propellant mass
[hb,ub,tb,h,t] = rckeqn(Mp,q);

% write results into a table
Results = [Mp;M0;Mb;a_max;g_max;R;tb;t;hb;h;ub;I];
Rows = {'Propellant Mass';'Wet Mass';'Dry Mass';'Max Accel';'Max Gs';...
    'R';'burnout time';'apex time';'burnout alt';'apex alt';...
    'Burnout vel';'Total Impulse'};

Units = {'kg';'kg';'kg';'m/s^2';'g';'';'sec';'sec';'m';'m';'m/s';'N-sec'};

ResultsTable = table(Results,Units,'RowNames',Rows);

% print out table
ResultsTable

% write table to file
writetable(ResultsTable,'Results.txt','Delimiter','\t','WriteRowNames',true);



