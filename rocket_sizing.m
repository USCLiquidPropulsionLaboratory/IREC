%   Rahul Rughani
%   USC SERC
%   LPL IREC - Altitude Calcs
%   November 12, 2016

clc;
format long;


%% Define Inputs

h =     30000;  %   altitude [ft]
T =     4400;   %   thrust [N]
Isp =   295;    %   Isp [sec]
mdot =  1.5;    %   mass flow rate [kg/s]
Ms_0 =  50;     %   fixed mass [kg]
Ml = 4;         %   payload mass [kg]
alpha = 0.5;    %   fuel mass overhead
g = 9.8066;     %   gravitational accel at Earth's surface [m/s^2]
I_max = 4e4;    %   max total impulse [N-sec]

%rho = 1;    
d = 0.15;       %   rocket diameter [m]
Cd = 0.5;       %   drag coefficient 

h = h*0.3048;   % convert to meters
A = 0.25*pi*d^2; % cross sectional area [m^2]

%% Calculations

Mp = 10;    % initial guess of propellant mass [kg]

f = @(x) rckeqn_solve(x,mdot,Ms_0,Ml,g,Isp,alpha,h,A,Cd);

Mp = fzero(f,Mp);

I = Mp*g*Isp;

Ms = Ms_0 + alpha*Mp;
M0 = Ms + Mp + Ml;
Mb = Ms + Ml;

a_max = T/Mb;
g_max = a_max/g;

R = M0/Mb;

[hb,ub,tb,h,t] = rckeqn(Mp,mdot,Ms_0,Ml,g,Isp,alpha,A,Cd);

Results = [Mp;M0;Mb;a_max;g_max;R;tb;t;hb;h;ub;I];
Rows = {'Propellant Mass';'Wet Mass';'Dry Mass';'Max Accel';'Max Gs';...
    'R';'burnout time';'apex time';'burnout alt';'apex alt';...
    'Burnout vel';'Total Impulse'};

Units = {'kg';'kg';'kg';'m/s^2';'g';'';'sec';'sec';'m';'m';'m/s';'N-sec'};

ResultsTable = table(Results,Units,'RowNames',Rows);

ResultsTable



