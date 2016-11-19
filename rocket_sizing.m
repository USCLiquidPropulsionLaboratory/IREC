%   Rahul Rughani
%   USC Liquid Propulsion Lab
%   LPL IREC - Altitude Calcs
%   November 12, 2016

%   TO DO LIST:
%   - Shift atmospheric density profiles to local conditions, as we are not
%     starting at sea level
%   - Compute weight of rocket body as function of tank diameters and
%     heights rather than overhead percentage on propellant mass
%   - Modify fixed mass to be more accurate than an arbitrary number.
%     Perhaps set it to a sum of smaller fixed masses, like bulkhead masses
%     and nosecone, etc.
%   - Highlight max values on plots (probably steal code snippets from
%     Becca or something)
%   - Size piping mass accurately
%   - Get masses of valves currently used on Blue Steel to size valves and
%     ports
%   - Figure out the actual drag coefficient of a similarly sized rocket,
%     or implement some sort of code to solve for it
%   - Use a variable Cd coefficient that varies with Mach number (will
%     definietly make results more accurate)
%   - Can we use Carbon Fiber fuel tanks? This would definitely make the
%     vehicle lighter and use less fuel, making it even lighter, etc.


clc;
format long;

%% User Specific - You may want to change this to match your preferences

dir = '~/Desktop/IREC/'; % location of output files
% NOTE: WINDOWS USERS NEED TO USE \ instead of / IN DIRECTORY
% SPECIFICATION. / IS RESERVED FOR UNIX SYSTEMS


%% Define Inputs

% Requirements
h =     30000;  %   altitude [ft]
I_max = 4.04e4; %   max total impulse [N-sec]
Ml = 4;         %   payload mass [kg]

% Engine Properties
T =     4400;   %   thrust [N]
Isp =   295;    %   Isp [sec]
mdot =  1.5;    %   mass flow rate [kg/s]

% Rocket Properties
Ms_0 =  30;     %   fixed mass [kg]
M_bulkhead = 2; %   bulkhead mass [kg]
alpha = 0.2;    %   fuel mass overhead [-]
of_ratio = 1.8; %   Oxidizer to fuel mass ratio
rho_f = 810;    %   propellant density [kg/m^3]
P_f = 500;      %   Propellant pressure [Psi]
P_ox = 2000;    %   ox tank pressure [Psi]
d = 0.20;       %   rocket diameter [m]
d_tank = d - 0.02;  %   fuel and ox tank diameters [m]

% Drag Properties    
Cd = 0.2;       %   drag coefficient 

% Physical Properties
g = 9.8066;     %   gravitational accel at Earth's surface [m/s^2]
Req = 6378.14e3;%   Earth's equatorial radius [m]
Rp = 6356.8e3;  %   Earth's polar radius [m]  
lat = 34.42132; %   Latitude of SpacePort America [deg]
mu = 3.986004418e14;    %   Earth's gravitational parameter [m^2/s^3]  
R_univ = 8.3144598;     %   Universal gas constant [J/mol-K]
MM_ox = 31.9988e-3;     %   Molar mass of oxygen gas (O2) [kg/mol]
MM_air = 28.964e-3;     %   Molar mass of air [kg/mol]
T_amb = 298.15; %   Ambient temperature (Sea Level) [K]
gam = 1.4;      %   Specific heat capacity ratio for air [-]

% fuel and ox tank material properties [kg/m^3] & [Pa]
rho_tank = 8000;        %   Density of 304 SS [kg/m^3]
sig_tank = 500e6;       %   Yield Stress of 304 SS [Pa]
% rho_tank = 1800;        % Carbon Fiber
% sig_tank = 6370e6;      % Carbon Fiber
FS_tank = 1.2;  %   Safety factor on tank failure

% Simulation Properties
dt = 0.01;      %   simulation time-step [sec]


%% Iterate to Solve for Propellant Mass

P_f = P_f * 6894.75729;     % convert to Pascals
P_ox = P_ox * 6894.75729;   % convert to Pascals
rho_ox = P_ox/(R_univ/MM_ox)/T_amb; % get density of oxygen in tank

R_sp = sqrt( ((Req^2*cosd(lat))^2 + (Rp^2*sind(lat))^2) / ...
             ((Req*cosd(lat))^2 + (Rp*sind(lat))^2) );

h = h*0.3048;   % convert to meters

Mp = 10;    % initial guess of propellant mass [kg]

% create function handle for solver input (solver requires a function
% with only one input variable, so need to specify other inputs in a 
% function handle)
q = [mdot,Ms_0,Ml,g,Isp,alpha,d,Cd,of_ratio,R_sp,mu,rho_f,rho_ox,...
    rho_tank,sig_tank,FS_tank,P_f,P_ox,d_tank,dt];
f = @(x) rckeqn_solve(x,q,h);

% solve for propellant mass required to achieve desired altitude
Mp = fzero(f,Mp);


%% Compute Properties for Converged Solution

I = Mp*g*Isp;   % compute total impulse [N-sec]

if I > I_max
    error(['Solution Error: Total impulse is larger than max permissible'...
           'value, aborting solution']);
end

Mf = Mp/(1+of_ratio);   % fuel mass
Mox = Mp - Mf;          % oxidizer mass

Vf = Mp/rho_f;      % volume of propellant
Vox = Mox/rho_ox;   % volume of oxidizer

% we are using cylindrical tanks with domed ends, so lets find the height
% of the cylindrical portion of the tanks
W_f = (Vf - pi*d_tank^3/6)/(pi*d_tank^2/4);  
W_ox = (Vox - pi*d_tank^3/6)/(pi*d_tank^2/4);

% determine total heigh of tanks
h_tank_f = W_f + d_tank;
h_tank_ox = W_ox + d_tank;

% determine masses of tanks, assuming a safety factor on yield stress
M_tank_f = pi*d_tank^2/2*(d_tank/2+W_f)*FS_tank*P_f*rho_tank/sig_tank;
M_tank_ox = pi*d_tank^2/2*(d_tank/2+W_ox)*FS_tank*P_ox*rho_tank/sig_tank;

% get mass values of rocket from converged propellant mass
Ms = Ms_0 + M_tank_f + M_tank_ox + alpha*Mp;
M0 = Ms + Mp + Ml;
Mb = Ms + Ml;

% get rocket mass ratio
R = M0/Mb;

% compute burnout and apex properties with converged propellant mass
[hb,ub,tb,h,t] = rckeqn(Mp,q);

[alt,vel,temp,time] = rckeqn_hist(Mp,q);%   get altitude and velocity, and 
                                        %   air temp vs time

tb_index = find(time==tb);  %   find index of burnout time

% compute max acceleration, using numerical differentiation over timestep
a_max = (vel(tb_index) - vel(tb_index-1))/dt; % [m/s^2]
g_max = a_max/g;    % [G's]


%% Tabulate Results

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
writetable(ResultsTable,[dir,'Results.txt'],'Delimiter','\t','WriteRowNames',true);

%% Plot Data for Visualization

%   Generate plots of altitude vs time, velocity vs time, Mach number vs
%   time, and Mach number vs altitude

% obtain local speed of sound at each time-step
a = (gam*(R_univ/MM_air)*temp).^0.5;    %   [m/s]

% obtain Mach number at each time-step
Mach = vel./a;                          %   [-]

figure(1);
subplot(2,2,1);                     %   begin a subplot (2x2, first cell)
plot(time,alt);                     %   plot primary curve
hold on;                               
plot([tb,tb],[0,max(alt)],'--');    %   plot vertical line marking burnout
xlabel('time (sec)');               %   label X-axis
ylabel('altitude (m)');             %   label Y-axis
title('Altitude vs Time');          %   make title
grid on;                            %   show gridlines on plot

subplot(2,2,2);
plot(time,vel);
hold on;
plot([tb,tb],[0,max(vel)],'--');
xlabel('time (sec)');
ylabel('velocity (m/s)');
title('Velocity vs Time');
grid on;

subplot(2,2,3);
plot(time,Mach);
hold on;
plot([tb,tb],[0,max(Mach)],'--');
xlabel('time (sec)');
ylabel('Mach number');
title('Mach number vs Time');
grid on;

subplot(2,2,4);
plot(alt,Mach);
hold on;
plot([hb,hb],[0,max(Mach)],'--');
xlabel('altitude (m)');
ylabel('Mach number');
title('Mach number vs Altitude');
grid on;

saveas(gcf,[dir,'Ascent_Profile.png']); %   save figure as PNG file












