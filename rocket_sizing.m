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
%   - Create some simulated curves to see what happens if we miss our burn
%     target (i.e. go over or under by 1/10th of a second)
%   - Get actual radius of Earth at Spaceport America using Earth geodesic
%     model rather than oblate spheroid model
%   - Vary thrust with ambient pressure (rather than using Isp, use thrust
%     equation from nozzle


clc;
format long;

%% User Specific - You may want to change this to match your preferences

dir = '~/Documents/Data Files/USC/LPL/IREC/'; % location of output files
% NOTE: WINDOWS USERS NEED TO USE \ instead of / IN DIRECTORY
% SPECIFICATION. / IS RESERVED FOR UNIX SYSTEMS


%% Define Inputs

% Requirements
h =     30000;  %   altitude [ft]
I_max = 40960;  %   max total impulse [N-sec]
Ml = 4;         %   payload mass [kg]

% Engine Properties
T =     4400;   %   thrust [N]
Isp =   250;    %   Isp [sec]
mdot =  1.5;    %   mass flow rate [kg/s]

% Rocket Properties
Ms_0 =  40;     %   fixed mass [kg]
M_bulkhead = 2; %   bulkhead mass [kg]
alpha = 0.2;    %   fuel mass overhead [-]
of_ratio = 1.8; %   Oxidizer to fuel mass ratio
rho_f = 810;    %   propellant density [kg/m^3]
P_f = 15;       %   Propellant pressure [Psi]
P_ox = 3500;    %   ox tank pressure [Psi]
P_He = 2000;    %   pressurant tank pressure [Psi]
P_u = 100;      %   ullage pressure [Psi]
P_comb = 1200;  %   combustion pressure [Psi]
d = 8;          %   rocket diameter [in]
d_tank = 6.5;   %   fuel, ox and pressurant tank diameters [in]

% Drag Properties    
Cd = 0.2;       %   drag coefficient 
visc = 1.81e-5; %   kinematic viscosity of air [Pa-s]

% Physical Properties
g = 9.8066;     %   gravitational accel at Earth's surface [m/s^2]
Req = 6378.14e3;%   Earth's equatorial radius [m]
Rp = 6356.8e3;  %   Earth's polar radius [m]  
lat = 34.42132; %   Latitude of SpacePort America [deg]
mu = 3.986004418e14;    %   Earth's gravitational parameter [m^2/s^3]  
R_univ = 8.3144598;     %   Universal gas constant [J/mol-K]
MM_ox = 31.9988e-3;     %   Molar mass of oxygen gas (O2) [kg/mol]
MM_He = 4.002602e-3;    %   Molar mass of helium gas (He) [kg/mol]
MM_air = 28.964e-3;     %   Molar mass of air [kg/mol]
T_amb = 298.15; %   Ambient temperature (Sea Level) [K]
gam = 1.4;      %   Specific heat capacity ratio for air [-]
gam_He = 1.66;  %   Specific heat capacity ratio for Helium [-]
Zuf = 1;        %   compressibility of gas (assumed ideal)

% fuel and ox tank material properties [kg/m^3] & [Pa]
rho_tank = 8000;        %   Density of 304 SS [kg/m^3]
sig_tank = 500e6;       %   Yield Stress of 304 SS [Pa]
% rho_tank = 1800;        % Carbon Fiber
% sig_tank = 6370e6;      % Carbon Fiber
FS_tank = 1.2;  %   Safety factor on tank failure

% Simulation Properties
dt = 1e-4;      %   simulation time-step [sec]

% Conversion Factors
psi2pa = 6894.75729;    %   PSI to Pascals
ft2m = 0.3048;          %   feet to meters
ft2in = 12;             %   feet to inches
in2m = ft2m/ft2in;      %   inches to meters
lbf2N = 4.44822;        %   pound-force to Newtons
lbm2kg = 0.453592;      %   pound-mass to kilograms

%% Iterate to Solve for Propellant Mass

d = d * in2m;           % convert to meters
d_tank = d_tank * in2m; % convert to meters

P_f = P_f * psi2pa;     % convert to Pascals
P_ox = P_ox * psi2pa;   % convert to Pascals
P_He = P_He * psi2pa;   % convert to Pascals
P_u = P_u * psi2pa;     % convert to Pascals
P_comb = P_comb*psi2pa; % convert to Pascals
rho_ox = P_ox/(R_univ/MM_ox)/T_amb; % get density of oxygen in tank
rho_He = P_He/(R_univ/MM_He)/T_amb; % get density of helium in tank

% Radius of Earth at Spaceport America (approximated using oblate spheroid
% model, not geodesic patch model)
R_sp = sqrt( ((Req^2*cosd(lat))^2 + (Rp^2*sind(lat))^2) / ...
             ((Req*cosd(lat))^2 + (Rp*sind(lat))^2) );

h = h * ft2m;   % convert to meters

Mp = 10;    % initial guess of propellant mass [kg]

% create function handle for solver input (solver requires a function
% with only one input variable, so need to specify other inputs in a 
% function handle)
% q = [mdot,Ms_0,Ml,g,Isp,alpha,d,Cd,of_ratio,R_sp,mu,rho_f,rho_ox,rho_He,...
%     rho_tank,sig_tank,FS_tank,P_f,P_ox,P_He,d_tank,dt];
q1 = [mdot,g,Isp,d,Cd,R_sp,mu,dt,visc];
q2 = [Ms_0,Ml,alpha,of_ratio,rho_f,rho_ox,rho_He,rho_tank,sig_tank,...
        FS_tank,P_f,P_ox,P_He,d_tank,R_univ,MM_He,gam_He,Zuf,P_u,...
        P_comb,T_amb];
f = @(x) rckeqn_solve(x,q1,q2,h);

% solve for propellant mass required to achieve desired altitude
Mp = fzero(f,Mp);


%% Compute Properties for Converged Solution

I = Mp*g*Isp;   % compute total impulse [N-sec]

if I > I_max
    error(['Solution Error: Total impulse is larger than max permissible'...
           'value, aborting solution']);
end

% get rocket properties from converged propellant mass
[Mf,Mox,Mpress,Vf,Vox,Vpress,W_f,W_ox,d_tank_press,M_tank_f,M_tank_ox,...
    M_tank_press,M_frame,Ms,M0,Mb] = getMassAndVolume(q2,Mp);

% determine total height of tanks
h_tank_f = W_f + d_tank;
h_tank_ox = W_ox + d_tank;
h_tank_press = d_tank_press;

height = h_tank_f + h_tank_ox + h_tank_press;

% get rocket mass ratio
R = M0/Mb;

% compute burnout and apex properties with converged propellant mass
[hb,ub,tb,h,t] = rckeqn(Mp,q1,q2);

%   get altitude and velocity, and air temp vs time
[alt,vel,temp,time] = rckeqn_hist(Mp,q1,q2);
                                        

tb_index = find(time==tb);  %   find index of burnout time

% compute max acceleration, using numerical differentiation over timestep
a_max = (vel(tb_index) - vel(tb_index-1))/dt; % [m/s^2]
g_max = a_max/g;    % [G's]

% obtain local speed of sound at each time-step
a = (gam*(R_univ/MM_air)*temp).^0.5;    %   [m/s]

% obtain Mach number at each time-step
Mach = vel./a;                          %   [-]

% obtain Mach number at burnout
Mach_b = Mach(tb_index);                %   [-]


%% Tabulate Results

% write results into a table
Results_SI = [Mp;M0;Mb;a_max;g_max;R;tb;t;hb;h;ub;Mach_b;I;height];
Rows = {'Propellant Mass';'Wet Mass';'Dry Mass';'Max Accel';'Max Gs';...
    'R';'burnout time';'apex time';'burnout alt';'apex alt';...
    'Burnout vel';'Burnout Mach #';'Total Impulse';'Tankage Height'};

Units_Si = {'kg';'kg';'kg';'m/s^2';'g';'';'sec';'sec';'m';'m';'m/s';'';...
            'N-sec';'m'};
Units_ENG = {'lb';'lb';'lb';'ft/s^2';'g';'';'sec';'sec';'ft';'ft';'ft/s';...
             '';'lbf-sec';'ft'};

% conversion factors between English and SI units for each result type
ENG2SI = [lbm2kg;lbm2kg;lbm2kg;ft2m;1;1;1;1;ft2m;ft2m;ft2m;1;lbf2N;ft2m];

% convert results from SI to English units for secondary output
Results_ENG = Results_SI./ENG2SI;

ResultsTable = table(Results_SI,Units_Si,Results_ENG,Units_ENG,'RowNames',Rows);

% print out table
ResultsTable

% write table to file
writetable(ResultsTable,[dir,'Results.txt'],'Delimiter','\t','WriteRowNames',true);

%% Plot Data for Visualization

%   Generate plots of altitude vs time, velocity vs time, Mach number vs
%   time, and Mach number vs altitude


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












