function [ rho ] = rho_alt( h )
%gives density of air as a function of altitude

%% Constants
P0 = 101325;    % Sea level std press [Pa]
T0 = 288.15;    % Sea level std temp [K]
g = 9.80665;    % gravitational acc at surface [m/s^2]
L = 0.0065;     % Temperature lapse rate [K/m]
R = 8.31447;    % Universal gas constant [J/mol-K]
M = 0.0289644;  % Molar mass of dry air [kg/mol]

h0 = 1401;      % Altitude of Launch Site [m]


%% Density Variation Calculations
T = T0 - L*(h+h0);   % Temperature at altitude [K]
P = P0*(1 - L*(h+h0)/T0)^(g*M/R/L);  % Pressure at altitude [Pa]

rho = P*M/R/T;


end

