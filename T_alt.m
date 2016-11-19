function [ T ] = T_alt( h )
% determines atmospheric temperature given altitude above sea level

T0 = 288.15;    % Sea level std temp [K]
L = 0.0065;     % Temperature lapse rate [K/m]

T = T0 - L*h;   % Temperature at altitude [K]

end

