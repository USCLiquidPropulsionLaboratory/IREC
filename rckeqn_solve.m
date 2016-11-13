function [ F ] = rckeqn_solve(Mp,mdot,Ms_0,Ml,g,Isp,alpha,h,A,Cd)

% get altitude of apex with a specified propellant mass
[~,~,~,alt,~] = rckeqn(Mp,mdot,Ms_0,Ml,g,Isp,alpha,A,Cd);

% output difference between calculated altitude and desired altitude.
% Solver function will attempt to make this difference zero by adjusting
% propellant mass
F = alt - h;

end
