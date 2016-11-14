function [ F ] = rckeqn_solve(Mp,q,h)

% get altitude of apex with a specified propellant mass
[~,~,~,alt,~] = rckeqn(Mp,q);

% output difference between calculated altitude and desired altitude.
% Solver function will attempt to make this difference zero by adjusting
% propellant mass
F = alt - h;

end

