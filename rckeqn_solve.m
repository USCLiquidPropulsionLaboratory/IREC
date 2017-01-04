function [ F ] = rckeqn_solve(Mp,q1,q2,h)

% get altitude of apex with a specified propellant mass
[~,~,~,alt,~] = rckeqn(Mp,q1,q2);

% output difference between calculated altitude and desired altitude.
% Solver function will attempt to make this difference zero by adjusting
% propellant mass
F = alt - h;

end

