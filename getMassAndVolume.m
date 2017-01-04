function [Mf,Mox,Mpress,Vf,Vox,Vpress,W_f,W_ox,d_tank_press,M_tank_f,...
    M_tank_ox,M_tank_press,M_frame,Ms,M0,Mb] = getMassAndVolume( q2,Mp )
% compute rocket masses and volumes

% extract data from array
Ms_0 = q2(1);
Ml = q2(2);
alpha = q2(3);
of_ratio = q2(4);
rho_f = q2(5);
rho_ox = q2(6);
rho_press = q2(7);
rho_tank = q2(8);
sig_tank = q2(9);
FS_tank = q2(10);
P_f = q2(11);
P_ox = q2(12);
P_press = q2(13);
d_tank = q2(14);
R_univ = q2(15);
MM_He = q2(16);
gam_He = q2(17);
Zuf = q2(18);
P_u = q2(19);
P_comb = q2(20);
T_uf = q2(21);



% split propellant mass into fuel mass and oxidizer mass
Mf = Mp/(1+of_ratio);   % fuel mass
Mox = Mp - Mf;          % oxidizer mass

% get volumes
Vf = Mf/rho_f;      % volume of propellant
Vox = Mox/rho_ox;   % volume of oxidizer

% lets calculate the pressurant volume based off the fuel properties
R_He = R_univ/MM_He;
Mpress = P_u*Vf/R_He/T_uf/Zuf/(1 - (P_comb/P_press)^(1/gam_He));
Vpress = Mpress/rho_press;

% we are using cylindrical tanks with domed ends, so lets find the height
% of the cylindrical portion of the tanks
W_f = (Vf - pi*d_tank^3/6)/(pi*d_tank^2/4);  
W_ox = (Vox - pi*d_tank^3/6)/(pi*d_tank^2/4);


% get tank thicknesses
t_tank_f = FS_tank*P_f*d_tank/2/sig_tank;  
t_tank_ox = FS_tank*P_ox*d_tank/2/sig_tank;

% get volume reduction of usable space due to non-zero tank thickness
dV_f = 0.25*pi*(d_tank^2 - (d_tank-t_tank_f)^2)*W_f + pi*(d_tank^3 - ...
    (d_tank - t_tank_f)^3)/6;   
dV_ox = 0.25*pi*(d_tank^2 - (d_tank-t_tank_ox)^2)*W_f + pi*(d_tank^3 - ...
    (d_tank - t_tank_ox)^3)/6;   

% get new tank heights required to store specified fluid volume within the
% specified tank diameter (tanks now have thickness)
W_f = W_f + dV_f/(0.25*pi*(d_tank-t_tank_f)^2);
W_ox = W_ox + dV_ox/(0.25*pi*(d_tank-t_tank_ox)^2);

% determine masses of tanks, assuming a safety factor on yield stress
M_tank_f = 0.5*pi*d_tank^2/2*(d_tank/2+W_f)*FS_tank*P_f*rho_tank/sig_tank;
M_tank_ox = 0.5*pi*d_tank^2/2*(d_tank/2+W_ox)*FS_tank*P_ox*rho_tank/sig_tank;

            
% for the pressurant tank, we don't need much volume so we'll use a sphere
% instead of domed cylinders
d_tank_press = (3*Vpress/4/pi)^(1/3);   % inner diameter of tank
t_tank_press = FS_tank*P_press*d_tank_press/4/sig_tank; % press tank thickness
M_tank_press = 3/2*FS_tank*P_press*Vpress*rho_tank/sig_tank;
d_tank_press = d_tank_press + 2*t_tank_press;   % convert tank inner diameter 
                                                % to be outer diameter

% determine mass of airframe and housing (estimated as an overhead to the
% propellant mass
M_frame = alpha*Mp;

% Assume structural mass of rocket is the base structural mass, plus
% tankage masses, and an overhead for the rest of the fuselage
Ms = Ms_0 + M_tank_f + M_tank_ox + M_tank_press + M_frame;

% calculate wet and dry masses of rocket
M0 = Ms + Mp + Mpress + Ml;
Mb = Ms + Ml;

end

