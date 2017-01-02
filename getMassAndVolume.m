function [Mf,Mox,Mpress,Vf,Vox,Vpress,W_f,W_ox,W_press,M_tank_f,...
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
W_press = (Vpress - pi*d_tank^3/6)/(pi*d_tank^2/4);

% determine masses of tanks, assuming a safety factor on yield stress
M_tank_f = 0.5*pi*d_tank^2/2*(d_tank/2+W_f)*FS_tank*P_f*rho_tank/sig_tank;
M_tank_ox = 0.5*pi*d_tank^2/2*(d_tank/2+W_ox)*FS_tank*P_ox*rho_tank/sig_tank;
M_tank_press = 0.5*pi*d_tank^2/2*(d_tank/2+W_press)*FS_tank*P_press*...
                rho_press/sig_tank;

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

