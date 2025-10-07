%Sensitivity analysis of capital and energy production cost for off-grid building integrated photovoltaic systems
% DOI: https://doi.org/10.1016/j.renene.2022.01.003
function [c,Area_PV]=Cost_PV(Total_electricity)
fg=1.5;
S=6; %hours
El=Total_electricity*1000*24; %Wh/day It is assumed that the total electricity is provided by solar PV for the entire day
ED=El*fg;
Wp=ED/S;        % S is solar hours, assuming 6 hours of solar energy is available
Cbar_array=4;   % $/W
Cbar_bat=1.25;  % $/Ahr
Vbat=12;        % Single battery nominal voltage 
D=3;            % Autonomy days
Cbar_cc=5;      % $/A
Vo=110*1000;    % Volts
Cbar_inv=2;     % $/W
finv=1.1;       

Carray=Wp*Cbar_array;
Cbat=((ED*D)/Vbat)*Cbar_bat;
Ccc=(Wp/Vo)*Cbar_cc;
Cinv=finv*Wp*Cbar_inv;


Ccap=Carray+Cbat+Ccc+Cinv;
Ir=0.07;
life=25;
Capital_recovery_factor=(Ir*((1+Ir)^life))/(((1+Ir)^life)-1);
fmaint_PV=0.1;
TAC_PV=(1+fmaint_PV)*Capital_recovery_factor*(Ccap);
Coperating=0.05*Ccap;       % Techno-economic analysis of solar photovoltaic (PV) and solar photovoltaic thermal (PVT) systems using exergy analysis (https://doi.org/10.1016/j.seta.2021.101520)
c=TAC_PV;
% Capacity of one solar PV of 156*156 cm2 has a capacity of 210 Wp
% (Ecoinvent)
Area_PV=(((156*156)/210)*(ED))*(10^(-4));
% Carray=Capital cost of PV array
% Cbat=Capital cost of batteries
% Cinv=Capital cost of inverter
% Cc=Capital cost of charge controller
end