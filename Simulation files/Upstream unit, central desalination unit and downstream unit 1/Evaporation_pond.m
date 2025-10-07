function[Apond,Capital_cost,Operating_cost]=Evaporation_pond(Xbn,Mcw_total,plant_life,Tcw,Tf,X)
%--------------------------------------Latent heat of vaporisation Correlation-------------------------------------------------------------%
Lamda=@(T) (2501.897149-2.407064037*T+(1.192217*10^-3)*T^2-(1.5863*10^-5)*T^3); % Latent heat of vaporisation kJ/kg
rhol=@(T,x) (1*10^3*(((4.032219*0.5+0.115313*((0.002*x-150)/150)+(3.26*10^-4)*(2*(((0.002*x-150)/150))^2-1))*(0.5))+...
            (-0.108199*0.5+(1.571*10^-3)*((0.002*x-150)/150)-(4.23*10^-4)*(2*(((0.002*x-150)/150))^2-1))*((2*T-200)/160)+...
            (-0.012247*0.5+(1.74*10^-3)*((0.002*x-150)/150)-(9*10^-6)*(2*(((0.002*x-150)/150))^2-1))*(2*(((2*T-200)/160)^2)-1)+...
            ((6.92*10^-4)*0.5-(8.7*10^-5)*((0.002*x-150)/150)-(5.3*10^-5)*(2*(((0.002*x-150)/150))^2-1))*...
            (4*(((2*T-200)/160)^3)-3*(2*T-200)/160))); %Density of liquid kg/m3
% x is in Ppm
%% --------------------------------------------------------------------------------------------------------------
% Methods for calculating brine evaporation rates during salt production (https://doi.org/10.1016/j.jas.2007.10.013)
% Each scenario involves evaporation of a fixed volume of brine to dryness resulting in the 
% precipitation of salt crystals. For any batch evaporation the amount of
% salt produced can be determined by:
T=(Tcw+Tf)/2;
Mcw_kg=Mcw_total*8760*3600*plant_life;
S=10^2*Xbn;
Msalt_kg=Mcw_kg*(1.52*10^(-4)*S^2+9.50*10^(-3)*S);
Latent_heat=Lamda(T)*10^(-3); % MJ/kg % Assuming ambient temperature to be 25 degree celcius
% Average solar radiation is 4.829 (kwh/m2*day)*6 hours/day*365 day; average solar radiation of west bengal 
% (Modelling an off-grid integrated renewable energy system for rural electrification in India using photovoltaics and anaerobic digestion) 
Solar_radiation=4.829*6;  % kW/m2day % Assuming 6 hours of solar radiation (https://doi.org/10.1016/j.renene.2014.08.055)
Rn=(Solar_radiation*8760*3600*25)*10^(-3); %MJ/m2day
u=2.6; % Average wind speed in Kolkata 
fu=6.43*(1+0.536*u);
es=Psaturation(Tf); %kPa
e=Psaturation(Tcw);%kPa
delta=(4098*es)/(237.3+T)^2;
gamma=0.000655*101.325; % gamma=0.000655*Patm
A=(delta/(delta+gamma))*Rn;
B=(gamma/(delta+gamma))*fu;

Qevaporation=(((Latent_heat)^(-1))*(A+B*(es-e)))*10^(-3); % m/day
brine_m3=Mcw_kg/rhol(Tf,X*10^6); % Total volume of brine
Qbrine=((brine_m3)/(8760*25))*24; % m3/day
Apond=Qbrine/Qevaporation; % m2

% Construction cost of Pond for handling 6000 m3/d of brine is 15600000 $. scaling up this value depending on the value of the brine flowrate.
% Reverse Osmosis Membrane Zero Liquid Discharge for Agriculture Drainage Water Desalination: Technical, Economic, and Environmental Assessment
%(https://doi.org/10.3390/membranes12100923)
Apond_acres=Apond*0.000247105;

% A feasibility study of a small-scale photovoltaic-powered reverse osmosis desalination plant for 
% potable water and salt production in Madura Island: A techno-economic evaluation
% (https://doi.org/10.1016/j.tsep.2022.101450)
Capital_cost=9632*Apond_acres;
Operating_cost=0.05*Capital_cost;

end
function y=Psaturation(T) % El-Dessouky correlation Appendix A1
Tc=647.286;
Pc=22089;
a=0;
f=[-7.419242 0.29721 -0.1155286 0.008685635 0.001094098 -0.00439993 0.002520658 -0.000521868];
for i=1:numel(f)
a=a+f(i)*(0.01*(T+273.15-338.15))^(i-1);
end
y=Pc*(exp((Tc/(T+273.15)-1)*a));
end