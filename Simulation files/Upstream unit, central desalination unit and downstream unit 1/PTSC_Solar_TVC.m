function [collector_CC, collector_OC,TES_cost, BHX_CC,BHX_OC, A_solar,q_ra0] = PTSC_Solar_TVC(Qprocess_1,Mms_1,Pms_1,Tv_1_1,Water_utility)

%--------------------------------------------------------- HTC Correlations ---------------------------------------------------------%
Ucon=@(T)(1.7194+(3.2063*10^(-2))*T+(-1.5971*10^(-5))*(T*T)-(1.9918*10^(-7))*(T*T*T)); %HTC of condenser @Tv kW/m2-deg C

Qprocess=Qprocess_1; %kJ/s
Mms=Mms_1; %kg/s
Pms=Pms_1; %kPa 
Tv_1=Tv_1_1; % oC

% T = 340; %in K
% X = [0.133 0.0067 0.6543 0.215];
%Thermal model of parabolic trough solar collector PTSC
% input parameters of the receiver
D2=0.047; % meter, inner diameter of the absorber
D3=0.051; %meter, outer diameter of the absorber
D4=0.070; %meter, inner diameter of the glass envelope
D5=0.074; %meter, outer diameter of the glass envelope
Dri = D2;
Dro = D3;
Dgi = D4;
Dgo = D5;
L_collector = 6.1; %meter, length of the collector
W = 2.3; % width of the collector
%% -----------
% Ar = D5*L_collector ; %unit: m2; reciver aperture area
Aa = (W-Dgo)*L_collector ;
Ar = pi*Dro*L_collector ;
Ac = pi*Dgo*L_collector ;
%%
In = 0.48; % intercept factor
TRS_env = 0.85; % transmittance of the glass envelope
ABS_abs = 0.95; % absorbance of the absorber
ABS_env = 0.02; % Absorptance of the glass envelope
refl_mirror = 0.9; % mirror reflectance
%% Input variables
% FL=1; % gpm, flow rate of the HTF
% FL = 7.742; % kg/sec
% ma = 0.03; % unit kg/sec.m2
%  mc = 32; % unit kg/sec
 mc=Water_utility; %kg/s
%  Tin = 312;% in K
 Tin = 31.5+273;% in K
%  Tout = 348; %in K
%
u5 = 5; % m/s, wind speed
T6 = 31.5+273; % K, ambient temperature
Tam = T6;
T7 = T6-8; %K, sky temperature
Month = 6; % month's number
n = 166; % day's number
% T1 = 45.3; % Inlet temperature
Sstand = 12.00; % local time (h)
%%
% Solar radiation code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sstand1 = Sstand*60;
B = (n-1)*360/365;
E = 229.2*(0.000075+0.001868*cosd(B)-0.032077*sind(B)-0.014615*cosd(2*B)-0.04089*sind(2*B));
% M = 4*(60-81)+E; % it is depends on longitude & altitude
L_st = 5.5 * 15; % Unit in o; standard meridian of local zone = 5.5(for khargapur) it is the time difference from GMT; 
L_loc = 87 ; % longitude of the location (kharagpur 87.23 W)
M = 4*(L_st-L_loc)+E; % it is depends on latitude & longitude 
ST = Sstand1+M;
HA = 15*(ST-12*60)/60; % Hour angle
SD = 23.45*sind(360*(284+n)/365); %Solar declination angle
SA = asind(sind(29.1735)*sind(SD)+cosd(29.1735)*cosd(SD)*cosd(HA));
%Solar altitude
ZA = 90-SA;%Zenith angle
AOI = ZA; % angle of incidence
%
A = [1202 1187 1164 1130 1106 1092 1093 1107 1136 1136 1190 1204];
Z = [0.141 0.142 0.149 0.164 0.177 0.185 0.186 0.182 0.165 0.152 0.144 0.141];
C = [0.103 0.104 0.109 0.120 0.130 0.137 0.138 0.134 0.121 0.111 0.106 0.103];

Ibn = A(Month)*exp(-Z(Month)/cosd(AOI));
Id = C(Month)*Ibn;      % diffuse radiation
Ib = Ibn*cosd(AOI);     % normal radiation
q_ra0 = (Ib+Id);  % global radiation (W/m2)
% q_ra = q_ra0*1.1; % W/m solar irradiation multiplied by the collector width
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optical model
f = 0.8; % focal distance
K = 3*10^-5*AOI^2-0.0072*AOI+1.2257; % Incidence angle modifier
Endloss = 1-(f/L_collector )*tand(AOI); % End loss
EEF_abs = In*TRS_env*ABS_abs*refl_mirror*K*Endloss; %effective optical efficiency

%% ==================  Properties by Group contribution   ====================== %%
%% -------------- density calculation ------------%
% Ro = DENS*X';
Ro = 1000 ; %water as HTF density in kg/m3
% Ro = 100.1 - 0.086*T - 0.0003*T^2;
%% speed of fluid calculation
% u1 = ma/Ro; % unit: m/sec
u1 = mc/(Ro*pi*D2^2);

%% ---- Specific heat ---- %
CP = 4.187; %4.187 ; % water as HTF in  kJ/kgK
% CP = (4.214 - (2.286*10^-3)*T + (4.991*10^-5)*T^2 - (4.519*10^-7)*T^3 + (1.857*10^-9)*T^4)*1000;
%% calculation of viscosity [Ref: Group Contribution Model for Predicting Viscosity of Fatty CompoundsRoberta Ceriani, Cintia B. Gonçalves, Juliana Rabelo, Marcel Caruso, Ana C. C. Cunha, Flavio W. Cavaleri, Eduardo A. C. Batista, and Antonio J. A. Meirelles; Journal of Chemical & Engineering Data 2007 52 (3), 965-972]
MU = 0.001002; %water as HTF in kg/(m·s)
% MU = (1.684*10^-3) - (4.264*10^-5)*T + (5.062*10^-7)*T^2 - (2.244*10^-9)*T^3;
%% Thermal conductivity %
% % Ke-Jun Wu, Chun-Xia Zhao, Chao-Hong He; Development of a group contribution method for determination of thermal conductivity of ionic liquids; Fluid Phase Equilibria; Volume 339,2013, Pages 10-14; 
% K_fi = [0.0964 0.1327 0.1393 0.1170]; % (Watt/m/K)calculated at steady state reboiler temperature (@400) in the Excel 
% K_FI = K_fi*X';
% K_FI = 0.598 ; % water as HTF in W/m·K
K_FI=0.598*10^(-3); % water as HTF in kW/mK
% K_FI = 0.5636 + (1.946*10^-3)*T - (8.151*10^-6)*T^2;
%%
% u1=0.000063*FL/(pi/4*(D2.^2)); % unit: m/s, RD outlet fluid speed
Re=(Ro*u1*D2)/MU; % Reynolds number of th HTH fluid
Pr=(MU*CP)/K_FI;
%% ----------
% %convection from HTH to the absorber
% Ro1=1001.1-0.0867*T1-0.0035*T1^2; % kg/m3, density of the fluid
% MU1=1.684*10^-3-4.264*10^-5*T1+5.062*10^-7*T1^2-2.244*10^-9*T1^3; %kg/m.s, viscosity of the fluid
% Cp1=(4.214-2.286*10^-3*T1+4.991*10^-5*T1^2-4.519*10^-7*T1^3+1.857*10^-9*T1^4)*1000; %J/kg.K, heat capacity of the fluid
% K1=0.5636+1.946*10^-3*T1-8.151*10^-6*T1^2; % W/m.K, conductivity of the fluid
% u1=0.000063*FL/(pi/4*(D2.^2)); %m/s, HTF fluid speed
% Re_D2=(Ro1*u1*D2)/MU1; % Reynolds number of th HTH fluid
% f2=(1.85*log10(Re_D2)-1.64).^-2; % friction factor
% Pr1=(MU1*Cp1)/K1;
%% --------------
Vair = 3; % m/sec; values for specific location(kharagpur)
h_ca = (8.6 * (Vair)^0.6)/L_collector ^0.4; %convection coefficient between cover and air
% K_fi = ; % conductivity of the fluid (here reboiler liquid)
Nu = 0.023*Re^0.8*Pr^0.4;
h_fi = (Nu*K_FI)/Dri ;%convective heat transfer co-efficient of fluid in the tube(absorber or,receiver)

%%
UL = (Ar/(h_ca * Ac))^-1; % loss-coefficient neglecting radiation heat loss terms
F_ce = (1/UL)/((1/UL)+(Dro/(h_fi*Dri)) + (Dro/2*K_FI)*(log(Dro/Dri))); % collector efficiency factor

% % ma=Ro*u1*(pi/4*(D2.^2)); % mass flow rate
%  mc = ma*Aa; % collector flow rate (kg/sec)

% FR = ((32*74941)/(Ar*UL))*(1-exp(-(Ar*UL*F_ce)/(32*74941))); % collector heat removal factor
FR = ((mc*CP)/(Ar*UL))*(1-exp(-(Ar*UL*F_ce)/(mc*CP))); % collector heat removal factor
%--------------
Qu = FR*(Ib * EEF_abs *Aa - Ar*UL*(Tin - Tam)); %  (J/sec) useful energy gain by the collector
T_eff = Qu/(q_ra0*Aa);
% dfdd = ((Ib+Id)*Aa)*0.29
% end
%% -------------------molten salt (HitecXL)
% Qms = mc*CP*(Tout-Tin);
% N_Trough = (Qms*0.0000167)/(Qu*0.001); % number of solar trough collector
% A_solar = Aa*N_Trough;  % total aperture area needed

%% -----------------------------------------------------
% Qprocess = 23776.6523596945*1000 ; % convert kJ/sec to J/sec (3rd work: saturated water is converted to saturated vapour. Sensible plus latent heat)
% Qprocess = 5177.82476658137*1000 ;
N_Trough = (Qprocess*1000)/(Qu); % number of solar trough collector
A_solar = Aa*N_Trough;  % total aperture area needed
% Tout = Tin + (Qprocess/(A_solar*UL*FR))*(1-FR); %receiver tube outelet temperature
% Tout = Tin + (Qprocess/(mc*CP*1000)); %receiver tube outelet temperature
Tout = Tin + (Qprocess/(mc*CP)); %receiver tube outelet temperature
% % Tout = Tenv + ((((Q*0.0167)/mc)-lamda_water)/CP);
%% --------------cost----------
A_collector = (W/2*(1+W^2/(16*f^2))^0.5 + 2*f*(W/(4*f) + (1+(W^2/(16*f^2))^0.5)))*L_collector  ; % area of the collector in m2
total_collector_area = 4*A_collector;
% collector_cost = 150*(A_collector*(N_Trough)); % solar collector capital cost (in $)
collector_CC = 150*((total_collector_area*(N_Trough)))^0.944; % solar collector capital cost (in $)
collector_OC = collector_CC *0.15;
collector_TCC = collector_CC + collector_OC;

% A_solar = (Aa/W)*L_collector*(N_Trough-1)+W;
% collector_cost = 150*A_solar;

% Total_collector_cost = collector_TCC*4;
%% ------------effective cost for whole day
% F_in_MED = 74.87; % kg/sec
% MF = F_in_MED/mc;  % multiplication factor for whole day and supply Qprocess
% total_collector_cost = collector_cost * MF; % solar collector capital cost (in $)
%% Thermal energy storage 
m_htf_LCA = Qprocess/(4.2*60); % (kg/sec)//mass flow of HTF in 1st effect // CP water = 42 kJ/kg-K //temp difference for sensible heat transfer = 15 K
m_htf = mc; % (kg/sec)//mass flow of HTF in 1st effect // CP water = 42 kJ/kg-K //temp difference for sensible heat transfer = 15 K
MST = 18; % max storage time in hours
mc_volume_total=(mc/1000)*3600*8760*25
mc_total=mc_volume_total*1000
V_TES=(MST*m_htf)/1000 % m3 // volume of TES// density of water =1000 kg/m3
V_TES_LCA = (MST*m_htf_LCA*3600)/1000 % m3 // volume of TES// density of water =1000 kg/m3
% TES_cost = 3458.6 * (V_TES)^0.6589 * 1.05; % euro to usd exchange = 1.05
TES_cost=2*V_TES_LCA*60;


%% BHX
A_BHX=Qprocess/(Ucon(Tv_1)*(Tsaturation(Pms)-Tv_1));
delT1=Tv_1-31.5;
delT2=Tsaturation(Pms)-36;
LMTD_BHX=(delT1-delT2)/log((delT1)/(delT2))
Area=Qprocess/(Ucon(Tv_1)*LMTD_BHX)
BHX_CC=150*(Area^0.8);
BHX_OC=0.25*BHX_CC;
%% cost of cooling water in the heat exchanger
F_cw = 6.7945; % kg/sec
cw_kg = F_cw*3600*8760; %(kg) total cooling water required per year 
cw_ton = cw_kg/1000; %(ton) total cooling water required per year 
cost_cw = cw_ton*0.03; % $ 
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
%% ------------------------------------------Tsaturation function-----------------------------------
function T=Tsaturation(P)
% P is in kPa
% t is in degree C
T=(42.6776-(3892.7/(log(P/1000)-9.48654)))-273.15;
end