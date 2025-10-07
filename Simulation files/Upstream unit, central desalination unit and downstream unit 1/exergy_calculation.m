function [ef1]=exergy_calculation(P1,T1,S1,Tcw,Xf)

uT='C';
uS='w';
uP='kPa';

%Table 1
% T0=25;      %deg C
% P0=101.325; %kPa
% S0=0.035;  %kg/kg
% 
% P=101.325;
% T=90;
% S=0.035;

% Table 2
In=[P1 T1 S1];
T0=Tcw;      %deg C
P0=101.325; %kPa
S0=Xf;  %kg/kg

P=In(1);
T=In(2);
S=In(3);

% SW_ChemPot_s(T,uT,S,uS,P,uP)
% SW_ChemPot_w(T,uT,S,uS,P, uP)
% SW_Enthalpy(T,uT,S,uS,P,uP)
% SW_Entropy(T,uT,S,uS,P,uP)

h=SW_Enthalpy(T,uT,S,uS,P,uP);
h_star=SW_Enthalpy(T0,uT,S,uS,P0,uP);

s=SW_Entropy(T,uT,S,uS,P,uP);
s_star=SW_Entropy(T0,uT,S,uS,P0,uP);
% 
% uw_star=SW_ChemPot_w(T0,uT,S,uS,P0, uP);
% uw_not=SW_ChemPot_w(T0,uT,S0,uS,P0, uP);
% 
% us_star=SW_ChemPot_s(T0,uT,S,uS,P0,uP);
% us_not=SW_ChemPot_s(T0,uT,S0,uS,P0,uP);

%% Chemical Potential for zero S value
% Specific volume
a=[8.020*10^(2) -2.001  1.677*10^(-2) -3.060*10^(-5)    -1.613*10^(-5)];
rho_w=@(T,ws) (9.999*10^(2)+2.034*10^(-2)*T-6.162*10^(-3)*(T^2)+2.261*10^(-5)*(T^3)-...
               4.657*10^(-8)*T^(4));
rho_sw=@(T,ws) (rho_w(T,ws)+ws*(a(1)+a(2)*T+a(3)*(T^2)+a(4)*(T^3)+a(5)*ws*(T^2)));
v_sw=@(T,ws) (1/rho_sw(T,ws));


% Specific enthalpy
b1=-2.348*10^(4);
b2=3.152*10^(5);
b3=2.803*10^(6);
b4=-1.446*10^(7);
b5=7.826*10^(3);
b6=-4.417*10;
b7=2.139*10^(-1);
b8=-1.991*10^(4);
b9= 2.778*10^(4);
b10= 9.728*10;
h_w=@(T)(141.355+4202.070*T-0.535*T^(2)+0.004*T^(3)); %J/kg
h_sw=@(T,ws)(h_w(T)-ws*(b1+b2*ws+b3*ws^(2)+b4*ws^(3)+...
    b5*T+b6*T^(2)+b7*T^(3)+b8*ws*T+b9*ws^(2)*T+b10*ws*T^(2)));  %Without the effect of pressure
% hsw(T,p,sw)=hsw(T,po,ws)+v*(p-po) where po=101.325kpa atmospheric
% pressure
h_s=@(T,ws)(h_sw(T,ws)+v_sw(T,ws)*(P-101.325)*1000);   %With the effect of pressure

% Specific Entropy
c1=-4.231*10^(2);
c2=1.463*10^(4);
c3=-9.880*10^(4);
c4=3.095*10^(5);
c5=2.562*10;
c6=-1.443*10^(-1);
c7=5.879*10^(-4);
c8=-6.111*10;
c9=8.041*10;
c10=3.035*10^(-1);
s_w=@(T)(0.1543+15.383*T-2.996*10^(-2)*T^(2)+8.193*10^(-5)*T^(3)-1.370*10^(-7)*T^(4)); %J/kg
s_sw=@(T,ws)(s_w(T)-ws*(c1+c2*ws+c3*ws^(2)+c4*ws^(3)+c5*T+c6*T^(2)+c7*T^(3)+...
    c8*ws*T+c9*ws^(2)*T+c10*ws*T^(2)));     %J/kg

% Chemical Potential
g_sw=@(T,ws) (h_s(T,ws)-(T+273.15)*s_sw(T,ws));

dhsw_ws=@(T,ws) ((-1)*(b1+2*b2*ws+3*b3*(ws^2)+4*b4*(ws^3)+b5*T+b6*(T^2)+b7*(T^3)+2*b8*ws*T+...
                 3*b9*T*(ws^2)+2*b10*ws*(T^2)));
dsw_ws=@(T,ws)  ((-1)*(c1+2*c2*ws+3*c3*ws^(2)+4*c4*ws^(3)+c5*T+c6*T^(2)+c7*(T^3)+2*c8*ws*T+...
                  3*c9*(ws^2)*T+2*c10*ws*(T^2)));
              
dg_sw=@(T,ws) (dhsw_ws(T,ws)-(T+273.15)*dsw_ws(T,ws));
chempot_w=@(T,ws) (g_sw(T,ws)-ws*dg_sw(T,ws));
chempot_s=@(T,ws) (g_sw(T,ws)+(1-ws)*dg_sw(T,ws));
SW_Gibbs(T,uT,S,uS,P,uP);


%% Exergy Calculation using the chemical potential relation from the paper:
% On exergy calculations of seawater with applications in desalination systems

ef1=((h-h_star)-(T0+273.15)*(s-s_star)+...
    S*(chempot_s(T0,S)-chempot_s(T0,S0))+(1-S)*(chempot_w(T0,S)-chempot_w(T0,S0)))/1000;

% ef2=((h_s(T,S)-h_s(T0,S0))-(T0+273.15)*(s_sw(T,S)-s_sw(T0,S0))+...
%     S*(chempot_s(T0,S)-chempot_s(T0,S0))+(1-S)*(chempot_w(T0,S)-chempot_w(T0,S0)))/1000
%% Exergy Calculation using the chemical potential relation from the paper:
% Thermophysical properties of seawater: 
%A review and new correlations that include pressure dependence

% ef2=((h-h_star)-(T0+273.15)*(s-s_star)+...
%     S*(us_star-us_not)+(1-S)*(uw_star-uw_not))/1000

end


