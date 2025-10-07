function [Out,c]=Call_file(ABC,Tcw,Xf)
Ts=ABC(1);
Tn=ABC(2);
Tf=ABC(3);
n=ABC(4);
Mms=ABC(5);
Pms=ABC(6);
Xn=ABC(7);
% Dtotal=ABC(8);
Ftotal=600;
Tcw=31.5;
Xf=0.035;
dT=(Ts-Tn)/n;
%% ---------------------------------------------- Co-relations----------------------------------------------
%% ----------------------------- Co-relations -----------------------------------------------------------
rhol=@(T,x) (1*10^3*(((4.032219*0.5+0.115313*((0.002*x-150)/150)+(3.26*10^-4)*(2*(((0.002*x-150)/150))^2-1))*(0.5))+...
            (-0.108199*0.5+(1.571*10^-3)*((0.002*x-150)/150)-(4.23*10^-4)*(2*(((0.002*x-150)/150))^2-1))*((2*T-200)/160)+...
            (-0.012247*0.5+(1.74*10^-3)*((0.002*x-150)/150)-(9*10^-6)*(2*(((0.002*x-150)/150))^2-1))*(2*(((2*T-200)/160)^2)-1)+...
            ((6.92*10^-4)*0.5-(8.7*10^-5)*((0.002*x-150)/150)-(5.3*10^-5)*(2*(((0.002*x-150)/150))^2-1))*...
            (4*(((2*T-200)/160)^3)-3*(2*T-200)/160))); %Density of liquid kg/m3
% x is in Ppm             
rhov=@(T) (0.005059+0.00023748*T+(1.777*10^-5)*T^2-(4.327*10^-8)*T^3+(4.342*10^-9)*T^4); %Desnity of vapour kg/m3
%--------------------------------------------------------- HTC Correlations ---------------------------------------------------------%
Ucon=@(T)(1.7194+(3.2063*10^(-2))*T+(-1.5971*10^(-5))*(T*T)-(1.9918*10^(-7))*(T*T*T)); %HTC of condenser @Tv kW/m2-deg C

Ue=@(T)(10^(-3)*(1939.4+1.40562*T-0.0207525*T*T+0.0023186*T^3)); %HTC of evaporator @Tb kW/m2-deg C

% Ucon=@(T) (1.7194+(3.2063*10^-2)*T-(1.5971*10^-5)*T^2+(1.9918*10^-7)*T^3);  %HTC of condenser @Tv kW/m2-deg C
% 
% Ue=@(T) (1.9695+(1.2057*10^-2)*T-(8.5989*10^-5)*T^2+(2.565*10^-7)*T^3); %HTC of evaporator @Tb kW/m2-deg C
%----------------------------------------------Specific heat capacity at constant pressure Correlation-----------------------------------%
Cp=@(T,x)(1*(10^-3)*(((4206.8-6.6197*x+(1.2288*10^-2)*(x^2))+(-1.1262+(5.4178*10^-2)*x-(2.2719*10^-4)*x^2)*T)+...
          (((1.2026*10^-2)-(5.3566*10^-4)*x+(1.8906*10^-6)*x^2)*T^2)+...
        (((6.8777*10^-7)+(1.517*10^-6)*x-(4.4268*10^-9)*(x^2))*T^3)));
% x is in gm/kg
%-----------------------------------------------Enthalpy Correlations-----------------------------------------------------------------%
hb=@(T) (-0.33635409+4.207557011*T-(6.200339*10^-4)*T^2+(4.459374*10^-6)*T^3); % Enthalpy of saturated water kJ/kg
hv=@(T) (2501.689845+1.806916015*T+(5.087717*10^-4)*T^2-(1.221*10^-5)*T^3);    % Enthalpy of saturated water vapor kJ/kg
%--------------------------------------Latent heat of vaporisation Correlation-------------------------------------------------------------%
Lamda=@(T) (2501.897149-2.407064037*T+(1.192217*10^-3)*T^2-(1.5863*10^-5)*T^3); % Latent heat of vaporisation kJ/kg

%--------------------------------------BPE Correlation-----------------------------------------------------------------------------------%
BPE=@(T,x) (x*((8.325*10^-2+1.883*10^-4*T+4.02*10^-6*T^2)+(-7.625*10^-4+9.02*10^-5*T-5.2*10^-7*T^2)*x+...
            (1.522*10^-4-3*10^-6*T-3*10^-8*T^2)*x^2));
% x is in salt percentage.
%-----------------------------------LMTD relation--------------------------------------- %
lmtde=@(Tvi_1,Tfd,Tbi) ((Tvi_1-Tfd)-(Tvi_1-Tbi))/log((Tvi_1-Tfd)/(Tvi_1-Tbi));
lmtdc=@(Tvn,Tfd,Tcond,Tcw)((Tvn-Tf)-(Tcond-Tcw))/log((Tvn-Tf)/(Tcond-Tcw));
%-----------------------------------Specific entropy of --------------------------------------- %
sf=@(T) (-0.00057846+0.015297489*T-2.63129*10^(-5)*T^(2)+4.11959*10^(-8)*T^(3)); % Specific entropy of saturated liquid kJ/kg
sv=@(T) (9.149505306-2.581012*10^(-2)*T+9.625687*10^(-5)*T^2-1.786615*10^(-7)*T^(3));%Specific entropy of saturated vapour kJ/kg (Simplified Equations for Saturated Steam Properties for Simulation Purpose)
%% -----------------------------------------------------------------------------------------------------------
%------------------------------Partial Differentiation--------------------%
dpldT=@(T,x)(0.00000026183333333333331897907481078391*x - 0.010575000000000000903183777767325*(0.000013333333333333333333333333333333*x - 1.0)^2 - 1000.0*(0.000625*T - 0.0625)*(0.000018000000000000000456015433747403*(0.000013333333333333333333333333333333*x - 1.0)^2 - 0.0000000232*x + 0.0078544999999999996029227721308641) - 1000.0*(0.15*(0.0125*T - 1.25)^2 - 0.0375)*(0.0000000011599999999999998260997539449117*x + 0.00010599999999999998800351980188239*(0.000013333333333333333333333333333333*x - 1.0)^2 - 0.00048599999999999998859549418805948) - 0.69059375000000002039319234314885);
dbpedT=@(T,x)(x*(0.0000080399999999999992550099928156904*T - 1.0*x*(0.0000010399999999999999952896836616367*T - 0.000090199999999999983424543714694011) - 1.0*x^2*(0.00000006000000000000000787279855023193*T + 0.0000030000000000000000760025722912339) + 0.00018829999999999999624362978511982));
dpvdT=@(T)(0.000000017368000000000000408381661820556*T^3 - 0.00000012981000000000000895951688490923*T^2 + 0.000035539999999999995151846871044299*T + 0.00023748000000000000987883386205368);
dbpedX=@(T,x)(0.00018829999999999999624362978511982*T - 1.0*x*(0.00000051999999999999999764484183081836*T^2 - 0.000090199999999999983424543714694011*T + 0.0007625) - 1.0*x*(0.00000051999999999999999764484183081836*T^2 - 0.000090199999999999983424543714694011*T + 2.0*x*(0.000000030000000000000003936399275115965*T^2 + 0.0000030000000000000000760025722912339*T - 0.00015220000000000001322969511718952) + 0.0007625) - 1.0*x^2*(0.000000030000000000000003936399275115965*T^2 + 0.0000030000000000000000760025722912339*T - 0.00015220000000000001322969511718952) + 0.0000040199999999999996275049964078452*T^2 + 0.08325);
dpldX=@(T,x)(0.00000000011591111111111111444918167181742*x - 1000.0*(0.0125*T - 1.25)*(0.00000000000030080000000000002569056078982612*x - 0.000000043506666666666667445118044099672) - 1000.0*(0.0000000000000064000000000000001621388208879656*x - 0.000000023680000000000000012160411566597)*(2.0*(0.0125*T - 1.25)^2 - 1.0) - 1000.0*(0.000000000000037688888888888884623473707335961*x - 0.0000000016666666666666665206607741052854)*(4.0*(0.0125*T - 1.25)^3 - 0.0375*T + 3.75) + 0.00076005999999999999198104413030327);
dhbdT=@(T)(0.00001337812200000000174371448530275*T^2 - 0.0012400677999999999180646970131647*T + 4.2075570109999995693783603201155);
dhvdT=@(T)(- 0.000036629999999999999962373847806063*T^2 + 0.0010175434000000000128033361690427*T + 1.8069160150000000975012426351896);
dhbdX=@(T)(0);
%--------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------------------


%% -------------------------------- Steady state code solution -------------------------------------------------
[Area_evap,Steam,beta,Cv,beta_con,Acond]=Steady_state(ABC,Tcw,Xf,Ftotal);

%% -------------------------------- Dynamic code solution -------------------------------------------------
tspan=1:1:4000;
for i=1:n
Temp(i)=Ts-i*dT;
end
C0=[];
% Temp=40*ones(1,n);
Lb=0.1*ones(1,n);
% Lb=L;
% Temp=Tv;
% Xb=X;
% Tcon=(Tf+Tcw)/2;
Lcon=0.1;
Xb=Xf*ones(1,n);
C0=[Temp Lb Xb Lcon];

[t,c]=ode15s(@Dynamic_PCF,tspan,C0,[],ABC,Area_evap,Steam,beta,Cv,Tcw,beta_con,Acond,Xf,Ftotal);
% [~,Brine,Distillate]=Dynamic_PCF(tspan,C0,ABC,Area_evap,Ftotal,Steam,beta);
%% -------------------------- Steady state solutions from dynamic code -------------------------------------
Tv=[]; L=[]; X=[];
for i=1:1:n
    Tv=[Tv c(:,i)];
    L=[L c(:,i+n)];
    X=[X c(:,i+2*n)];
end

%% ------------------------- Calculating brine and vapour flowrate at steady state --------------------------

dim_c=size(c); % Dimension of output vector c
row_c=dim_c(1);% Row of output vector c
Cs=c(row_c,:); % Extracting the steady state (last row) of output vector c
Tvs=Cs(1:n);
Ls=Cs(n+1:2*n);
Xs=Cs(2*n+1:3*n);

Wf=Ftotal/n;
x_gkg=10^3*[Xs];
Xppm=10^6*[Xs];
x=10^2*[Xs];

for k=1:n
Tbs(k)=Tvs(k)+BPE(Tvs(k),x(k));
P(k)=Psaturation(Tvs(k))*1000; % in Pa
rho(k)=rhol(Tbs(k),Xppm(k));
end

Db=18; %Diameter of the brine pipe in inches (Mazini et al., 2014)
Ab=((0.0254*Db)^2)*(pi/4); % Cross sectional area of brine pipes

for i=1:n
    Br(i)=beta(i)*rho(i)*Ab*(Ls(i))^0.5;
end

Tcon=(Tf+Tcw)/2;
for i=1:n-1
Di(i)=Cv(i)*(abs(Psaturation(Tvs(i))*1000-Psaturation(Tvs(i+1))*1000))^0.5;
end
Di(n)=Cv(n)*(abs(Psaturation(Tvs(n))*1000-Psaturation(Tcon)*1000))^0.5;
for j=1:n

% Dtotal=Distillate;
Devaporator=sum(Di);
Brine=Br(n);
Dev=Steam-Mms;
Dcon=Di(n)-Dev;
Qcon=Dcon*Lamda(Tvs(n))+Dcon*4.314*(Tvs(n)-Tcw);
% LMTD_con=(Tf-Tcw)/log((Tvs(n)-Tcw)/(Tvs(n)-Tf));
% U=Ucon(Tvs(n));
% Acond=Qcon/(U*LMTD_con);
Mcw=Qcon/(Cp(Tf,Xf*1000)*(Tf-Tcw))
Msw=Mcw+Ftotal
Daspen=694.319941804669;
Dtotal=Devaporator+Daspen
Recovery=Dtotal/Msw

%% -------------------------------------- CO2 emission -------------------------------------------------------------
%Nano-catalytic heterogeneous reactive distillation for algal biodiesel production: Multi-objective optimization and heat integration
Tftb=1800;                                      %Flame temperature of the boiler (deg C)
Ta=Tcw;                                          %Ambient temperature (deg C)
Tstack=160;                                     %Stack temperature (deg C)
alpha=3.67;                                     %Ratio of molar masses of CO2 to C
NHV=5.61e4;                                     %Net heating value of the fuel (kJ/kg)
C_percent=75.38;                                %Natural gas carbon percentage
%----------------------- MED-PCF-TVC---------------------------------------
Qproc=Mms*Lamda(Tsaturation(Pms))
% Qproc=1681736.901
% Qproc=S*Lamda(Ts);
Tms_TVC=Tsaturation(Pms);
Lamda_proc_TVC=Lamda(Tsaturation(Pms));             %Latent heat of steam supplied to the process
Hproc_TVC=hv(Tsaturation(Pms));                     %Enthalpy of steam supplied to the process
Hwater_TVC=hb(Tvs(1));                           %Enthalpy of water to the boiler feed
% Hwater=hb(Tv(1));
Q_flash_aspen=1461449.76; % kJ/s
Q_crystallizer_aspen=91981.7701; % kJ/s
Qtotal=Qproc+Q_flash_aspen+Q_crystallizer_aspen
% Qfuel_evap=((Qproc/Lamda_proc_TVC)*(Hproc_TVC-Hwater_TVC))*((Tftb-Ta)/(Tftb-Tstack))
Qfuel_evap=((Qtotal/Lamda_proc_TVC)*(Hproc_TVC-Hwater_TVC))*((Tftb-Ta)/(Tftb-Tstack))
Mms_total=Qtotal/Lamda(Tsaturation(Pms))
% ------------------------- Estimating the electricity to be derived from turbine -------------------------------------------------
Tinlet=1027;
Toutlet=720;
To=Tcw;
nc=(Tinlet-Toutlet)/(Tinlet+273);
ngt=(Toutlet-Tstack)/(Toutlet-To);
fdelP=3571;     %Pressure losses
rho_d=1000;     %density of distilled water
mu=0.75;        %efficiency of power generation and pumps
Eevaporator=(fdelP*Devaporator)/(rho_d*mu)
E_pump1_aspen=6.82947325
E_pump2_aspen=32.330615300000005
Etotal=Eevaporator+E_pump1_aspen+E_pump2_aspen
% Eevaporator=7993.494596
% Qfuel_turbine=(Eevaporator/ngt)*(1/(1-nc));
Qfuel_turbine=(Etotal/ngt)*(1/(1-nc))
% Qfuel_turbine=0;
CO2_emission=((Qfuel_evap+Qfuel_turbine)/NHV)*(C_percent/100)*alpha %kg/sec
% CO2_emission_per_distillate=CO2_emission/Dtotal
%% --------------------------------- Solar calculation --------------------------------------------------------
% Qprocess=Qtotal; %kJ/s
% [collector_CC, collector_OC,TES_cost, BHX_CC,BHX_OC,Asolar,q_ra0] = PTSC_Solar_TVC(Qprocess,Mms,Pms,Tvs(1))
%% --------------------------------- Open raceway pond --------------------------------------------------------
THY=8760;       %Total hours per year (hr/yr)
CO2_emission_solar=CO2_emission;
[Pond_area,NS_growth,Total_cost_ORP,SP_lipids,energy_lipid_extraction_kW,cost_lipid_extraction_per_year]=Open_raceway_pond(CO2_emission_solar,THY)
%% --------------------------------- Cost of PV --------------------------------------------------------
% Total_electrical=Etotal; %kJ/s
% [TAC_PV,Area_PV]=Cost_PV(Total_electrical);
% Capital_PV=Cost(1)
% Operating_PV=Cost(2)
% %% --------------------------------- Evaporation pond --------------------------------------------------------
% life=25;        %life of plant (yr)
% Mcw_total=Mcw_crystallizer+Mcw;
% X=(Mcw_crystallizer*Xf+Mcw*Xf)/Mcw_total;
% [Area_of_evaporation_pond,Capital_cost_evap_pond,Operating_cost_evap_pond,Msalt_pond]=Evaporation_pond(Xn,Mcw_total,life,Tcw,Tf,X)
%% -------------------------------------- FWPC -------------------------------------------------------------
% A multi-objective optimisation framework for MED-TVC seawater desalination process based on particle swarm optimisation
Cmat_MED=3644;  %Material of MED ($/m2) 
Ir=0.07;        %Interest rate
Cmat_Cond=500;  %Material of Condenser ($/m2) 
fdelP=3571;     %Pressure losses
mu=0.75;        %efficiency of power generation and pumps
life=25;        %life of plant (yr)
Csteam=0.004;   %Cost of steam ($/kg)
Kmed=1.4;       %Coefficient for MED
Clab=0.05;      %Coefficient for labour ($/m3)
THY=8760;       %Total hours per year (hr/yr)
Cchem=0.024;    %Cost of chemical treatment ($/m3)
Cpow=0.09;      %Cost of power ($/kWh)
% Cpow=0.049;
Kintake=50;     %Seawater intake ($day/m3)
Kcond=2.8;      %Coefficient for condenser
rho_d=1000;     %density of distilled water
%----------------------------------- Total capital cost--------------------
Ctvc=(7912*Dev*((Tvs(n)/Psaturation(Tvs(n))))^0.005)*(Pms^(+0.75));
Cmed=Kmed*Cmat_MED*((n*Area_evap(1))^0.64);
Ccondenser=Kcond*Cmat_Cond*(Acond^0.8);
Cintake=(Kintake*24*3600*Msw)/rhol(Tf,10^6*Xf);
Cap_aspen=42127500;
Capex_equipment=Cintake+Cmed+Ccondenser;
Capex_civilwork=0.15*Capex_equipment;
Capex_direct=Capex_equipment+Capex_civilwork
Capex_indirect=0.25*Capex_direct
Tcc=Capex_direct+Capex_indirect+Cap_aspen        %Total capital cost ($)
% ------------------------ Csteam calculation -----------------------------
if (Pms<=101.325)
    Csteam_TVC=Csteam;
else 
    Csteam_TVC=Csteam*(0.05*Pms*0.001+0.95);
end
% if (Psaturation (Ts_MVC)<=101.325)
%     Csteam_MVC=0.004;
% else 
%     Csteam_MVC=0.004*(0.05*Pms*0.001+0.95);
% 
% end
%----------------------------------- Total operating cost--------------------
AOClab=(Clab*THY*3600*Devaporator)/rho_d                   %Cost of human labour
AOCman=0.002*Tcc                                      %Cost of manutenation
Tms=Tsaturation(Pms);
AOCsteam_TVC=((Csteam_TVC*THY*(Tsaturation(Pms)-40)*Mms)/80)    %Cost of external steam for TVC
% AOCsteam_TVC=0
% AOC_pumps=((Cpow*THY*1)/(rho_d*mu))*fdelP*Devaporator;       %Cost of pumps
AOC_pumps=((Cpow*THY*1)/(rho_d*mu))*fdelP*Dtotal;       %Cost of pumps
AOCpow=AOC_pumps
AOCchem=(Cchem*THY*3600*Devaporator)/rho_d         %Annual operating cost of chemicals
Cooling_water=12577.4; %kg/s from aspen
AOC_ASPEN=16859500;
AOC=AOCchem+AOClab+AOCpow+AOCman+AOCsteam_TVC+AOC_ASPEN
%----------------------------------- Total Annual cost--------------------

Capital_recovery_factor=(Ir*((1+Ir)^life))/(((1+Ir)^life)-1)   %$(1/yr)
Salt_SP=(46.3988593287717*3600*THY)*0.03;
% TAC=AOC+Tcc*Capital_recovery_factor+Total_cost_ORP+cost_lipid_extraction_per_year-Salt_SP-SP_lipids        %Total annual cost ($/yr)
TAC=AOC+Tcc*Capital_recovery_factor        %Total annual cost ($/yr)
% FWC=TAC*((86.4*sum(Wv)*0.9*365)^(-1))     %Freshwater cost in $/m3
FWC=TAC*(Dtotal*THY*3600)^(-1);          %Freshwater cost in $/kg
tac=FWC*1000                               %Freshwater cost in $/m^3

%% ------------------------------------------- Exergy efficiency ----------------------------------------------
% Assuming the efficiency of pump is 75% and the maximum energy wasted in
% pumps is 2kWh (The feasibility of integrating ME-TVC+MEE with Azzour
% South Power Plant: Economic evaluation)
% With this we are calculating the outlet pressure for MSW_in (seawater
% flowing in from the sea using pump)

% Tdistillate=(Devaporator*Tcon+Dcrystallizer*Tvcrystallizer)/Distillate;
Tdistillate=(Devaporator*(Tcw+273)+Daspen*360.194340092048)/Dtotal;
Exergy_steam_in=Lamda(Tsaturation(Pms))-(Tcw+273)*Lamda(Tsaturation(Pms))/(Tsaturation(Pms)+273);
Exergy_steam_out=Lamda(Tvs(1))-(Tcw+273)*Lamda(Tvs(1))/(Tvs(1)+273);
% Exergy_distillate_out=exergy_calculation(101.325,Tcon,0,Tcw,Xf);
Exergy_brine_out=exergy_calculation(101.325,Tcw,Xs(n),Tcw,Xf);
Exergy_seawater_in=exergy_calculation(101.325,Tcw,Xf,Tcw,Xf);
Exergy_seawater_out=exergy_calculation(105,Tcw,Xf,Tcw,Xf);
Exergy_coolseawater_out=exergy_calculation(105,Tf,Xf,Tcw,Xf);
Exergy_brine_in=exergy_calculation(Psaturation(Tvs(n)),Tvs(n),Xs(n),Tcw,Xf);
Exergy_distillate_in=exergy_calculation(Psaturation(Tcon),Tcon,0,Tcw,Xf);
% E_in=Mms*(Exergy_steam_in-Exergy_steam_out)+(1/0.75)*(Msw*(Exergy_seawater_out-Exergy_seawater_in)+...
%      Distillate*(Exergy_distillate_out-Exergy_distillate_in)+...
%      Brine*(Exergy_brine_out-Exergy_brine_in));
Exergy_distillate_out=exergy_calculation(101.325,Tdistillate,0,Tcw,Xf);
% Tambient=Tcw+273.15;
% Tsun=5770;  % Temperature of sun (K)
% Exergy_solar_in=Qtotal*(1-(4/3)*(Tambient/Tsun)+(1/3)*((Tambient/Tsun)^4)); %kW
% ST=q_ra0;
% Apv=Area_PV;
% Exergy_PV_in=(1-Tambient/Tsun)*(ST*(10^(-3)))*Apv; %kW
Tsteam=Tsaturation(Pms)+273.15;
% % Thcrystallizer=Ts+273.15;
Tout=Tcw+273;
% Q=(Mms_total*Lamda(Tsaturation(Pms)))*(1-Tout/Tsteam);
% E_in=Q+(1/0.75)*(Msw*(Exergy_seawater_out-Exergy_seawater_in)+...
%      Distillate*(Exergy_distillate_out-Exergy_distillate_in))+Ecrystallizer+Eevaporator;
% Tot_Brine=Brine+Mcw;
Tot_Brine=46.3988593287717;
Xsalt=0.0624707463713841;
N=(Tot_Brine*THY*3600)/(Xsalt*58.44+(1-Xs(n))*18.01528);    %kmoles
R=8.31447;    % kJ/kmolK
T=325.624310430863; %K
moles_of_salt=((Tot_Brine*THY*3600)*Xsalt)/58.44;
moles_of_water=((Tot_Brine*THY*3600)*(1-Xsalt))/18.01528;
xa=moles_of_water/(moles_of_water+moles_of_salt);
xb=moles_of_salt/(moles_of_water+moles_of_salt);
Wsep_total=-N*R*T*(xa*log(xa)+xb*log(xb)); % kJoules
Wsep=Wsep_total/(THY*3600); % kJ/s
% xa=water; xb=salt; n=number of moles of solution; R=gas constant; 
% Q=Mms_total*(Exergy_steam_in-Exergy_steam_out);
% Q=Qtotal*(1-Tout/Tsteam);
Q=Mms_total*(1-Tout/Tsteam)+Qtotal
E_in=Q+Etotal;

% Wmin=Distillate*Exergy_distillate_out+Brine*Exergy_brine_out+Mcw*Exergy_coolseawater_out;

Wmin=Dtotal*Exergy_distillate_out+Wsep;
Efficiency=(Wmin/(E_in))*100;

%% ----------------------------------------------------------------------------------------
% A=[Efficiency Devaporator CO2_emission tac];
A=[Efficiency Dtotal CO2_emission tac]; %Solar
k=length(A);
E=isnan(A);
R=any(E);
T_check=abs(((Tn-Tvs(n))/Tn));
X_check=abs(((Xn-Xs(n))/Xn));
Balance_check=Ftotal-(Devaporator+Brine);
if (isreal(c))&&(isreal(A)&&(Devaporator>0)&&(CO2_emission>0)&&(tac>0)&&R==0&&T_check<10^(-1)&&Devaporator<0.6*Ftotal&&...
        Balance_check<=0.01&&Efficiency<100)
% Out=[CO2_yearly TAC GOR];
Out=A;
else
Out=(-1)*ones(1,length(A));
end

end
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