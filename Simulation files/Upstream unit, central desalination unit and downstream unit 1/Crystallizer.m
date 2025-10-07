function c=Crystallizer(Xbn,mbn,Tbn,Ts,Tf,Tcw)

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
%% ------------------------------------------Crystallizer----------------------------------------------------------------
% A zero liquid discharge system integrating multi-effect distillation and evaporative crystallization for desalination brine treatment
% DOI: https://doi.org/10.1016/j.desal.2020.114928
Rc=0.9;                     % Crystallization ratio
% Xbn=0.053;                  % Brine salinity from last effect
% mbn=110;                    % Brine flowrate from last effect
ms=mbn*Xbn;                 % Massflowrate of salt (Extracting all the salt from the brine)
Dc=mbn*(1-Xbn);             % Distillate flowrate from crystallizer
mcb=10*mbn;                 % Assuming that the recirculated brine flowrate is 10 times the feed flowrate
mcf=mbn+mcb;                % Mixing of feed and recirculation streams
% Tbn=45.3;                   % Temperature of brine from last effect of MED
TH=Ts;                          % Temperature of the heating source (value taken from paper)
Tcb=TH;                         % Assuming that the temperature of the crystallizer is equal to the temperature of the heating source

UA=25000/1000;              % Converting W/K to kW/K since Cp is in kJ/kgK 
Xs=@(T) (0.2628+(62.75e-6)*T+(1.08e-6)*T*T);

Xcf=(ms/((mcf-Dc)*Rc))+Xs(Tcb);  % Salinity of input stream in crystallizer
Xcb=(mcf*Xcf-mbn*Xbn)/mcb;       % Salinity of recirculating stream
Tcf0=(mbn*Tbn+mcb*Tcb)/mcf;      % Temperature of the recirculated stream after mixing mbn and mcb
Tcf=TH-(TH-Tcf0)*exp(-(UA/(mcf*Cp(Tcf0,Xcf*10^3))));        % Temperature of the brine after getting heated        
Qcrystallizer=mcf*Cp(Tcf0,Xcf*10^3)*(Tcf-Tcf0);             % Heat utility of the crystallizer
mhs=Qcrystallizer/Lamda(TH);                                % Flowrate of heat source
Heat=mhs*Lamda(TH);
Xcf0=((mcf-Dc)/mcf)*Xcf;                                    % Circulating stream before heating
Tvcrystal=Tcb-BPE(Tcb,Xcf*100);                             % Temperature of the vapors leaving crystallizer

% It is assumed that for cooling, the seawater feed is utilized. The inlet
% feed conditions are Tcw=31.5 oC and Xf=0.035. In the present scenario,
% Tf= 41.5. Assuming that only 5 oC is achieved across the heat exchanger.
Tcw=31.5;
Xf=0.035;
delT=Tf-Tcw; % This is modelled for MED-TVC-MVC
mcw=Dc*Lamda(Tvcrystal)/(Cp(Tcw,Xf*1000)*delT);

Qcondenser=Dc*Lamda(Tvcrystal);
LMTDcryst=((Tvcrystal-Tf)-(Tvcrystal-Tcw))/log((Tvcrystal-Tf)/(Tvcrystal-Tcw));
Acrystallizer=Qcondenser/(Ucon(Tvcrystal)*LMTDcryst);

delPcf=1.5*10^5;              % bar to Pa; Pressure drop for circulation of brine
delPcw=2*10^5;                % bar to Pa; Pressure drop for cooling water
delPdc=0.5*10^5;              % bar to Pa; Pressure drop for distillate extraction
eff_pump=0.75;
Electricity=((mcf*delPcf)/(rhol(Tcf,Xcf*10^6)*eff_pump)+(mcw*delPcw)/(rhol(Tcw,Xf*10^6)*eff_pump)+...
             (Dc*delPdc)/(1000*eff_pump))/1000;    % in kW

c=[ms Qcrystallizer Electricity Dc mhs mcw Tvcrystal Acrystallizer];

