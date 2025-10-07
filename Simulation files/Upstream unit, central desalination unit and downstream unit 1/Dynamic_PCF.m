function [y]=Dynamic_PCF(t,c,ABC,Area_evap,Steam,beta,Cv,Tcw,beta_con,Acond,Xf,Ftotal)

Ts=ABC(1);
Tn=ABC(2);
Tf=ABC(3);
n=ABC(4);
% Xf=ABC(5);
% Xn=ABC(6);
Mms=ABC(5);
Pms=ABC(6);
% Ftotal=ABC(7);

Wf=Ftotal/n;
%% --------------------------- Sorting the input variables ---------------------------------------------
Tv=[]; L=[]; X=[];
for i=1:1:n
    Tv=[Tv c(i)];
    L=[L c(i+n)];
    X=[X c(i+2*n)];
end
Lc=c(3*n+1);
% Tcon=c(3*n+2);
x_gkg=10^3*[X];
Xppm=10^6*[X];
x=10^2*[X];

%% -------------------------------------- Physical attributes of the evaporator ------------------------------
% ------------------------Area of Evaporator Effects------------------------------*%
% Width of effects
we=3.6;
% Height of effect
% He=6.8;
He=4;
%Length of effects
Le=6.1;
% Area of effects
% Ace=Le*we;
dia=4.8;
% dia=1;
Ace_evap=(pi/4)*(dia^2);
% Ace=we*Le;
% Ace=18;
%*------------------------Area of condenser------------------------------*%
% Height of the condenser
Hcon=3.2;
% Length of the condenser
Lcond=8;
% Width of the condenser
Wcond=2.2;
% Cross sectional area of the condenser
Ac=Wcond*Lcond;
%Number of condenser tubes
Ncon=1890;
%Diameter of condenser tubes
Dcon=19.05e-3;
%Length of condenser tubes
Lcon=6.5;
%Area of the condenser
% Acon=Ncon*pi*Dcon*Lcon;
% Acon=184.5457;
% Vcon=Ncon*(pi/4)*(Dcon^2)*Lcon;
Vcon=0.31124;

%*------------------------Heat Transfer Area of evaporator------------------------------*%
% Number of evaporator tubes 
Nt=3500;
% Evaporator tube diameter
De=28.575e-3;
% Evaporator tube length
Le=3.9;
% Ae=Nt*pi*De*Le;
% Ae=59.8637;
% Ae=415.1255;
%*------------------------Cross sectional area of brine pipes------------------------------*%
Db=18; %Diameter of the pipe in inches
Ab=((0.0254*Db)^2)*(pi/4);
% ---------------------------------------------------------------------------------------------------
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
for k=1:length(Tv)
Tb(k)=Tv(k)+BPE(Tv(k),x(k));
P(k)=Psaturation(Tv(k))*1000; % in Pa
rho(k)=rhol(Tb(k),Xppm(k));
end
% v=vb;
g=9.8;
% Wb=Brine;
Wb=[];
for i=1:n
    Br(i)=beta(i)*rho(i)*Ab*(L(i))^0.5;
end
% for i=1:n
% if (i==1)
%     d(i)=0;
%     Di(i)=((Steam*Lamda(Ts))-Wf*(hb(Tv(i))-hb(Tf)))/Lamda(Tv(i));
% else
%     d(i)=Br(i-1)*(hv(Tb(i-1))-hv(Tv(i)))/Lamda(Tv(i));
%     Di(i)=(((Di(i-1)+d(i-1))*Lamda(Tv(i-1)))-Wf*(hb(Tv(i))-hb(Tf)))/Lamda(Tv(i))+...
%           Br(i-1)*Cp(Tb(i-1),X(i-1)*1e3)*(Tb(i-1)-Tv(i))/Lamda(Tv(i));
% end
% end
Tcon=(Tf+Tcw)/2;
for i=1:n-1
Di(i)=Cv(i)*(abs(Psaturation(Tv(i))*1000-Psaturation(Tv(i+1))*1000))^0.5;
end
Di(n)=Cv(n)*(abs(Psaturation(Tv(n))*1000-Psaturation(Tcon)*1000))^0.5;
Wb=Br;
Wv=Di;

Ae=Area_evap;
Ace=Ace_evap;
for i=1:1:n
if i==1
k1=Ace*(rhol(Tb(i),Xppm(i))-rhov(Tv(i)));
k2=Ace*L(i)*dpldT(Tb(i),Xppm(i))*(1+dbpedT(Tv(i),x(i)))+...
   Ace*(He-L(i))*dpvdT(Tv(i));
k3=Ace*L(i)*(dpldT(Tb(i),Xppm(i))*dbpedX(Tv(i),x(i))+dpldX(Tb(i),Xppm(i)));
k4=Wf-Wb(i)-Wv(i);

k5=Ace*(hb(Tb(i))*rhol(Tb(i),Xppm(i))-hv(Tv(i))*rhov(Tv(i)));
k6=(Ace*L(i)*rhol(Tb(i),Xppm(i))*dhbdT(Tb(i))*(1+dbpedT(Tv(i),x(i)))+...
   rhov(Tv(i))*Ace*(He-L(i))*dhvdT(Tv(i)))+...
   (Ace*L(i)*hb(Tb(i))*dpldT(Tb(i),Xppm(i))*(1+dbpedT(Tv(i),x(i)))+...
   hv(Tv(i))*Ace*(He-L(i))*dpvdT(Tv(i)));
k7=Ace*L(i)*rhol(Tb(i),Xppm(i))*(dhbdT(Tb(i))*dbpedX(Tv(i),x(i))+dhbdX(Tv(i)))+...
   Ace*L(i)*hb(Tb(i))*(dpldT(Tb(i),Xppm(i))*dbpedX(Tv(i),x(i))+dpldX(Tb(i),Xppm(i)));
k8=Wf*hb(Tf)+Ue(Tv(i))*Ae(1)*(Ts-Tv(i))-Wb(i)*hb(Tb(i))-Wv(i)*hv(Tv(i));

k9=Ace*rhol(Tb(i),Xppm(i))*X(i);
k10=Ace*X(i)*L(i)*dpldT(Tb(i),Xppm(i))*(1+dbpedT(Tv(i),x(i)));
k11=Ace*L(i)*(X(i)*dpldT(Tb(i),Xppm(i))*dbpedX(Tv(i),x(i))+...
    X(i)*dpldX(Tb(i),Xppm(i))+rhol(Tb(i),Xppm(i)));
k12=Wf*Xf-Wb(i)*X(i);

A=k1*k11-k3*k9;
B=k2*k11-k3*k10;
C=k4*k11-k3*k12;
D=k5*k11-k7*k9;
E=k6*k11-k7*k10;
F=k8*k11-k7*k12;

% y(i)=((C*E)-(B*F))/((A*E)-(B*D));
% y(i+n)=((A*F)-(C*D))/((A*E)-(B*D));
% y(i+2*n)=(k12-k9*(y(i))-k10*(y(i+n)))/k11;
y(i)=((A*F)-(C*D))/((A*E)-(B*D));
y(i+n)=((C*E)-(B*F))/((A*E)-(B*D));
y(i+2*n)=(k12-k9*(y(i+n))-k10*(y(i)))/k11;

else
k1=Ace*(rhol(Tb(i),Xppm(i))-rhov(Tv(i)));
k2=Ace*L(i)*dpldT(Tb(i),Xppm(i))*(1+dbpedT(Tv(i),x(i)))+...
   Ace*(He-L(i))*dpvdT(Tv(i));
k3=Ace*L(i)*(dpldT(Tb(i),Xppm(i))*dbpedX(Tv(i),x(i))+dpldX(Tb(i),Xppm(i)));
k4=Wf-Wb(i)-Wv(i)+Wb(i-1);

k5=Ace*(hb(Tb(i))*rhol(Tb(i),Xppm(i))-hv(Tv(i))*rhov(Tv(i)));
k6=(Ace*L(i)*rhol(Tb(i),Xppm(i))*dhbdT(Tb(i))*(1+dbpedT(Tv(i),x(i)))+...
   rhov(Tv(i))*Ace*(He-L(i))*dhvdT(Tv(i)))+...
   (Ace*L(i)*hb(Tb(i))*dpldT(Tb(i),Xppm(i))*(1+dbpedT(Tv(i),x(i)))+...
   hv(Tv(i))*Ace*(He-L(i))*dpvdT(Tv(i)));
k7=Ace*L(i)*rhol(Tb(i),Xppm(i))*(dhbdT(Tb(i))*dbpedX(Tv(i),x(i))+dhbdX(Tv(i)))+...
   Ace*L(i)*hb(Tb(i))*(dpldT(Tb(i),Xppm(i))*dbpedX(Tv(i),x(i))+dpldX(Tb(i),Xppm(i)));
k8=Wf*hb(Tf)+Ue(Tv(i))*Ae(i)*(Tv(i-1)-Tv(i))-Wb(i)*hb(Tb(i))-Wv(i)*hv(Tv(i))+Wb(i-1)*hb(Tb(i-1));

k9=Ace*rhol(Tb(i),Xppm(i))*X(i);
k10=Ace*X(i)*L(i)*dpldT(Tb(i),Xppm(i))*(1+dbpedT(Tv(i),x(i)));
k11=Ace*L(i)*(X(i)*dpldT(Tb(i),Xppm(i))*dbpedX(Tv(i),x(i))+...
    X(i)*dpldX(Tb(i),Xppm(i))+rhol(Tb(i),Xppm(i)));
k12=Wf*Xf-Wb(i)*X(i)+Wb(i-1)*X(i-1);

A=k1*k11-k3*k9;
B=k2*k11-k3*k10;
C=k4*k11-k3*k12;
D=k5*k11-k7*k9;
E=k6*k11-k7*k10;
F=k8*k11-k7*k12;


y(i)=((A*F)-(C*D))/((A*E)-(B*D));
y(i+n)=((C*E)-(B*F))/((A*E)-(B*D));
y(i+2*n)=(k12-k9*(y(i+n))-k10*(y(i)))/k11;    
end

% y=y';
end
%% ---------------------------------Condenser-----------------------------------------
Tcon=(Tf+Tcw)/2;
Dc=Wv(n)-(Steam-Mms);
Dcon=beta_con*(Lc^0.5*1000*Ab);
U=Ucon(Tv(n));
LMTD_con=(Tf-Tcw)/log((Tv(n)-Tcw)/(Tv(n)-Tf));
Qcon=Acond*(U*LMTD_con)+Dcon*4.314*(Tv(n)-Tcon);

d11=Ac*(rhol(Tcon,0)-rhov(Tv(n)));
d12=Ac*Lc*dpldT(Tcon,0)+(Hcon-Lc)*Ac*dpvdT(Tv(n));
d13=Dc-Dcon;
d21=Ac*(rhol(Tcon,0)*hb(Tcon)-rhov(Tv(n))*hv(Tv(n)));
d22=Ac*Lc*hb(Tcon)*dpldT(Tcon,0)+(Hcon-Lc)*hv(Tv(n))*Ac*dpvdT(Tv(n))+...
       Ac*Lc*rhol(Tcon,0)*dhbdT(Tcon)+(Hcon-Lc)*rhov(Tv(n))*Ac*dhvdT(Tv(n));
d23=Dc*hv(Tv(n))-Dc*hb(Tcon)-Qcon;

y(3*n+1)=((d23*d12-d13*d22)/(d21*d12-d11*d22)); %Lcon
% y(3*n+2)=((d13*d21-d23*d11)/(d12*d21-d22*d11)); %Tcon


y=y';

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