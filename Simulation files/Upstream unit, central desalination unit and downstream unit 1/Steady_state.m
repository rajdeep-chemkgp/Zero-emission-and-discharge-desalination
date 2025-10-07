function [A,Steam,beta,Cv,beta_con,Acond]=Steady_state(ABC,Tcw,Xf,Ftotal)


Ts=ABC(1);
Tn=ABC(2);
Tf=ABC(3);
n=ABC(4);
% Xf=ABC(5);
% Xn=ABC(6);
Mms=ABC(5);
Pms=ABC(6);
% Xn=ABC(7);
% Ftotal=ABC(8);




% Ftotal=167;
F=Ftotal/n;
dT=(Ts-Tn)/n;
delA_max=0.1;
delT=dT*ones(1,n);

for i=1:n
Tv(i)=Ts-i*dT;
end


%*------------------------Area of Evaporator Effects------------------------------*%
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
Ace=(pi/4)*(dia^2)*Le;
% Ace=we*Le;
% Ace=18;
%*------------------------Cross sectional area of brine pipes------------------------------*%
Db=18; %Diameter of the pipe in inches
Ab=((0.0254*Db)^2)*(pi/4);
%% ------------------------- Co-relations ------------------------------------------------------------------
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
%% ------------------------------------------------------------------------------------------------------------------------------


% S=Steam;
% for i=1:numel(Tv)
% Tb(i)=Tv(i)+BPE(Tv(i),x(i));
% end

while (delA_max>0.00001)
[~,S]=TVC_Hassan_Darwish(Mms,Pms,Ts,Tn);
Tt=Ts-Tn;
d(1)=0;
D(1)=((S*Lamda(Ts))-F*Cp(Tf,Xf*1e3)*(Tv(1)-Tf))/Lamda(Tv(1));
% D(1)=((S*Lamda(Ts))-F*(hb(Tv(1))-hb(Tf)))/Lamda(Tv(1));
B(1)=F-D(1);
X(1)=(F*Xf)/B(1);
Tb(1)=Tv(1)+BPE(Tv(1),X(1)*10^2);
% fprintf('\nD(1)=%f\tB(1)=%f\tX(1)=%f\tTb(1)=%f\tTv(1)=%f\n',D(1),B(1),X(1),Tb(1),Tv(1));
for i=2:n
% if (i==n)
d(i)=B(i-1)*(hv(Tb(i-1))-hv(Tv(i)))/Lamda(Tv(i));
% D(i)=(((D(i-1)+d(i-1))*Lamda(Tv(i-1)))-F*(hb(Tv(i))-hb(Tf)))/Lamda(Tv(i))+B(i-1)*Cp(Tb(i-1),X(i-1)*1e3)*(Tb(i-1)-Tv(i))/Lamda(Tv(i));
% D(i)=(((D(i-1)+d(i-1))*Lamda(Tv(i-1)))-F*Cp(Tf,Xf*1000)*(Tv(i)-Tf))/Lamda(Tv(i))+B(i-1)*Cp(Tb(i-1),X(i-1)*1e3)*(Tb(i-1)-Tv(i))/Lamda(Tv(i));
D(i)=(((D(i-1)+d(i-1))*Lamda(Tv(i-1)))-F*Cp(Tf,Xf*1000)*(Tv(i)-Tf))/Lamda(Tv(i))+0;
% D(i)=(((D(i-1)+d(i-1))*Lamda(Tv(i-1)))-F*(hb(Tv(i))-hb(Tf)))/Lamda(Tv(i))+(B(i-1)*(hb(Tb(i-1))-hb(Tv(i))))/Lamda(Tv(i));
B(i)=F-D(i)+B(i-1)-d(i);
X(i)=(B(i-1)*X(i-1)+F*Xf)/B(i);
Tb(i)=Tv(i)+BPE(Tv(i),X(i)*10^2);
% fprintf('\nD(%d)=%f\tB(%d)=%f\tX(%d)=%f\tT(%d)=%f\tTv(%d)=%f\n',i,D(i),i,B(i),i,X(i),i,Tb(i),i,Tv(i));
% else
% d(i)=B(i-1)*(hv(Tb(i-1))-hv(Tv(i)))/Lamda(Tv(i));
% D(i)=(((D(i-1)+d(i-1))*Lamda(Tv(i-1)))-F*(hb(Tv(i))-hb(Tf)))/Lamda(Tv(i))+B(i-1)*Cp(Tb(i-1),X(i-1)*1e3)*(Tb(i-1)-Tv(i))/Lamda(Tv(i));
% % D(i)=(((D(i-1)+d(i-1))*Lamda(Tv(i-1)))-F*(hb(Tv(i))-hb(Tf)))/Lamda(Tv(i));
% B(i)=F-D(i)+B(i-1)-d(i);
% X(i)=(B(i-1)*X(i-1)+F*Xf)/B(i);
% Tb(i)=Tv(i)+BPE(Tv(i),X(i)*10^2);
% fprintf('\nD(%d)=%f\tB(%d)=%f\tX(%d)=%f\tT(%d)=%f\tTv(%d)=%f\n',i,D(i),i,B(i),i,X(i),i,Tb(i),i,Tv(i));
% end
end
A(1)=(S*Lamda(Ts))/(Ue(Tv(1))*delT(1));
for j=2:n
A(j)=((D(j-1)+d(j-1))*Lamda(Tv(j-1)))/(Ue(Tv(j))*delT(j));
delA(j)=abs(A(j)-A(j-1));
end
delA_max=max(delA);
Am=sum(A)/n;
%new temperature profile
for i=1:n
delT_adj(i)=delT(i).*(A(i)/Am);
end
delT_sum=sum(delT_adj);
for i=1:n
delT(i)=delT_adj(i)*(Tt/delT_sum);
end
Q=D>0;
W=B>0;
if (all(Q)&&all(W))
    continue
else
    break
end
end


Tcw=31.5;
Tcon=(Tcw+Tf)/2;
Dev=S-Mms;
Dcon=D(n)-Dev;
Qcon=Dcon*Lamda(Tv(n))+Dcon*4.314*(Tv(n)-Tcon);
% Qcon=Dcon*Lamda(Tv(n));
LMTD_con=(Tf-Tcw)/log((Tv(n)-Tcw)/(Tv(n)-Tf));
U=Ucon(Tv(n));
Acond=Qcon/(U*LMTD_con);

Dt=sum(D)+sum(d);
GOR=(Dt)/Mms;
Steam=S;
x_gkg=10^3*[X];
Xppm=10^6*[X];
x=10^2*[X];
% velocity(rho,brine,Area_of_brinepipe)
for k=1:n
v(k)=velocity(rhol(Tb(k),X(k)),B(k),Ab);
P(k)=Psaturation(Tv(k))*1000; %in kPa
rho(k)=rhol(Tb(k),Xppm(k));
end
vb=v;

for l=1:n
if (l==1)
    M(l)=B(l);
else
    M(l)=F+B(l-1);
end
    hs(l)=M(l)/(rho(l)*Ab);
    beta(l)=B(l)/((hs(l))^0.5*rho(l)*Ab);
end
Tcon=(Tf+Tcw)/2;
hcon=Dcon/(1000*Ab);
beta_con=Dcon/(hcon^0.5*1000*Ab);
for k=1:n-1
    Cv(k)=D(k)/((abs(Psaturation(Tv(k))-Psaturation(Tv(k+1)))*1000)^0.5);
end
Cv(n)=D(n)/((abs(Psaturation(Tv(n))-Psaturation(Tcon))*1000)^0.5);
% Cv(n)=D(n)/((abs(Psaturation(Tv(n))-Psaturation(Ts))*1000)^0.5);

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