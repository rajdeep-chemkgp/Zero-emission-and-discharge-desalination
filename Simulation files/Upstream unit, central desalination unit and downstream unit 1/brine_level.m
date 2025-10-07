function [dL,v]=brine_level(L,n,P,rho,v,brine,Abb)
g=9.8;
brine;
Ab=Abb;
for i=1:n
if i==n
P1=((P(i)+rho(i)*g*L(i))/rho(i));
P2=101325/rho(i);
v(i)=(2*(abs(P1-P2)))^0.5;
B(i)=v(i)*rho(i)*Ab;
dL(i,1)=brine(i)-B(i);
else
P1=((P(i)+rho(i)*g*L(i))/rho(i));
P2=((P(i+1)+rho(i+1)*g*L(i+1))/rho(i));
v(i)=(2*(abs(P1-P2)))^0.5;
B(i)=v(i)*rho(i)*Ab;
dL(i,1)=brine(i)-B(i);
end
end
end