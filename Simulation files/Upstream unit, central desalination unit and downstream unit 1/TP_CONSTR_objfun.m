function [y, cons] = TP_CONSTR_objfun(ABC)
% Objective function : Test problem 'CONSTR'.
%*************************************************************************
In=[ABC]
y = [0,0,0];
cons = [0,0,0];
% [RR,CO2_yearly,tac,PR]=MED_FF(ABC);
[Out]=Call_file(ABC)
% [~,B]=MED_PCF(ABC)
%fprintf('%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n',TIME,TAC,XBBIODIESEL,CO2total,NT,Qt,Ratio,felt1,NF1,NF2)
% [TAC,precovery,totalbatch]=Binary_Batch_Distillation_MW(ABC);
y(1) =(-1)*Out(1); %Distillate flowrate      
y(2) =(-1)*Out(2);
y(3) =Out(3);      %CO2 emission
y(4)= Out(4);      %CO2 emission 
%end

% calculate the constraint violations
c1 = Out(1)-0;
if(c1<0)
    cons(1)=abs(c1);
end
c2=Out(2)-0;
if(c2<0)
    cons(2)=abs(c2);
end
c3=Out(3)-0;
if(c3<0)
    cons(3)=abs(c3);
end
c4=Out(4)-0;
if(c4<0)
    cons(4)=abs(c4);
end

% 
% c = -x(2) + 9*x(1) - 1;
% if(c<0)
%     cons(2) = abs(c);
% end