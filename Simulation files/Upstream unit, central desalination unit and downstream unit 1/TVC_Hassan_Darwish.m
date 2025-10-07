function [Mr,Md,Ts]=TVC_Hassan_Darwish(Mms,Pms,Ts,Tn)
% Md is the flowrate of the mixed stream = S
Er=Pms/Psaturation(Tn);
% Cr=Psaturation(Ts)/Psaturation(Tn);
% Cr=Compression_ratio;
Cr=Psaturation(Ts)/Psaturation(Tn);
% Ts=Tsaturation(Cr*Psaturation(Tn))
if Er>=100
    
    Mr=-1.93422581403321+2.152523807931*Cr+113.490932154749/Er...
       -0.522221061154973*(Cr^2)-14735.9653361836/(Er^2)...
       -31.8519701023059*(Cr/Er)+0.047506773195604*(Cr^3)...
       +900786.044551787/(Er^3)-495.581541338594*(Cr/(Er^2))...
       +10.0251265889018*((Cr^2)/Er); 
elseif (10<=Er) && (Er<100)
    
    Mr=-3.20842210618164+3.93335312452389*Cr+27.2360043794853/Er...
       -1.19206948677452*(Cr^2)-141.423288255019/(Er^2)...
       -22.5455184193569*(Cr/Er)+0.125812687624122*(Cr^3)...
       +348.506574704109/(Er^3)+41.7960967174647*(Cr/(Er^2))...
       +6.43992939366982*((Cr^2)/Er);
elseif (2<=Er) && (Er<10)
    
    Mr=-1.61061763080868+11.0331387899116*log(Cr)+13.5281254171601/Er...
       -14.9338191429307*((log(Cr))^2)-34.4397376531113/(Er^2)...
       -48.4767172051364*(log(Cr)/Er)+6.46223679313751*(log(Cr)^3)...
       +29.9699902855834/(Er^3)+70.8113406477665*(log(Cr)/(Er^2))...
       +46.9590107717394*((log(Cr))^2)/Er;
else Er<2
        disp('Wrong Er')
end
Md=(Mms*(1+Mr))/Mr;
Mev=Md-Mms;
%% ------------------------------------------Psaturation function-----------------------------------
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
end