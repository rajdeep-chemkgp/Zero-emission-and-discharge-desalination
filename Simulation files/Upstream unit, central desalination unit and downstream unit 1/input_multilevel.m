% Excel Workbook name: multilevelfactorialconv.xlsx
% Worksheet name: InputData
% MULTILEVEL FACTORIAL matrix range: A3:E752
clear
clc
% filename = 'multilevelfactorialconv.xlsx';
% sheetread = 'InputData';
OutPut=[];
% ABC=xlsread('multilevelfactorialconv.xlsx','InputData');
ABC=xlsread('co2calculation for final ndsCopy.xlsx','InputData');
for i=1:length(ABC(:,1))
%     output=MOO_reactive_middle_ok_REDUCED_FEED_TEMP_final(ABC(i,:));
output=CO2emission(ABC(i,:));
    %output=[TAC,totalbatch,stime,XD1,PTIME,Annualprod];
%     output=[XBBIODIESEL,TAC,CO2total];
    output=[CO2total];
    OutPut=[OutPut;output];   
end
xlswrite( 'co2calculation for final ndsCopy.xlsx',OutPut,'outputdata')
% NT_range = 'A3:A752';
% Q_range='B3:B752';
% P_range='C3:C752';
% WH_range='D3:D752';
% DCOL_range='E3:E752';
% NT_range = xlsread(filename,sheetread, NT_range);
% Q_range = xlsread(filename,sheetread, Q_range);
% P_range = xlsread(filename,sheetread, P_range);
% WH_range = xlsread(filename,sheetread, WH_range);
% DCOL_range = xlsread(filename,sheetread, DCOL_range);
% xlswrite( 'multilevelfactorialconv.xlsx',stime,'outputdata' ,'B2:ABW2') 
% xlswrite( 'multilevelfactorialconv.xlsx',XD(1),'outputdata' ,'B3:ABW3')
% xlswrite( 'multilevelfactorialconv.xlsx',PTIME,'outputdata' ,'B4:ABW4')
% xlswrite( 'multilevelfactorialconv.xlsx',Annualprod,'outputdata' ,'B5:ABW5')
% xlswrite( 'multilevelfactorialconv.xlsx',TAC,'outputdata' ,'B6:ABW6')
% xlswrite( 'multilevelfactorialconv.xlsx',totalbatch,'outputdata' ,'B7:ABW7')