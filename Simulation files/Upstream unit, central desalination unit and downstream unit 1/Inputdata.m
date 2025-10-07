%function [output]=Inputdata()
% Excel Workbook name: optimumSelection.xlsx
% Worksheet name: InputData
% Decision making matrix range: C4:D103
clear
clc
filename = 'SupplierSelection.xlsx';
sheetread = 'InputData';
decisionMakingMatrix_range = 'C4:D355';
criteriaSign_range= 'C356:D356';
lambdaWeight_w1='C357:D357'; 
output=[filename,sheetread,decisionMakingMatrix_range,criteriaSign_range,lambdaWeight_w1];
decisionMakingMatrix = xlsread(filename,sheetread, decisionMakingMatrix_range);
criteriaSign = xlsread(filename, sheetread, criteriaSign_range);
lambdaWeight = xlsread(filename, sheetread,lambdaWeight_w1);
%end



 








