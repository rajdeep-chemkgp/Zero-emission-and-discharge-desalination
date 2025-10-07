function [ cc ] = topsis (decisionMakingMatrix,lambdaWeight,criteriaSign )
 %-Technique for Order of Preference by Similarity to Ideal Solution

 %-This function implements TOPSIS method with Information entropy
 % weighting Method
% - criteriaSign is a vector specifying whether an criterion has to be maximized
% or minized . +1 is for positive criterion and -1 for negetive criterion
% -Before executing the function you have to define "decisionMakingMatrix"
%%%%%%%%%%%%%%%%%%%%%%%
sumDmm=sum(decisionMakingMatrix());
sumDmmMatrix=repmat(sumDmm,size(decisionMakingMatrix,1),1);
pij=decisionMakingMatrix./sumDmmMatrix; %Normalizing Decision Making Matrix  
% you can use this code if you like: pij=decisionMakingMatrix./repmat(clsm,size(evl,1),1)
lnm= -1 / log(size(decisionMakingMatrix,1)); 
lnNormDmm = log(pij);
%Step 3: Calculate weight of criteria by entropy technique
E=lnm .* sum(pij .* lnNormDmm);
dj=ones(1,size(E,2))-E ;%Calculating the information Entropy of Criterion j 
weightEntropy=dj ./sum(dj) ;% computing Entropy weight
wt=lambdaWeight .*weightEntropy ./sum(lambdaWeight .*weightEntropy); 
sqrtxij=sqrt(sum(decisionMakingMatrix().^2)) ;
%Step 4-Construct a normalized decision matrix:
N =decisionMakingMatrix./repmat(sqrtxij,[size(decisionMakingMatrix,1) 1]);										    
Wj=eye(size(wt,2)) .* repmat(wt.*criteriaSign,size(wt,2),1);
%Step 5: Construct the weighted normalized decision matrix by building the diagonal matrix
V=N*Wj;	
%Step 6: Compute the positive ideal solution (PIS) A+ and the negative ideal solution (NIS) A? of the alternatives:
Apositive=max(V); 
Anegetive= min(V);
ApositivMtrix=repmat(Apositive,size(V,1),1); 
AnegetiveMtrix=repmat(Anegetive,size(V,1),1);
s1=(V-ApositivMtrix).^2;
s2=(V-AnegetiveMtrix).^2;
for (j=1:1:size(s1,1))
sumAPositive(j)=sum(s1(j,:));
end
for (j=1:1:size(s2,1))
sumANegetive(j)=sum(s2(j,:));
end
%Step 7: Compute the distance of each alternative from PIS(dPositive) and
%NIS (dNegative)
dPositive=sqrt(sumAPositive); 
dNegetive=sqrt(sumANegetive);
sumD=dNegetive+dPositive;
%Step 8: Compute the closeness coefficient of each alternative:
cc=dNegetive./sumD;   
xlswrite( 'SupplierSelection.xlsx',cc,'OutputData CC' ,'B3:MO3') 
xlswrite( 'SupplierSelection.xlsx',N,'OutputForN' ,'B2:C353') 
xlswrite( 'SupplierSelection.xlsx',V,'OutputForV' ,'B2:C353') 
xlswrite( 'SupplierSelection.xlsx',N,'TOPSIS OUTPUT Variables' ,'B3:C354') 
xlswrite( 'SupplierSelection.xlsx',V,'TOPSIS OUTPUT Variables' ,'G3:H354')
xlswrite( 'SupplierSelection.xlsx',Apositive,'TOPSIS OUTPUT Variables' ,'B364:MO364')
xlswrite( 'SupplierSelection.xlsx',Anegetive,'TOPSIS OUTPUT Variables' ,'B365:MO365')
xlswrite( 'SupplierSelection.xlsx',dPositive,'TOPSIS OUTPUT Variables' ,'B361:MO361')
xlswrite( 'SupplierSelection.xlsx',dNegetive,'TOPSIS OUTPUT Variables' ,'B362:MO362')
xlswrite( 'SupplierSelection.xlsx',cc,'TOPSIS OUTPUT Variables' ,'B363:MO363')
xlswrite( 'SupplierSelection.xlsx',sumDmmMatrix,'Entropy Output' ,'B2:C353')
xlswrite( 'SupplierSelection.xlsx',pij,'Entropy Output' ,'G2:H353')
xlswrite( 'SupplierSelection.xlsx',wt,'Entropy Output' ,'L2:M353')
end

