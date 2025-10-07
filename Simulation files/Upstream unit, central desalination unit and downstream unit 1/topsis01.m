function [ cc ] = topsis (decisionMakingMatrix,lambdaWeight,criteriaSign )
 %-Technique for Order of Preference by Similarity to Ideal Solution
 % -Author:Omid Ameri Sianaki
 %-This function implements TOPSIS method with Information entropy
 % weighting Method
% - criteriaSign is a vector specifying whether an criterion has to be maximized
% or minized . +1 is for positive criterion and -1 for negetive criterion
% -Before executing the function you have to define "decisionMakingMatrix"
% variable based on size of decision making matrix that you have.
% -It is very importand to observe the compatibility of dimension in decision making and
% weighing matrix.
% - I offer you study Page 79 of my PhD thesis in below address:
% http://espace.library.curtin.edu.au/R?func=dbin-jump-full&local_base=gen01-era02&object_id=240088
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

end

