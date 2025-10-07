function [MWA,DENSA] = MWDENS(X)
% UNIT: DENSA=ld/ft3  & MWA=ldm/ldmol
% X = [0.25 0.3 0.25 0.2];
% Methanol -- TGA -- Biodiesel -- Glycerol %
MW = [32.04  831.62  280.18 92.094]; % molecular wt. (unit: - gm/mol) %
DENS = [0.792  0.892  0.860 1.261]; % Density (unit: - gm/cm3) %

for i = 1:4
    MWa(i) = MW(i)*X(i);
    DENSa(i) = DENS(i)*X(i);
end
MWA = sum(MWa);
DENSA = (sum(DENSa))*62.428;%(Unit: lb/ft3)
end