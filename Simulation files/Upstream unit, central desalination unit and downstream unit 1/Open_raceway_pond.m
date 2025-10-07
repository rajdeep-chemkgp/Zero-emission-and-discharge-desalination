function[Pond_area,NS_growth,Total_cost_ORP,SP_lipids,energy_lipid_extraction_kW,cost_lipid_extraction_per_year,NS_growth_m3_year]=Open_raceway_pond(CO2_emission,THY)

% Converting CO2 emission to kg/day
Co2_emission_day=CO2_emission*86400 % one day=86400 secs/day (kg/day)
% For 1 m2 of algae culture, the CO2 consumption in m3/day
algaegrowthrate=20*(10^(0))*0.9; % algae growth rate of 20 gm/m2day with 90% harvest efficiency
CO2_consumption_m3_day=0.6565*algaegrowthrate+5.0784; % Optimal integration of a self sustained algae based facility with solar and/or wind energy 
CO2_consumption_kg_day=CO2_consumption_m3_day*1.98 % 1.98 Density of CO2
% The above is the amount of CO2 consumed in 1 day per m2
Pond_area_m2=Co2_emission_day/CO2_consumption_kg_day; 
Pond_area_km2=Pond_area_m2/(1000*1000); %m2

% growth rate of Nannochloropsis salina is 20 g/m2.day
% 0.9 is the plant load factor and 365 is the number of days
NS_growth=(20*Pond_area_m2*0.9*365)/1000;  % (kg) per year
desnity=900; %kg/m3; Density of algal lipids: 0.9 g/ml: Composition of Algal Oil and Its Potential as Biofuel (https://onlinelibrary.wiley.com/doi/10.1155/2012/285185) 
NS_growth_m3_year=NS_growth/desnity
volume=(NS_growth/desnity)*264.172; %gallons conversion factor
Lipids=NS_growth*0.3
SP_lipids=(0.3*(NS_growth/desnity)*264.172)*8.52
energy_lipid_extraction=(NS_growth*0.3)*2.16*1000; %kJ
energy_lipid_extraction_kW=energy_lipid_extraction/(THY*3600*25);
energy_lipid_extraction_kWh_year=energy_lipid_extraction_kW*THY;
cost_lipid_extraction_per_year=energy_lipid_extraction_kWh_year*0.049
NS_growth_tonnes=NS_growth*0.001; %tonne/y
Total_cost_ORP=NS_growth_tonnes*673.65; % $/year Techno-economic assessment of open microalgae production systems (DOI: https://doi.org/10.1016/j.algal.2017.01.005)
Pond_area=Pond_area_m2*(10^(-6)); 
% Cost_of_land=7.4*Pond_area; % $
% Cost_of_pond=34*Pond_area;

end