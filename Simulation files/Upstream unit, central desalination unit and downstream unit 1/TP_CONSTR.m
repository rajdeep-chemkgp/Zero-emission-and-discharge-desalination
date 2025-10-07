%*************************************************************************
% Test Problem : 'CONSTR'
% Description:
%   (1)constrained
% Reference : [1] Deb K, Pratap A, Agarwal S, et al. A fast and elitist 
%   multiobjective genetic algorithm NSGA-II[J]. Evolutionary Computation.2002, 6(2): 182-197.
%*************************************************************************
clear all
clc
options = nsgaopt();                     % create default options structure
options.popsize = 100;                   % populaion size
options.maxGen  = 50;                   % max generation
options.numObj = 4;                      % number of objectives
options.numVar = 6;                      % number of design variables
options.vartype=[1 1 1 2 1 1 ];         % TYPE OF VARIABLES 2 FOR INTEGER AND 1 FOR REAL
options.numCons = 3;                     % number of constraints
options.lb = [65 52  45  4     9     1500    ];               % lower bound of x
options.ub = [80 55  50  12    20    3000    ];               % upper bound of x 95%AVG 
% options.lb = [65 45.4  41.5  4  8.89  2100 171];               % lower bound of x
% options.ub = [85 45.4  41.5  4  8.89  2100 171]; 
options.objfun = @TP_CONSTR_objfun;     % objective function handle
more off
%options.plotInterval = 50;               % interval between two calls of "plotnsga". 
%NT,Q,TAC,PTIME,precovery,RefluxRatio,totalbatch,Annualprod
result = nsga2(options);                % begin the optimization!
