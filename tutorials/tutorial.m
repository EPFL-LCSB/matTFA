clear
restoredefaultpath

%% Load solver

cplex_loaded = load_cplex;

if ~cplex_loaded
    error('CPLEX not found')
end

%% Load the COBRA Toolbox
addpath(genpath(fullfile('..','ext')))

%%  Load dependencies
addpath(genpath(fullfile('..','matTFA')))
addpath(genpath(fullfile('..','thermoDatabases')))
addpath(genpath(fullfile('..','models')))
addpath(genpath(fullfile('..','plotting')))

% switch to cplex solver
% changeCobraSolver('ibm_cplex')
changeCobraSolver('cplex_direct')

tmp = load(fullfile('..','models','small_ecoli.mat'));
mymodel = tmp.model_red;
clear tmp

tmp = load('thermo_data.mat');
ReactionDB = tmp.DB_AlbertyUpdate;
clear tmp

%Fix big stoichiometries
[~,is_big_stoich] = find(sum(mymodel.S,1) > 100);
mymodel.S(:,is_big_stoich) = 1/10 * mymodel.S(:,is_big_stoich);


%% Test simple fba
solFBA = optimizeCbModel(mymodel);

%% Prepare for TFA

%need field for description
prepped_m = prepModelforTFA(mymodel,ReactionDB,mymodel.CompartmentData);

%% Convert to TFA
min_obj = solFBA.f*0.95;
tmp = convToTFA(prepped_m,ReactionDB,[],'DGo', min_obj);

% add net flux variables, which are equal to forwards flux - backwards flux
% NF_rxn = F_rxn - B_rxn
this_tmodel = addNetFluxVariables(tmp);

%% Solve tFA
save tutorial_01.mat
soltFA = solveTFBAmodelCplex(this_tmodel);

%% Perform FVA 

fva = runMinMax(this_tmodel);

%% Perform TVA

% Get the variables representing the net fluxes
NF_ix = getAllVar(this_tmodel,{'NF'});
tva = runTMinMax(this_tmodel, this_tmodel.varNames(NF_ix));
save tutorial_02.mat

%% Plot the difference
f = @(x) abs(x(:,2) - x(:,1));
m = @(x) abs(x(:,2) + x(:,1))/2;
n = @(x) (x(:,2) .* x(:,1)) < 0;
s = @(x,y) (f(x)-f(y))./f(x) .* (f(x) > 0.01)  + n(x);

loss = s(fva,tva);

[~,rix] = sort(loss)
is_bd = (n(fva))

to_plot = intersect(find(is_bd),rix)
xlabels = this_tmodel.rxns(to_plot)
figure
hold on
p1 = p_patch_va(fva(to_plot,:),0.3,[20 43 140]./255)
p2 = p_patch_va(tva(to_plot,:),0.25,[222 125 0]./255)
plot(xlim,[0 0],'k:')
set(gca,'XTick',1:sum(is_bd),'XTickLabel',xlabels,'XTickLabelRotation',-45,'TickLabelInterpreter','none')
legend([p1 p2] , 'Initial constraints', 'Thermodynamics','Location','southwest')
ylabel ('flux [mmol.gDw^{-1}]')
title('FVA difference after constraints on bidirectional reactions')

%% Get Thermodynamic displacements

% Add thermo_disp as variables
basalFlux = 1e-6;
gamma_model = convToTFA(prepped_m, ReactionDB, [], 'DGo', min_obj, [], 1, true);
gamma_model = addNetFluxVariables(gamma_model);
gamma_model = addMinFluxRequirements(gamma_model, basalFlux, [], [], false);

solTFA = solveTFAmodelCplex(gamma_model);
NFids = getAllVar(gamma_model,{'NF'});
LnGammaids = getAllVar(gamma_model,{'LnGamma'});
lngammaValues = solTFA.x(getAllVar(gamma_model,{'LnGammaids'}));
save tutorial_03.mat

%% Sampling
% basalFluxValue = 1e-6;
% sampled_model = sampleDPforMCA(gamma_model, soltFA, 1, [], tva, basalFluxValue)

% Set the directionalities of the FBA model According to the solution
% obtained by the TFA model.
gamma_model = addTFASolBoundsToFBAmodel(gamma_model, solTFA.x, NFids);



% Settings of the sampler (COBRA)
soloptions.nWarmupPoints = 500;
soloptions.nFiles = 1;
soloptions.nPointsPerFile = 5000;
soloptions.nStepsPerPoint = 100;
soloptions.nPointsReturned = 5000;
soloptions.NvZeroTolerance = 1e-8;
soloptions.nFilesSkipped = 0;
soloptions.removeLoopsFlag = false;
soloptions.removeLoopSamplesFlag = true;

% Sample Fluxes




