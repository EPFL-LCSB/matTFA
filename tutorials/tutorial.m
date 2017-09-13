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

soltFA = solveTFAmodel(this_tmodel);

%% Perform FVA 

fva = runMinMax(this_tmodel);

%% Perform TVA

% Get the variables representing the net fluxes
NF_ix = getAllVar(this_tmodel,{'NF'});
tva = runTMinMax(this_tmodel, this_tmodel.varNames(NF_ix));

%% Plot the difference
%%
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
tmp = prepped_m;
tmp = convToTFA(prepped_m,ReactionDB,[],'DGo', min_obj, [], 1);
gamma_model = addNetFluxVariables(tmp);

lngamma = soltFA.x(getAllVar(gamma_model,{'LnGamma'}));


%% Sampling
basalFluxValue = 1e-6;
sampled_model = sampleDP(gamma_model, soltFA, 1, [], tva, basalFluxValue)
