function [gurobiModel, gurobiParams] = changeToGurobi_WithOptions_selections(tModel, options)


% Source for Gurobi MATLAB model and solution documentation: https://www.gurobi.com/documentation/5.0/refman/node650.html#sec:MATLAB
% Source for Gurobi MATLAB parameters documentation: https://www.gurobi.com/documentation/current/refman/parameters.html

% Set default values for options and merge with provided options
mergedOptions = mergeOptionsWithDefaults(options);

% Convert the model to GUROBI format using the merged options
gurobiModel = setupModelProperties(tModel);

% Set up optimization parameters for GUROBI
gurobiParams = setupOptimizationParameters(mergedOptions);

end

%%

function mergedOptions = mergeOptionsWithDefaults(options)
    defaultOptions = struct(...
        'manualScalingFactor', [], ...
        'mipTolInt', 1e-6, ...
        'mipTolRel', 1e-4, ...
        'emphPar', 1, ...
        'feasTol', 1e-6, ...
        'scalPar', 0, ...
        'TimeInSec', 3600, ...
        'barrierDisplay',2, ...
        'mipDisplay', 5, ...
        'CPXPARAMdisp', 2, ...
        'timePolishing', 1e75, ...
        'preSolve', 1, ...
        'varSelec', 2, ...
        'hFreq', 10, ...
        'treeMem', 1e75 ...
    );
    
    mergedOptions = defaultOptions;
    if isempty(options)
        return;
    end
    
    optionNames = fieldnames(options);
    for i = 1:numel(optionNames)
        if ~isempty(options.(optionNames{i}))
            mergedOptions.(optionNames{i}) = options.(optionNames{i});
        end
    end

end


%%
function gurobiModel = setupModelProperties(model)

% Convert the constraint types into the right format
contypes = '';
num_constr = length(model.constraintType);
for i=1:num_constr
    contypes = strcat(contypes, model.constraintType{i,1});
end

% convert variable types into the right format
vtypes = '';
num_vars = length(model.vartypes);
for i=1:num_vars
    vtypes = strcat(vtypes, model.vartypes{i,1});
end

% Set up model properties
gurobiModel.A = model.A;
gurobiModel.obj = model.f;
gurobiModel.lb = model.var_lb;
gurobiModel.ub = model.var_ub;
gurobiModel.rhs = model.rhs;
gurobiModel.sense = contypes;
gurobiModel.vtype = vtypes;
gurobiModel.varnames = model.varNames;

% Set up model objective
if model.objtype == -1
    gurobiModel.modelsense = 'max';
elseif model.objtype == 1
    gurobiModel.modelsense = 'min';
else
    error(['No objective type specified ' ...
        '(model.objtype should be in {-1,1})']);
end


end

%%
function params = setupOptimizationParameters(options)

% Set default values for optimization parameters
params.LogToConsole = 1;
params.IntFeasTol = options.mipTolInt;
params.MIPGap = options.mipTolRel;
params.DisplayInterval = options.mipDisplay;
params.FeasibilityTol = options.feasTol;
params.TimeLimit = options.TimeInSec;

% params.ScaleFlag = options.manualScalingFactor;

% Uncomment the following lines if these options are needed

% params.LogFile = "";
% params.preprocessing_presolve = options.preSolve;
% params.mip_limits_treememory = options.treeMem;
% params.barrier_display = options.barrierDisplay;
% params.mip_strategy_nodeselect = options.nodeSelec;
% params.mip_strategy_variableselect = options.varSelec;
% params.mip_strategy_heuristicfreq = options.heurFreq;
%params.mip_polishafter_time = options.timePolishing;
%params.emphasis_numerical = options.emphPar;
%params.emphasis_mip = options.emphMIPswitch;

end

