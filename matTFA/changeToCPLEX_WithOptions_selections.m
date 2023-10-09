function cplex = changeToCPLEX_WithOptions_selections(tModel, options)

% Set default values for options and merge with provided options
mergedOptions = mergeOptionsWithDefaults(options);

% Convert the model to CPLEX format using the merged options
cplex = convertModelToCPLEX(tModel, mergedOptions);

% Set up optimization parameters for CPLEX
cplex = setupOptimizationParameters(cplex, mergedOptions);

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
function cplex = convertModelToCPLEX(model, options)
% Manual scaling of problem (if specified)
if ~isempty(options.manualScalingFactor)
    model.A = options.manualScalingFactor * model.A;
    model.rhs = options.manualScalingFactor * model.rhs;
end

% Initialize and set up the CPLEX model
cplex = initializeCplexModel(model);
setupModelProperties(cplex, model);

end

%%
function cplex = initializeCplexModel(model)
vtypes = cell2mat(model.vartypes');
cplex = Cplex(model.description);
cplex.Model.ctype = vtypes;
if (model.objtype == -1)
    cplex.Model.sense = 'maximize';
else
    cplex.Model.sense = 'minimize';
end
end

%%
function setupModelProperties(cplex, model)
% Set up model properties
cplex.Model.A = model.A;
cplex.Model.obj = model.f;
cplex.Model.lb = model.var_lb;
cplex.Model.ub = model.var_ub;

% Set up constraint bounds
lhs = [];
rhs = [];

for i = 1:size(model.A, 1)
    switch model.constraintType{i, 1}
        case '='
            lhs(i, 1) = model.rhs(i);
            rhs(i, 1) = model.rhs(i);
        case '>'
            lhs(i, 1) = model.rhs(i);
            rhs(i, 1) = inf;
        case '<'
            lhs(i, 1) = -Inf;
            rhs(i, 1) = model.rhs(i);
        otherwise
            error('Constraint type not recognized.');
    end
end

cplex.Model.lhs = lhs;
cplex.Model.rhs = rhs;

% Set up optional fields
if isfield(model, 'Q'), cplex.Model.Q = model.Q; end
if isfield(model, 'varNames'), cplex.Model.colname = char(model.varNames); end
if isfield(model, 'constraintNames'), cplex.Model.rowname = char(model.constraintNames); end

end

%%
function cplex = setupOptimizationParameters(cplex, options)
% Set default values for optimization parameters
params = {
    'output.clonelog', 0;
    'mip.tolerances.integrality', options.mipTolInt;
    'mip.tolerances.mipgap', options.mipTolRel;
    'mip.display', options.mipDisplay;
    'mip.polishafter.time', options.timePolishing;
    'emphasis.numerical', options.emphPar;
    'emphasis.mip', options.emphMIPswitch;
    'simplex.tolerances.feasibility', options.feasTol;
    'read.scale', options.scalPar;
    'timelimit', options.TimeInSec;
    'preprocessing.presolve', options.preSolve;
    'mip.strategy.variableselect', options.varSelec;
    'mip.limits.treememory', options.treeMem;
    'barrier.display' , options.barrierDisplay;
    'mip.strategy.nodeselect', options.nodeSelec;
    'mip.strategy.variableselect',options.varSelec;
    'mip.strategy.heuristicfreq',options.heurFreq;
};


% Set the parameter values in CPLEX
for i = 1:size(params, 1)
    paramNameParts = split(params{i, 1}, '.');
    switch numel(paramNameParts)
        case 3
            cplex.Param.(paramNameParts{1}).(paramNameParts{2}).(paramNameParts{3}).Cur = params{i, 2};
        case 2
            cplex.Param.(paramNameParts{1}).(paramNameParts{2}).Cur = params{i, 2};
        case 1
            cplex.Param.(paramNameParts{1}).Cur = params{i, 2};
        otherwise
            warning('Invalid parameter name: %s. Skipping this parameter.', params{i, 1});
    end
end


% Set display function
if ~exist('CPXPARAMdisp', 'var') || isempty(options.CPXPARAMdisp)
    CPXPARAMdisp = 1;
else
    if ~ismember(options.CPXPARAMdisp, [0 1])
        error('Parameter value out of range!')
    end
end
if options.CPXPARAMdisp == 0
    CPXPARAMdisp = [];
else
    CPXPARAMdisp = @disp;
end
cplex.DisplayFunc = CPXPARAMdisp;

end
