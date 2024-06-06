function [gurobiModel, gurobiParams] = changeToGurobi_FBA_WithOptions_selections(model, options)

% Set default values for options and merge with provided options
mergedOptions = mergeOptionsWithDefaults(options);

% Convert the model to GUROBI format using the merged options
gurobiModel = setupModelProperties(model);

% Set up optimization parameters for GUROBI
gurobiParams = setupOptimizationParameters(mergedOptions);

end

%%

function mergedOptions = mergeOptionsWithDefaults(options)
    defaultOptions = struct(...
        'scalPar', -1, ...
        'feasTol', 1e-9, ...
        'emphPar', 1, ...
        'TimeInSec', 3600);
    
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



%% Function to set up model properties in the CPLEX model
function gurobiModel = setupModelProperties(model)

    % Set the model name
    gurobiModel.modelname = char(model.description);
    
    % Set the optimization direction (maximize or minimize)
    if ~isfield(model, 'osense')
        gurobiModel.modelsense = 'max';
    else
        if model.osense == -1
            gurobiModel.modelsense = 'max';
        elseif model.osense == 1
            gurobiModel.modelsense = 'min';
        else
            error("Wrong value for the osense field of the model structure");
        end
    end

    % Set the constraint matrix and names of the rows/columns
    gurobiModel.A = sparse(model.S);
    gurobiModel.constrnames = model.mets;
    gurobiModel.varnames = model.rxns;

    % Set the objective function coefficients and name
    gurobiModel.obj = model.c;

    % Set bounds for the variables
    gurobiModel.lb = model.lb;
    gurobiModel.ub = model.ub;

    % Lefthand and righthand side values for the constraints in the constraint matrix A
    gurobiModel.rhs = model.b;

end

%%
function params = setupOptimizationParameters(options)

% Set default values for optimization parameters
params.LogToConsole = 1;
params.DisplayInterval = 5;
params.FeasibilityTol = options.feasTol;
params.TimeLimit = options.TimeInSec;

end
