function cplex = changeToCPLEX_FBA_WithOptions_selections(model, options)

% Set default values for options and merge with provided options
mergedOptions = mergeOptionsWithDefaults(options);

% Convert the model to CPLEX format using the merged options
cplex = setupModelProperties(model);

% Set up optimization parameters for CPLEX
cplex = setupOptimizationParameters(cplex, mergedOptions);

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
function cplex = setupModelProperties(model)
% SETUPMODELPROPERTIES Sets up the properties of the CPLEX model.
% This function assigns the name, constraint matrix, objective function,
% variable bounds, and constraint senses to the CPLEX model based on the input model structure.
% Source: https://www.ibm.com/docs/en/icos/12.10.0?topic=list-cplex#a051e4d035e053a4636efc58c1bde9b3e

    % Create CPLEX model object with the model description
    cplex = Cplex(model.description);

    % Set the model name
    cplex.Model.name = string(model.description);

    % Set the optimization direction (maximize or minimize)
    if ~isfield(model, 'osense')
        cplex.Model.sense = 'max';
    else
        if model.osense == -1
            cplex.Model.sense = 'max';
        elseif model.osense == 1
            cplex.Model.sense = 'min';
        else
            error("Wrong value for the osense field of the model structure");
        end
    end

    % Set the constraint matrix and names of the rows/columns
    cplex.Model.A = model.S;
    cplex.Model.rowname = char(model.mets);
    cplex.Model.colname = char(model.rxns);

    % Set the objective function coefficients and name
    cplex.Model.obj = model.c;
    cplex.Model.objname = strjoin(string(model.rxns(find(model.c))), " + ");

    % Set bounds for the variables
    cplex.Model.lb = model.lb;
    cplex.Model.ub = model.ub;

    % Lefthand and righthand side values for the constraints in the constraint matrix A
    if ~isfield(model, 'csense')
        % If csense is not declared in the model, assume that all constraints are equalities.
        csense(1:length(model.mets), 1) = 'E';
    else
        % If csense is in the model, ensure its length matches the number of metabolites
        if length(model.csense) ~= length(model.mets)
            warning('Length of csense is invalid! Defaulting to equality constraints.');
            csense(1:length(model.mets), 1) = 'E';
        else
            csense = columnVector(model.csense);
        end
    end

    % Fill in the RHS vector if not provided
    if ~isfield(model, 'b')
        b = zeros(size(model.S, 1), 1);
    else
        b = model.b;
    end
    b = full(b);

    if ~isempty(csense)
        % Set up constant vectors for CPLEX
        b_L(csense == 'E', 1) = b(csense == 'E');
        b_U(csense == 'E', 1) = b(csense == 'E');
        b_L(csense == 'G', 1) = b(csense == 'G');
        b_U(csense == 'G', 1) = Inf;
        b_L(csense == 'L', 1) = -Inf;
        b_U(csense == 'L', 1) = b(csense == 'L');
    else
        b_L = b;
        b_U = b;
    end

    % Set the lefthand and righthand side values in the CPLEX model
    cplex.Model.lhs = b_L;
    cplex.Model.rhs = b_U;
end


%%
function cplex = setupOptimizationParameters(cplex, options)
% Source: https://www.ibm.com/docs/en/icos/12.10.0?topic=list-cplex#a56d92885c62355169101500c45638e29

% Set default values for optimization parameters
params = {
    'output.clonelog', 0;
    'read.scale', options.scalPar;
    'simplex.tolerances.feasibility', options.feasTol;
    'emphasis.numerical', options.emphPar;
    'timelimit', options.TimeInSec;
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
    options.CPXPARAMdisp = 1;
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
