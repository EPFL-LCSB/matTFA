function [model, drains, MinMets, cDPs] = analysisIMM_save(model, ...
    flagUpt, minObj, maxObj, NumAlt, drainsForiMM, rxnNoThermo, ...
    ReactionDB, metabData, save_path)

% Save original version of the model 
model_original = model;

% Set default values for optional inputs
if nargin < 2, flagUpt = true; end
if nargin < 8, ReactionDB = load('DB_AlbertyUpdate_keng_NEW.mat').DB_AlbertyUpdate; end
if nargin < 9, metabData = []; end

% Check for valid minObj and maxObj
if minObj > maxObj, error('minObj is greater than maxObj'); end

% Preprocess the model
[model, flagChange, drains, ~] = preprocessModel(model, drainsForiMM, flagUpt);
if flagChange || (numel(rxnNoThermo) == numel(model.rxns)), model = convToTFA(model, ReactionDB, rxnNoThermo); end
if ~isempty(metabData), model = loadConstraints(model, metabData); end

% Set the objective and bounds
model.var_lb(model.f == 1) = minObj;
model.var_ub(model.f == 1) = maxObj;
model.objtype = -1; % maximize

% Check for feasibility
sol = solveTFAmodelCplex_selections(model);
if sol.val < minObj
    error("growth lower bound constraint not imposed well")
end

model.f = zeros(length(model.varNames), 1);

% Add integer variables and constraints
[model, ~, ~] = addMILPVarsAndConstraints(model, drains);
   
% Get alternative solutions
MinMets = cell(NumAlt, 1);
[cDPs, model, indUSE] = findDPMinMets_informationalCallback(model, model_original, NumAlt, [],30*60, drains,minObj,save_path); 
if isempty(cDPs), disp('no solution found'); return; end

save('./workspaces/iMM_over.mat')


% Retrieve metabolite names
MinMets = retrieveMetaboliteNames(model, cDPs, indUSE);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions

function [model, flagChange, drains, drainMets] = preprocessModel(model, drainsForiMM, flagUpt)
    [model, flagChange, drains, drainMets] = putDrainsForward(model);
    if ~isempty(drainsForiMM)
        if all(ismember(drainsForiMM, model.rxns)) || all(ismember(drainsForiMM, drains))
            drains = drainsForiMM;
            drainMets = printRxnFormula(model, drains, 0, 0, 1);
        else
            warning('Not all drainsForiMM were identified as drains or rxns');
        end
    end
    if flagUpt
        prefix = 'R_';
    else
        prefix = 'F_';
    end
    n = numel(drains);

    for i = 1:n
        drains{i} = strcat(prefix,drains{i});
    end

    drains = [drains, drainMets];
end

function [model, rowDrain, num_vars] = addMILPVarsAndConstraints(model, drains)
    [~, num_vars] = size(model.A);
    [~, rowDrain] = ismember(drains(:, 1), model.varNames);
    model = addIntegerVariables(model, drains, num_vars);
    model = addConstraints(model, rowDrain, num_vars);
end

function model = addIntegerVariables(model, drains, num_vars)
    intTag = {'BFUSE'};
    model.varNames(num_vars + 1 : num_vars + size(drains, 1)) = strcat(intTag, '_', drains(:, 1));
    model.vartypes(num_vars + 1 : num_vars + size(drains, 1)) = {'B'};
    model.var_lb(num_vars + 1 : num_vars + size(drains, 1)) = 0;
    model.var_ub(num_vars + 1 : num_vars + size(drains, 1)) = 1;
    model.f(num_vars + 1 : num_vars + size(drains, 1)) = 1;
end

function model = addConstraints(model, rowDrain, num_vars)
    [num_constr,~ ] = size(model.A);
    for i = 1 : length(rowDrain)
        model.constraintNames{num_constr + i, 1} = strcat('BFMER_', num2str(i));
        model.constraintType{num_constr + i, 1} = '<';
        model.A(num_constr + i, rowDrain(i)) = 1;
        model.A(num_constr + i, num_vars + i) = 100;
        model.rhs(num_constr + i, 1) = 100;        
    end
end

function MinMets = retrieveMetaboliteNames(model, cDPs, indUSE)
    MinMets = cell(size(cDPs,2), 1);
    for i = 1 : size(cDPs, 2)
        rxnsBFUSE = model.varNames(indUSE(cDPs(indUSE, i) < 0.1));
        if isempty(rxnsBFUSE)
            error('no solution was found')
        else
            MinMets{i, 1} = strrep(rxnsBFUSE, 'BFUSE_R_EXC_BOTH_', '');
        end
    end
end
