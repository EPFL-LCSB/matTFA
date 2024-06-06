function minmax = runMinMax_selections(model, rxnNames, verbose, scalPar, feasTol, emphPar)

if ~exist('scalPar','var') || isempty(scalPar)
    scalPar = [];
end
if ~exist('feasTol','var') || isempty(feasTol)
    feasTol = [];
end
if ~exist('emphPar','var') || isempty(emphPar)
    emphPar = [];
end
if ~exist('rxnNames','var') || isempty(rxnNames)
    rxnNames = model.rxns;
end
if ~exist('verbose','var') || isempty(verbose)
    verbose = true;
end

minmax = zeros(length(rxnNames),2);
rxn_id = find_cell(rxnNames, model.rxns);

for i = 1:length(rxn_id)
    
    if verbose && ~isfield(model,'CS_varNames')
        fprintf('minmax for %s\t',model.rxns{rxn_id(i)});
        model = changeObjective(model,model.rxns{rxn_id(i)});
    elseif verbose && isfield(model,'CS_varNames')
        fprintf('minmax for %s\t',model.CS_varNames{rxn_id(i)});
        model.c = zeros(size(model.S,2),1);
        model.c(rxn_id(i)) = 1;
    end
    
    % Minimization
    model.osense = 1;
    sol = solveFBAmodel_selections(model, 'scalPar', scalPar, 'feasTol', feasTol, 'emphPar', emphPar);
    minmax(i,1) = sol.f;
    
    % Maximization
    model.osense = -1;
    sol = solveFBAmodel_selections(model, 'scalPar', scalPar, 'feasTol', feasTol, 'emphPar', emphPar);
    minmax(i,2) = sol.f;
    
    if verbose
        fprintf('min: %d\t max: %d\n', minmax(i,1), minmax(i,2));
    end
    
end