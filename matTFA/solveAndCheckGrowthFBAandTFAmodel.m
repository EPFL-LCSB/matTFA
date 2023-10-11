function [solFBA, solTFA] = solveAndCheckGrowthFBAandTFAmodel(model,biomass_special,flagTFA)

    if ~exist('flagTFA','var') || isempty(flagTFA)
        flagTFA = 1;
    end

    % Set objective to biomass
    if ~exist("biomass_special",'var') || isempty(biomass_special)
        model = set_objective(model, {'biomass'}, 1, "max", flagTFA);
    else
        model = set_objective(model, biomass_special, 1, "max", flagTFA);
    end

    % Solve FBA model
    solFBA = solveFBAmodelCplex(model);

    % Solve TFA model
    if flagTFA
        solTFA = solveTFAmodelCplex_selections(model,'TimeInSec',2*60,'emphPar',0);
    end

    % Check if the solution is valid
    if isnan(solFBA.f) || isempty(solFBA.f) || solFBA.f < 1e-9
        error('Model not growing');
    end

    % Check if the thermo solution is valid
    if flagTFA
        if isnan(solTFA.val) || isempty(solTFA.val) || solTFA.val < 1e-9
            error('thermo Model not growing');
        end
    end


    % Print solFBA.val and solTFA.val results on separate lines
    fprintf('%s\n', repmat('#', 1, 20));
    fprintf('solFBA.val: %f\n', solFBA.f);
    if flagTFA
        fprintf('solTFA.val: %f\n', solTFA.val);
    end

end
