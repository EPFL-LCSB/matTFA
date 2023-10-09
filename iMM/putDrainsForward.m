function [model, flagChange, drains, drainMets] = putDrainsForward(model)
% Verifies that the model has drains defined as A => (and not => A),
% corrects if this is not the case, flags the change, and finds drained
% metabolites and reactions
%
% USAGE:
%
%       [model, flagChange, drains, drainMets] = putDrainsForward(model)
%
% INPUT:
%    model:           COBRA model structure
%
% OUTPUTS:
%    model:           Model with drains corrected from "=> A" to "A =>"
%    flagChange:      Defines if any exchange in the model was modified 
%                     ("0" means no change was made and "1" means at least 
%                     a change was made)
%    drains:          Exchange/Drain reactions
%    drainMets:       Metabolites exchanged
%
% .. Author:
% Meri? Ataman 2014
% 

% Preallocate drains and drainMets for efficiency
drains = cell(length(model.rxns),1);
drainMets = drains;
flagChange = 0;
counter = 0;

for i=1:length(model.rxns)
    met_indices = find(model.S(:, i));
    no_of_mets = length(met_indices);

    % Check for drain reactions
    if no_of_mets==1
        counter = counter + 1;
        drains{counter} = model.rxns{i};
        drainMets{counter} = model.mets{met_indices};

        % Check the stoichiometry and flip if needed
        if model.S(met_indices, i)==1
            flagChange = 1;
            model.S(:, i) = -model.S(:, i);
            [model.ub(i), model.lb(i)] = deal(-model.lb(i), -model.ub(i));
        end
    end
end

% Remove unused cells in drains and drainMets
drains(counter+1:end) = [];
drainMets(counter+1:end) = [];

end
