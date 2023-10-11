function model = addMetabolite(model, metId, metName, metFormula, metCharge)
% addMetabolite Add a metabolite to the model
%
% USAGE:
%   model = addMetabolite(model, metId, metName, metFormula, metCharge)
%
% INPUTS:
%   model       - The model to which the metabolite will be added.
%   metId       - String representing the ID of the new metabolite.
%   metName     - String representing the name of the new metabolite.
%   metFormula  - String representing the chemical formula of the new metabolite.
%
% OPTIONAL INPUT:
%   metCharge   - Integer representing the charge of the new metabolite. 
%                 If not provided, it defaults to 0.
%
% OUTPUT:
%   model       - Model with the new metabolite added.
%
% NOTE:
%   This function will add a new metabolite to the model, expanding all
%   relevant fields. The stoichiometric matrix (model.S) will be expanded with 
%   a new row of zeros.

% Check if metCharge is provided
if nargin < 5
    metCharge = 0;
end

% Validate that the input IDs and Names do not already exist
if any(strcmp(model.mets, metId))
    error('Metabolite ID already exists in the model.');
end

if any(strcmp(model.metNames, metName))
    error('Metabolite name already exists in the model.');
end

% Add new metabolite ID
model.mets{end+1,1} = metId;

% Add new metabolite name
model.metNames{end+1,1} = metName;

% Add new metabolite formula
model.metFormulas{end+1,1} = metFormula;

% Add new metabolite charge
model.metCharge(end+1,1) = metCharge;

% Expand stoichiometric matrix with a new row of zeros
model.S(end+1, :) = 0;

end
