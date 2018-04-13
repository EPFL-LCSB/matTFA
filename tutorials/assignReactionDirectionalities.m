function model = assignReactionDirectionalities(model, minFluxValue)
% Assigns a directionality upon arbitrary sequential fixation of the
% bi-directional reactions (BDR)
%
% INPUTS
%   model : FBA and TFA model structure that has bi-directional reactions
%   in both the FBA and the TFA part.
%   minFluxValue : Define a minimum flux value. This value will be used to
%   set the directionaity of the reactions instead of zero (to avoid zero
%   solution fluxes that do not indicate directionality)
% OUTPUT
%   model : FBA and TFA model structure that has no bi-directional
%   reactions in both the FBA and the TFA part.

% Make sure that the model structure contains the FBA and the TFA part
if ~isfield(model, 'S')
    error('This is not an FBA model')
end
if ~isfield(model, 'A')
    error('This is not a TFA model')
end
if isempty('NFids')
    error('This model has no net-flux variables')
end
if ~exist('minFluxValue','var') || isempty(minFluxValue)
    % Define a minimum flux value
    minFluxValue = 1e-6;
end
if isempty('NFids')
    error('This model has no net-flux variables')
end

% What is the tolerance of the solver? 1e-9?
SolTol = 1e-9;

% Get indices of BDRs
NFids = getAllVar(model, {'NF'});
tmm = runTMinMax(model, model.varNames(NFids));
id_BD = find(tmm(:,1)<-SolTol & tmm(:,2)>SolTol);
% Get directionalities of rest reactions (accounting for precision errors)
id_F = tmm(:,1)>=-SolTol & tmm(:,2)>SolTol;
id_R = tmm(:,1)<-SolTol & tmm(:,2)<=SolTol;

% Check if the number of forward, reverse and BDR sums up to the total
% number of the model
if ~ (sum(id_F) + sum(id_R) + size(id_BD, 1) == size(model.rxns, 1))
    error('Some of the directionality assignment of the reactions is incinsistent!')
end

% Set the directionality of the unidirectional reactions to avoid zeros in
% the solution
model.var_lb(NFids(id_F)) = minFluxValue;
model.var_ub(NFids(id_R)) = -minFluxValue;

solTFA = solveTFAmodelCplex(model);
if isempty(solTFA.val) || ~isnumeric(solTFA.val)
    error('Imposing the minimum flux value affects feasibility')
end

% By sequential attempt of each direction we define the directionalities
% of each of the BDRs
% Initialize a vector with the directinoalities of the BDR
directBD = zeros(size(id_BD));
for i = 1:length(id_BD)
    % Initialize the variable for attempting to impose the minimum flux
    failedToImposeMinFlux = true;
    % Assume forward directionality at first for the i-th bi-directinoal
    Directionality_i = 1;
    % Get the NF ranges of the ith net flux
    id_NFi_inVarNames = find_cell(['NF_',model.rxns{id_BD(i)}], model.varNames);
    range_NFi = [model.var_lb(id_NFi_inVarNames) model.var_ub(id_NFi_inVarNames)];
    while failedToImposeMinFlux
        if Directionality_i == 1;
            model.var_lb(id_NFi_inVarNames) = minFluxValue;
            model.var_ub(id_NFi_inVarNames) = range_NFi(2);
            directBD(i) = 1;
        elseif Directionality_i == -1;
            model.var_lb(id_NFi_inVarNames) = range_NFi(1);
            model.var_ub(id_NFi_inVarNames) = -minFluxValue;
            directBD(i) = -1;
        end
        solTFA_i = solveTFAmodelCplex(model);
        if ~isempty(solTFA_i.val) && isnumeric(solTFA_i.val)
            % If we got a solution we can exit the while loop
            failedToImposeMinFlux = false;
        else
            % If we didn't find a solution we stay in the while loop and
            % revert the directionality of the reaction
            failedToImposeMinFlux = true;
            Directionality_i = -1*Directionality_i;
        end
    end
end

% Get the signs of the net-fluxes of a solution of the TFA model
solTFA = solveTFAmodelCplex(model);
signs_NF_TFA = sign(solTFA.x(NFids));
% And check that that there are no zeros in this vector!
if ~isempty(find(signs_NF_TFA==0))
    error('There exist still zeros in the solution!')
end

% Impose minimum flux bounds also to the FBA.
% ATTENTION: Do not modify those bounds that already have a significant
% flux towards any of the two directions
idRxns_Significant_F_Flux = model.lb >  minFluxValue;
idRxns_Significant_R_Flux = model.ub < -minFluxValue;
model.lb(signs_NF_TFA==1  & ~idRxns_Significant_F_Flux) =  minFluxValue;
model.ub(signs_NF_TFA==-1 & ~idRxns_Significant_R_Flux) = -minFluxValue;

% Ensure that the signs of net-fluxes of both models are the same
solFBA = solveFBAmodelCplex(model);
signs_FBA = sign(solFBA.x);
if ~isequal(signs_NF_TFA, signs_FBA)
    error('The signs of FBA and TFA solutions are not the same!')
end

end