function [DPs, model, indUSE] = findDPMinMets_informationalCallback(model, model_original, NumAlt, intTag, time, drains, mu_lb, save_path, biomassRxn)
% Get alternative solutions for MILP (maximization)
%
% USAGE:
%       [DPs, model, objectives] = findDPMinMets(model, NumAlt, indUSE, time)
%
% INPUTS:
%    model:           Model with TFA structure and MILP formulation
%
% OPTIONAL INPUTS:
%    NumAlt:          Maximum number of alternatives to get (default = 1)
%    intTag:          Integer tag used in the MILP (default = 'BFUSE')
%    time:            Time in sec after which simulation will be stopped (default = empty / no time provided)
%
% OUTPUTS:
%    DPs:             Directionality profile matrix with alternatives in each column
%    model:           Model with integer cuts integrated to avoid repetition of same solution or supersolution
%    objectives:      Exchange/Drain reactions
%
% .. Author:
% Anush Chiappino 2017
% 

if nargin < 2 || isempty(NumAlt), NumAlt = 1; end
if nargin < 3 || isempty(intTag), intTag = 'BFUSE'; end
if nargin < 4, time = []; end

indUSE = getAllVar(model, {intTag});
[~,indDrainsR] = ismember(drains(:,1),model.varNames);
NumSols = 0;
sol = solveTFAmodel_selections(model);
stop_search = sol.val;
DPs = [];
objectives = [];
actUSEvecs = {};

if ~isempty(sol.x), model = addIntegerCutConstraint(model, indUSE, sol.val-0.5); end

while (NumSols < NumAlt) && ~isempty(sol.x) && (sol.val ~= 0)
    
    NumSols = NumSols + 1;
    objectives(NumSols, 1) = sol.val;
    DPs(:, NumSols) = sol.x;
    
    actUSEvec = indUSE(sol.x(indUSE) < 0.1);
    actUSEvecs = [actUSEvecs;{actUSEvec}];

    model = addIntegerCutConstraint(model, actUSEvec, 0.5);

    sol = solveTFAmodel_selections(model,'TimeInSec',time,'timePolishing',time*3/4,'stop_search_value',stop_search);
    
    if isempty(sol.x) || (sol.val == 0), break; end
    fprintf('Number of DPs:\t%d\n', NumSols);

    % Every 40 iMMs, make sure the model can grow on the latest one
    if rem(NumSols,40) == 0

        % Sanity check 1: make sure all the elements are same size or smaller size
        vectorSizes = cellfun(@numel, actUSEvecs);
        max_size = size(drains,1) - stop_search;
        isValid = all(vectorSizes <= max_size);
        
        if ~isValid
            error('Error: At least one vector in the cell array has a size larger than %d.\n',max_size);
        else
            fprintf('All vectors in the cell array have a size equal to or less than %d.\n',max_size);
        end

        % Sanity check 2 to ensure actUSEvec is different from previous iterations
        if size(unique(DPs(indUSE,:)','rows'),1) < size(DPs,2)
            error('New actUSEvec is the same as a previous one. Aborting search.');
        end

        % Sanity check 3 to ensure the model can grow on the medium
        model_test = set_objective(model_original,{biomassRxn},1,'max',1);
        [~,indDrainsR2open] = ismember(strrep(model.varNames(actUSEvec),'BFUSE_',''),model.varNames);
        indDrainsR2close = setdiff(indDrainsR,indDrainsR2open);
        model_test.var_lb(indDrainsR2close) = 0;
        model_test.var_ub(indDrainsR2close) = 0;
        
        sol_test =  solveTFAmodel_selections(model_test);

        if isempty(sol_test.val) || sol_test.val < 0.99*mu_lb
            error('iMM not working')
        end

        % Save
        save(save_path, 'DPs','model','actUSEvec');

    end

end

end


function model = addIntegerCutConstraint(model, indices, rhs_value)
    [NumCons, NumVars] = size(model.A);
    NewCons = zeros(1, NumVars);
    NewCons(indices) = 1;
    model.A(NumCons + 1, :) = NewCons;
    model.rhs(NumCons + 1) = rhs_value;
    model.constraintNames{NumCons + 1} = ['CUT_' num2str(NumCons)];
    model.constraintType{NumCons + 1} = '>';
end
