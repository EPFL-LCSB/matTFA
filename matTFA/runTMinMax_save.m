function TMinMax = runTMinMax_save(model, variables, TimeInSec, save_path)

if nargin < 4 || ~exist('save_path','var')
    save_path = [];
end
num_vars=length(model.var_lb);
model.f = zeros(num_vars,1);
[~,varList] = ismember(variables,model.varNames);
varNames = model.varNames;

% initialization
TMinMax_LB = zeros(size(varList,1),1);
TMinMax_UB = zeros(size(varList,1),1);

if exist(save_path,'file')
    load(save_path);
else
    k = 1;
end

for j = k:length(varList)


    i = ismember(varNames,variables{j});
    model.f = zeros(size(model.f));
    model.f(i) = 1;

    % Min

    model.objtype = 1;
    fprintf('Minimizing %s\n',model.varNames{i});        
    sol_min = solveTFAmodelCplex_selections(model,'TimeInSec',TimeInSec,'emphPar',0,'mipDisplay',0,'barrierDisplay',0,'CPXPARAMdisp',0);
                                                  
    if isempty(sol_min.val)
        TMinMax_LB(j) = NaN;
    else
        TMinMax_LB(j) = sol_min.val; 
    end

    fprintf('Min %d\n',TMinMax_LB(j));        
        
    % Max

    model.objtype = -1;
    fprintf('Maximizing %s\n',model.varNames{i});        
    sol_max = solveTFAmodelCplex_selections(model,'TimeInSec',TimeInSec,'emphPar',0,'mipDisplay',0,'barrierDisplay',0,'CPXPARAMdisp',0);
                                                  
    if isempty(sol_max.val)
        TMinMax_UB(j) = NaN;
    else
        TMinMax_UB(j) = sol_max.val; 
    end

    fprintf('Max %d\n',TMinMax_UB(j));       

    % Save
    
    if rem(j,100) == 0 && ~isempty(save_path)
        save(save_path)
    end

end

TMinMax = [TMinMax_LB TMinMax_UB];

end
