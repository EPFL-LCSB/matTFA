function TMinMax = runTMinMax_save_selections(model, variables, TimeInSec, save_path)

if nargin < 3 || ~exist('TimeInSec','var')
    TimeInSec = 3*60;
end
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
    k = j+1;
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
    sol_min = solveTFAmodel_selections(model,'TimeInSec',TimeInSec,'emphPar',1,'mipDisplay',5,'barrierDisplay',0,'CPXPARAMdisp',0);
                                                  
    if isempty(sol_min.val)
        sol_min = solveTFAmodel_selections(model,'TimeInSec',TimeInSec,'emphPar',0,'mipDisplay',5,'barrierDisplay',0,'CPXPARAMdisp',0);
        if isempty(sol_min.val)
            model.c = zeros(size(model.c));
            model.c(j) = 1;
            original_osense = model.osense; 
            model.osense = 1;
            sol_min_FBA = solveFBAmodel_selections(model);
            model.cense = original_osense;      
            TMinMax_LB(j) = sol_min_FBA.f; 
        else
            TMinMax_LB(j) = sol_min.val; 
        end            
    else
        TMinMax_LB(j) = sol_min.val; 
    end

    fprintf('Min %d\n',TMinMax_LB(j));        
        
    % Max

    model.objtype = -1;
    fprintf('Maximizing %s\n',model.varNames{i});        
    sol_max = solveTFAmodelCplex_selections(model,'TimeInSec',TimeInSec,'emphPar',1,'mipDisplay',5,'barrierDisplay',0,'CPXPARAMdisp',0);

    if isempty(sol_max.val)
        sol_max = solveTFAmodelCplex_selections(model,'TimeInSec',TimeInSec,'emphPar',0,'mipDisplay',5,'barrierDisplay',0,'CPXPARAMdisp',0);
        if isempty(sol_max.val)
            %TMinMax_UB(j) = NaN;
            model.c = zeros(size(model.c));
            model.c(j) = 1;
            original_osense = model.osense; 
            model.osense = 1;            
            sol_max_FBA = solveFBAmodelCplex(model);
            model.cense = original_osense;                  
            TMinMax_UB(j) = sol_max_FBA.f; 
        else
            TMinMax_UB(j) = sol_max.val; 
        end            
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
