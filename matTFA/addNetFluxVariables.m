function model = addNetFluxVariables(model)
% adds NetFluxes variables to the model and their associated constraints
%
% The variables will be prefixed with "NF_"
%
% INPUTS:
% - model: a TFA-ready model (has the .A matrix) with no NF_ variables
%
% OUTPUTS
% - model: a TFA-readay model with NF_ variables

[num_mets num_rxns] = size(model.S);
[num_constraints,num_vars] = size(model.A);

% NF_rxn_i: F_rxn - R_rxn - NF_rxn = 0

for i=1:num_constraints
    temp = regexp(model.constraintNames(i),'_','split');
    cons_prefix{i} = temp{1,1}{1};
end

mass_bal_cons = find(ismember(cons_prefix,'M'));
orig_A = model.A;

for i=1:num_rxns
    
    [num_constraints,num_vars] = size(model.A);
    model.constraintNames{num_constraints+1,1} = strcat('NF_',model.rxns{i});
    model.constraintType{num_constraints+1,1} = '=';
    model.rhs(num_constraints+1) = 0;
    
    model.varNames{num_vars+1} = ['NF_' model.rxns{i}];
    model.var_lb(num_vars+1) = -1000;
    model.var_ub(num_vars+1) = 1000;
    model.vartypes{num_vars+1} = 'C';
    model.A(num_constraints+1,num_vars+1) = -1;
    
    F_varIndex = find(ismember(model.varNames,strcat('F_', model.rxns{i})));
    R_varIndex = find(ismember(model.varNames,strcat('R_', model.rxns{i})));
    
    fprintf('replacing %s flux variables with net flux variable\n',model.rxns{i});
    
    if ~isempty(F_varIndex) && ~isempty(R_varIndex)
        % first we remove the separate flux variables from the mass balance constraints and replace
        % them with the net flux variable
        
        F_pos_indices = find(model.A(mass_bal_cons,F_varIndex) > 0);
        F_neg_indices = find(model.A(mass_bal_cons,F_varIndex) < 0);
        R_pos_indices = find(model.A(mass_bal_cons,R_varIndex) > 0);
        R_neg_indices = find(model.A(mass_bal_cons,R_varIndex) < 0);
        
        model.A(F_pos_indices,F_varIndex) = 0;
        model.A(F_neg_indices,F_varIndex) = 0;
        model.A(R_pos_indices,R_varIndex) = 0;
        model.A(R_neg_indices,R_varIndex) = 0;
        
        model.A(F_pos_indices,num_vars+1) = orig_A(F_pos_indices,F_varIndex);
        model.A(F_neg_indices,num_vars+1) = orig_A(F_neg_indices,F_varIndex);
        
        model.A(num_constraints+1,F_varIndex) = 1;
        model.A(num_constraints+1,R_varIndex) = -1;
        
    elseif isempty(F_varIndex)
        R_pos_indices = find(model.A(mass_bal_cons,R_varIndex) > 0);
        R_neg_indices = find(model.A(mass_bal_cons,R_varIndex) < 0);
        
        model.A(R_pos_indices,R_varIndex) = 0;
        model.A(R_neg_indices,R_varIndex) = 0;
        
        model.A(R_pos_indices,num_vars+1) = orig_A(R_pos_indices,R_varIndex);
        model.A(R_neg_indices,num_vars+1) = orig_A(R_neg_indices,R_varIndex);
        
        model.A(num_constraints+1,R_varIndex) = -1;
    elseif isempty(R_varIndex)
        F_pos_indices = find(model.A(mass_bal_cons,F_varIndex) > 0);
        F_neg_indices = find(model.A(mass_bal_cons,F_varIndex) < 0);
        
        model.A(F_pos_indices,F_varIndex) = 0;
        model.A(F_neg_indices,F_varIndex) = 0;
        
        model.A(F_pos_indices,num_vars+1) = orig_A(F_pos_indices,F_varIndex);
        model.A(F_neg_indices,num_vars+1) = orig_A(F_neg_indices,F_varIndex);
        
        model.A(num_constraints+1,F_varIndex) = 1;              
    end
    

end

[num_constraints,num_vars] = size(model.A);
model.f(num_vars,1) = 0;

NF_vi=getAllVar(model,{'NF'});
model.var_lb(NF_vi)=model.lb;
model.var_ub(NF_vi)=model.ub;

end