function [model,relaxedDGoVarsValues] = FBA2TFA(model,DB,M,M_thermo,flagInfeasibility,minObjSolVal)

    if ~exist("flagInfeasibility","var") || isempty(flagInfeasibility)
        flagInfeasibility = [];
    end
    if ~exist("minObjSolVal","var") || isempty(minObjSolVal)
        minObjSolVal = [];
    end


    %%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare model for thermo 
    ReactionDB = DB;
    CompartmentData = model.CompartmentData;
    replaceData = false;
    verboseFlag = true;
    writeToFileFlag = false;
    model = prepModelforTFA(model, ReactionDB, CompartmentData, replaceData, verboseFlag, writeToFileFlag);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Convert the model to thermo
    rxnNameListNoThermo = [];
    rxnNameListNoDGoRelax = [];
    flagToAddPotentials = false;
    flagToAddLnThermoDisp = false;
    verboseFlag = true;
    printLP = true;
    flagMCA_FarEquilibrium = false;
    bigM_value = M;
    bigM_thermo_value = M_thermo;
    [model, relaxedDGoVarsValues, ~] = convToTFA_selections(model, ReactionDB, rxnNameListNoThermo, flagInfeasibility, rxnNameListNoDGoRelax, minObjSolVal, flagToAddPotentials, flagToAddLnThermoDisp, verboseFlag, printLP, flagMCA_FarEquilibrium,bigM_value,bigM_thermo_value);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Add net flux variables
    model = addNetFluxVariables(model);

end
