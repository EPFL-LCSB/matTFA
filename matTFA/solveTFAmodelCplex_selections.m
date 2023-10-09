function sol = solveTFAmodelCplex_selections(tModel, varargin)
% SOLVETFAMODELCPLEX_V3 Solves a TFA model using specific solver settings
%
%   sol = solveTFAmodelCplex_v3(tModel, Name, Value)
%
% Inputs:
%   tModel             - TFA model
%   Name-Value pairs:
%   'manualScalingFactor' - manual scaling factor
%   'mipTolInt'          - mip tolerance integer
%   'emphPar'            - emphasis parameter
%   'feasTol'            - feasibility tolerance
%   'scalPar'            - scaling parameter
%   'TimeInSec'          - time in seconds
%   'mipDisplay'         - mip display
%   'CPXPARAMdisp'       - CPLEX parameter display
%   'timePolishing'      - time for polishing
%   'stop_search_value'  - stop search value
%
%   Note: Make sure the cplex object is saved as a .mst file
%
% Output:
%   sol                 - solution of the TFA model
%
% More details in changeToCPLEX_WithOptions

global TFA_MILP_SOLVER

% Set default values for input parameters if not provided
if ~exist('TFA_MILP_SOLVER', 'var') || isempty(TFA_MILP_SOLVER)
    TFA_MILP_SOLVER = 'cplex_direct';
end
solver = TFA_MILP_SOLVER;

% Define the input parameter names and their default values
paramNames = {'manualScalingFactor', 'mipTolInt','mipTolRel', 'emphPar', 'emphMIPswitch','feasTol', 'scalPar', 'nodeSelec', 'varSelec','heurFreq' , 'TimeInSec', 'mipDisplay', 'barrierDisplay', 'CPXPARAMdisp', 'timePolishing', 'stop_search_value','writeMIPstart','writeMIPpath','loadMIPstart','loadMIPpath'};
defaultValues = {[]                , 1e-9       , 1e-4      , 1        , 0              , 1e-9    , -1       , 1          , 0        , 0          ,[]          , 5          , 2               ,[]             , []              , []                 , 0             ,''            ,0             ,''};

% Parse the input parameters
p = inputParser;
for iParam = 1:numel(paramNames)
    p.addParameter(paramNames{iParam}, defaultValues{iParam});
end
p.parse(varargin{:});

% Choose the appropriate solver and solve the problem
switch solver
    case 'cplex_direct'
        if isempty(which('cplex.p'))
            error('You need to add CPLEX to the Matlab-path!!')
        end
        sol = x_solveCplex(tModel, p.Results);
    case 'gurobi_direct'
        if isempty(which('gurobi'))
            error('You need to add Gurobi to the Matlab-path!!')
        end
        sol = x_solveGurobi(tModel, p.Results);
end
end

% Private function for CPLEX solve
function sol = x_solveCplex(tModel, options)

if isempty(which('cplex.m'))
    error('cplex is either not installed or not in the path')
end

cplex = changeToCPLEX_WithOptions_selections(tModel, options);

if options.loadMIPstart == 1
    cplex.readMipStart(options.loadMIPpath);
end

if ~isempty(options.stop_search_value)
    cplex.InfoCallback.func = @stop_function;
    cplex.InfoCallback.data.IncValue = options.stop_search_value; 
end

try
    cplex.solve();
    if options.writeMIPstart
        cplex.writeMipStart(options.writeMIPpath);
    end
    sol = parseCplexSolution(cplex);

catch
    sol.x = NaN;
    sol.val = NaN;
    sol.cplexSolStatus = 'Solver crashed';
end

delete(cplex)

end

% Callback function
function stop = stop_function(info, data)
    stop = false;
    tolerance = 1e-6;
    
    if ~isempty(info.IncObj)
        if (abs(info.IncObj - data.IncValue) <= tolerance)
            stop = true;
        end
    end
end

% Private function for GUROBI solve
function sol = x_solveGurobi(tModel, options)

gmodel = prepareGurobiModel(tModel);

params = prepareGurobiParams(options);

try
    result = gurobi(gmodel, params);
    sol = parseGurobiSolution(result);
catch
    sol.x = NaN;
    sol.val = NaN;
    sol.status = 'Solver crashed';
end

end


function gmodel = prepareGurobiModel(tModel)
    num_constr = length(tModel.constraintType);
    num_vars = length(tModel.vartypes);

    contypes = '';
    vtypes = '';

    % convert contypes and vtypes into the right format
    for i=1:num_constr
        contypes = strcat(contypes, tModel.constraintType{i,1});
    end

    for i=1:num_vars
        vtypes = strcat(vtypes, tModel.vartypes{i,1});
    end

    gmodel.A = tModel.A;
    gmodel.obj = tModel.f;
    gmodel.lb = tModel.var_lb;
    gmodel.ub = tModel.var_ub;
    gmodel.rhs = tModel.rhs;
    gmodel.sense = contypes;
    gmodel.vtype = vtypes;

    gmodel.varnames = tModel.varNames;

    if tModel.objtype == -1
        gmodel.modelsense = 'max';
    elseif tModel.objtype == 1
        gmodel.modelsense = 'min';
    else
        error(['No objective type specified ' ...
            '(model.objtype should be in {-1,1})']);
    end
end

function params = prepareGurobiParams(options)
    params = struct();
    if ~isempty(options.TimeInSec)
        params.TimeLimit = options.TimeInSec;
    end
    if ~isempty(options.manualScalingFactor)
        params.ScaleFlag = options.manualScalingFactor;
    end
    if ~isempty(options.mipTolInt)
        params.MIPGap = options.mipTolInt;
    end
    if ~isempty(options.emphPar)
        params.MIPFocus = options.emphPar;
    end
    if ~isempty(options.feasTol)
        params.FeasibilityTol = options.feasTol;
    end
    if ~isempty(options.scalPar)
        params.ScaleFlag = options.scalPar;
    end
    if ~isempty(options.mipDisplay)
        params.DisplayInterval = options.mipDisplay;
    end
end

function sol = parseCplexSolution(cplex)
    if isfield(cplex.Solution,'x')
        x = cplex.Solution.x;
        if ~isempty(x)
            sol.x = cplex.Solution.x;
            sol.val = cplex.Solution.objval;
            sol.cplexSolStatus = cplex.Solution.status;
        else
            sol.x = [];
            sol.val = [];
            disp('Empty solution');
            warning('Cplex returned an empty solution!')
            sol.cplexSolStatus = 'Empty solution';
        end
    else
        sol.x = [];
        sol.val = [];
        disp('No field cplex.Solution.x');
        warning('The solver does not return a solution!')
        sol.cplexSolStatus = 'No field cplex.Solution.x';
    end
end

function sol = parseGurobiSolution(result)
    if isfield(result,'x')
        x = result.x;
        x(abs(x) < 1E-9) = 0;
    else
        warning('The solver does not return a solution!');
        result.x = [];
        result.objval = [];
    end

    sol.x = result.x;
    sol.val = result.objval;
    sol.status = result.status;
    % TODO: Add exitflag translation
end


