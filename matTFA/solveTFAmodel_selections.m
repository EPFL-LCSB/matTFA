function sol = solveTFAmodel_selections(tModel, varargin)
% SOLVETFAMODEL_SELECTIONS Solves a TFA model using specific solver settings
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
%
% Output:
%   sol                 - solution of the TFA model


%% Initialization

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


%% Private functions for CPLEX solve

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


function sol = parseCplexSolution(cplex)
    if isfield(cplex.Solution,'x')
        x = cplex.Solution.x;
        if ~isempty(x)
            sol.x = cplex.Solution.x;
            sol.val = cplex.Solution.objval;
            sol.SolStatus = cplex.Solution.status;
        else
            sol.x = [];
            sol.val = [];
            disp('Empty solution');
            warning('Cplex returned an empty solution!')
            sol.SolStatus = 'Empty solution';
        end
    else
        sol.x = [];
        sol.val = [];
        disp('No field cplex.Solution.x');
        warning('The solver does not return a solution!')
        sol.SolStatus = 'No field cplex.Solution.x';
    end
end

%% Private functions for GUROBI solve
function sol = x_solveGurobi(tModel, options)

if isempty(which('gurobi'))
    error('gurobi is either not installed or not in the path')
end

[gurobiModel, gurobiParams] = changeToGurobi_WithOptions_selections(tModel, options);


% if options.loadMIPstart == 1
%     gurobiModel.Start = readMIPStart(options.loadMIPpath, gurobiModel);
% end

try
    results = gurobi(gurobiModel, gurobiParams);
%     if options.writeMIPstart
%         writeMIPStart(result, options.writeMIPpath);
%     end
    sol = parseGurobiSolution(results);

catch
    sol.x = NaN;
    sol.val = NaN;
    sol.SolStatus = 'Solver crashed';
end


end



function sol = parseGurobiSolution(results)
    if isfield(results,'x')
        x = results.x;
        if ~isempty(x)
            sol.x = results.x;
            sol.val = results.objval;
            sol.SolStatus = lower(results.status);
        else
            sol.x = [];
            sol.val = [];
            disp('Empty solution');
            warning('Gurobi returned an empty solution!')
            sol.SolStatus = 'Empty solution';
        end
    else
        sol.x = [];
        sol.val = [];
        disp('No field results.Solution.x');
        warning('The solver does not return a solution!')
        sol.SolStatus = 'No field results.x';
    end
end