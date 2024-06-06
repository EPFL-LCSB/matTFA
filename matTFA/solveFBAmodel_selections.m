function sol = solveFBAmodel_selections(model, varargin)

%% Initialization

global TFA_LP_SOLVER

% Set default values for input parameters if not provided
if ~exist('TFA_LP_SOLVER', 'var') || isempty(TFA_LP_SOLVER)
    TFA_LP_SOLVER = 'cplex_direct';
end
solver = TFA_LP_SOLVER;

% Define the input parameter names and their default values
paramNames = {'scalPar', 'emphPar','feasTol', 'TimeInSec' };
defaultValues = { 1e-9  , 1        ,  1e-9    , 3600    };

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
        sol = x_solveCplex_FBA(model, p.Results);
    case 'gurobi_direct'
        if isempty(which('gurobi'))
            error('You need to add Gurobi to the Matlab-path!!')
        end
        sol = x_solveGurobi_FBA(model, p.Results);
end
end


%% Private functions for CPLEX solve

function sol = x_solveCplex_FBA(model, options)

if isempty(which('cplex.m'))
    error('cplex is either not installed or not in the path')
end

cplex = changeToCPLEX_FBA_WithOptions_selections(model, options);

try
    cplex.solve();
    sol = parseCplexSolution(cplex);

catch
    sol.x = NaN;
    sol.f = NaN;
    sol.SolStatus = 'Solver crashed';
end

delete(cplex)

end

function sol = parseCplexSolution(cplex)
    if isfield(cplex.Solution,'x')
        x = cplex.Solution.x;
        if ~isempty(x)
            sol.x = cplex.Solution.x;
            sol.f = cplex.Solution.objval;
            sol.SolStatus = cplex.Solution.status;
        else
            sol.x = [];
            sol.f = [];
            disp('Empty solution');
            warning('Cplex returned an empty solution!')
            sol.SolStatus = 'Empty solution';
        end
    else
        sol.x = [];
        sol.f = [];
        disp('No field cplex.Solution.x');
        warning('The solver does not return a solution!')
        sol.SolStatus = 'No field cplex.Solution.x';
    end
end



%% Private functions for GUROBI solve
function sol = x_solveGurobi_FBA(model, options)

if isempty(which('gurobi'))
    error('gurobi is either not installed or not in the path')
end

[gurobiModel, gurobiParams] = changeToGurobi_FBA_WithOptions_selections(model, options);


try
    results = gurobi(gurobiModel, gurobiParams);
    sol = parseGurobiSolution(results);

catch
    sol.x = NaN;
    sol.f = NaN;
    sol.SolStatus = 'Solver crashed';
end


end

function sol = parseGurobiSolution(results)
    if isfield(results,'x')
        x = results.x;
        if ~isempty(x)
            sol.x = results.x;
            sol.f = results.objval;
            sol.SolStatus = lower(results.status);
        else
            sol.x = [];
            sol.f = [];
            disp('Empty solution');
            warning('Gurobi returned an empty solution!')
            sol.SolStatus = 'Empty solution';
        end
    else
        sol.x = [];
        sol.f = [];
        disp('No field results.Solution.x');
        warning('The solver does not return a solution!')
        sol.SolStatus = 'No field results.x';
    end
end
