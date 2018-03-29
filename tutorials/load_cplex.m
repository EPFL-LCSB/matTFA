function ret = load_cplex()
%% Edit this with your actual CPLEX path
%addpath(genpath('C:\Program Files\IBM\ILOG\CPLEX_Studio127\cplex\matlab'));

cplexPath = [];
while isempty(cplexPath)
    cplexPath = input('Please provide your cplex path\n','s');
end

addpath(genpath(cplexPath))

[path_found,~] = which('cplex');

path_found
ret = ~isempty(path_found);