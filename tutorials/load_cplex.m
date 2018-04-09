function ret = load_cplex()
%% Edit this with your actual CPLEX path
% cplexPath = 'C:\Program Files\IBM\ILOG\CPLEX_Studio1271\cplex\matlab\x64_win64';
cplexPath = [];

while isempty(cplexPath)
    cplexPath = input('Please provide your cplex path\n','s');
end

addpath(genpath(cplexPath))

[path_found,~] = which('cplex');

path_found
ret = ~isempty(path_found);