
function ret = load_cplex()
%% Edit this with your actual CPLEX path
%addpath(genpath('C:\Program Files\IBM\ILOG\CPLEX_Studio127\cplex\matlab'));

addpath('/opt/ibm/ILOG/CPLEX_Studio1271/cplex/matlab/x86-64_linux/')

[path_found,~] = which('cplex');

path_found
ret = ~isempty(path_found);