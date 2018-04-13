function [modelSampling,samples] = sampleCbModel_LCSB(model,sampleFile,options)
%sampleCbModel Sample the solution-space of a constraint-based model
%
% [modelSampling,samples] = sampleCbModel(model,sampleFile,samplerName,options)
%
%INPUTS
% model         COBRA model structure
% sampleFile    File names for sampling output files
%
%OPTIONAL INPUTS
% options       Options for sampling and pre/postprocessing (default values
%               in parenthesis)
%   nWarmupPoints           Number of warmup points (5000)
%   nFiles                  Number of output files (10)
%   nPointsPerFile          Number of points per file (1000)
%   nStepsPerPoint          Number of sampler steps per point saved (200)
%   nPointsReturned         Number of points loaded for analysis (2000)
%   nFilesSkipped           Number of output files skipped when loading
%                           points to avoid potentially biased initial samples (2)
%   removeLoopsFlag         Attempt to remove loops from the model (false)
%   removeLoopSamplesFlag   Remove sampling data for reactions involved in
%                           loops (true)
%
%OUTPUTS
% modelSampling Cleaned up model used in sampling
% samples       Uniform random samples of the solution space
%
% Examples of use:
%
% 1)    Sample a model called 'superModel' using default settings and save the
%       results in files with the common beginning 'superModelSamples'
%
%           [modelSampling,samples] = sampleCbModel(superModel,'superModelSamples');
%
% 2)    Sample a model called 'hyperModel' using default settings except with a total of 50 sample files
%       saved and with 5000 sample points returned. 
%       Save the results in files with the common beginning 'hyperModelSamples'
%
%           options.nFiles = 50;
%           options.nPointsReturned = 5000;
%           [modelSampling,samples] = sampleCbModel(hyperModel,'hyperModelSamples');
%
% Markus Herrgard 8/14/06
% Alex Kyparissides has modified the original version of this function. Moreover,
% cobra functions that are modified by LCSB (that are called now within this function) are:
% ACHRSampler_LCSB, and createHRWarmup_LCSB Therefore we
% call it sampleCbModel_LCSB

% Default options
nWarmupPoints = 500;
nFiles = 1;
nPointsPerFile = 5000;
nStepsPerPoint = 100;
nPointsReturned = 5000;
NvZeroTolerance = 1e-8;
nFilesSkipped = 0;
removeLoopsFlag = false;
removeLoopSamplesFlag = true;


% Handle options
if (nargin == 3)
    if (isfield(options,'nWarmupPoints'))
        nWarmupPoints = options.nWarmupPoints;
    end
    if (isfield(options,'nFiles'))
        nFiles = options.nFiles;
    end
    if (isfield(options,'nPointsPerFile'))
        nPointsPerFile = options.nPointsPerFile;
    end
    if (isfield(options,'nStepsPerPoint'))
        nStepsPerPoint = options.nStepsPerPoint;
    end
    if (isfield(options,'nPointsReturned'))
        nPointsReturned = options.nPointsReturned;
    end
    if (isfield(options,'nFilesSkipped'))
        nFilesSkipped = options.nFilesSkipped;
    end
    if (isfield(options,'removeLoopsFlag'))
        removeLoopsFlag = options.removeLoopsFlag;
    end
    if (isfield(options,'removeLoopSamplesFlag'))
        removeLoopSamplesFlag = options.removeLoopSamplesFlag;
    end
        if (isfield(options,'NvZeroTolerance'))
        NvZeroTolerance = options.NvZeroTolerance;
    end
end


fprintf('Prepare model for sampling\n');
% Prepare model for sampling by reducing bounds
[nMet,nRxn] = size(model.S);
fprintf('Original model: %d rxns %d metabolites\n',nRxn,nMet);

modelSampling = model;

        % Use Artificial Centering Hit-and-run
        
        fprintf('Create warmup points\n');
        % Create warmup points for sampler
        warmupPts= createHRWarmup_LCSB(modelSampling,nWarmupPoints,NvZeroTolerance);

%         save modelSampling warmupPts

        fprintf('Run sampler for a total of %d steps\n',nFiles*nPointsPerFile*nStepsPerPoint);
        % Sample model
        ACHRSampler_LCSB(modelSampling,warmupPts,sampleFile,nFiles,nPointsPerFile,nStepsPerPoint,NvZeroTolerance);


fprintf('Load samples\n');
% Load samples
nPointsPerFileLoaded = ceil(nPointsReturned/(nFiles-nFilesSkipped));
if (nPointsPerFileLoaded > nPointsPerFile)
   error('Attempted to return more points than were saved'); 
end
samples = loadSamples(sampleFile,nFiles,nPointsPerFileLoaded,nFilesSkipped);
samples = samples(:,round(linspace(1,size(samples,2),nPointsReturned)));
% Fix reaction directions
[modelSampling,samples] = convRevSamples(modelSampling,samples);
