%% Demo 2: Specify parameters

% Load data
load data/boneSource.mat
load data/boneTarget.mat

% Specify bidirectional distance metric
Options.lambda = 100;
Options.bidirectional = 1;

% Perform non-rigid iterative closest point
[ pointsTransformed, X ] = nricp(Source, Target, Options);