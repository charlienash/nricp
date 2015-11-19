%% Demo 1: Use surface normals

% Load data
load data/faceSource.mat
load data/faceTarget.mat

% Specify that surface normals are available and can be used.
Options.useNormals = 1;

% Specify that the source deformations should be plotted.
Options.plot = 1;

% Perform non-rigid ICP
[pointsTransformed, X] = nricp(Source, Target, Options);
