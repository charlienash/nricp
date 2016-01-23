%% Demo 2: With missing data

load data/faceSource.mat
load data/faceTargetMissing.mat

% Specify that surface normals are available and can be used.
Options.useNormals = 1;

% Specify that the source deformations should be plotted.
Options.plot = 1;

% Perform non-rigid ICP
[pointsTransformed, X] = nricp(Source, TargetMissing, Options);