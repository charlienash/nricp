%% Complete meshes demo

%% Set working directory
cd /home/charlie/nricp

%% Load data
load data/faceSource.mat
load data/faceTarget.mat

%% Specify parameters
Options.gamm = 1;
Options.epsilon = 1e-8;

% Stiffness parameters
startValue = 100;
endValue = 1e-3;
nElements = 200;
Options.alphaSet = linspace(startValue, endValue, nElements);
A = triangulation2adjacency(Source.faces, Source.vertices);
Options.M = adjacency2incidence(A)';

%% nricp
correspondences = nricp( Source, Target, Options );


