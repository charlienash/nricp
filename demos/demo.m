%% Complete meshes demo

%% Set working directory
cd /home/charlie/nricp

%% Load data
load data/source.mat
load data/target.mat

%% Specify parameters
Options.gamm = 2;
Options.k = 8;
Options.lambda = 2;
Options.epsilon = 1e-1;

% Stiffness parameters
startValue = 10000;
endValue = 5;
nElements = 50;
Options.alphaSet = linspace(startValue, endValue, nElements);
A = triangulation2adjacency(Source.faces, Source.vertices);
Options.M = adjacency2incidence(A)';

%% nricp
points

