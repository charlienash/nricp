%% Complete meshes demo

%% Load data
load source.mat
load target.mat

%% Specify parameters
Options.gamm = 2;
Options.k = 8;
Options.lambda = 2;
Options.epsilon = 1e-1;

% Stiffness parameters
startValue = 100;
endValue = 5;
nElements = 50;
Options.alphaSet = linspace(startValue, endValue, nElements);
A = triangulation2adjacency(Source.faces, Source.vertices);
% Options.M = triangulation2adjacency(Source.faces, Source.vertex)


