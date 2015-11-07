function [ pointsTemplate, ittErr ] = nricp( Source, Target, Options )
% nricp performs an adaptive stiffness variant of non rigid ICP.
%
% This function deforms takes a dense set of landmarks points from a template
% template model and finds a deformation which matches a target shape. 
% The deformations are encouraged to be natural and smooth by means of a 
% stiffness constraint, which is relaxed in increments.
% 
% For details on the stiffness constraint and optimization procedure see: 
% 'Optimal Step Nonrigid ICP Algorithms for Surface Registration', 
% Amberg, Romandhani and Vetter, CVPR, 2007.
%
% Inputs:
%   nodesInit : N x 3 vertices of template nodes
%   Target : stuctured object as above for target model.
%   Options : structure object with fields:
%                   gamma : real valued, weights differences in the 
%                           rotational and skew part of the
%                           deformation against the translational part 
%                           of the deformation.
%                   alphaSet : decreasing vector of stiffness parameters.
%                   k : number of nearest neighbours used for 
%                       deformation graph. 
%
% Outputs:
%   nodeCor : J X 3 corresponding nodes on target.
%   ittErr : Vector of errors at each iteration.
%
% Example:

%% Get nPoints
normalsTemplate = Source.normals;
pointsTemplate = Source.vertices;
nPointsTemplate = size(pointsTemplate, 1);
pointsTarget = Target.vertices;

%% Set parameters 
G = diag([1 1 1 Options.gamm]);

%% Plot model
TargetPatch.faces = Target.faces;
TargetPatch.vertices = Target.vertices;
SourcePatch.faces = Source.faces;
SourcePatch.vertices = Source.vertices;
clf;
patch(TargetPatch, 'facecolor', 'b', 'EdgeColor', ...
      'none', 'FaceAlpha', 0.5);
hold on;
h = patch(SourcePatch, 'facecolor', 'r', 'EdgeColor', 'none', ... 
'FaceAlpha', 0.5);
material dull;
light;
grid on; 
hold on;
xlabel('x');
ylabel('y')
zlabel('z')
% h = scatter3(pointsTemplate(:,1), ...
%              pointsTemplate(:,2), ...
%              pointsTemplate(:,3), ...
%              'r', 'filled');
% scatter3(pointsTarget(:,1), ...
%          pointsTarget(:,2), ...
%          pointsTarget(:,3), ...
%          'g', 'filled');
view([90,0]);
axis equal; 
axis manual;
legend('Target', 'Source', 'Location', 'best')

%% Do pre-computations
D = sparse(zeros(nPointsTemplate, 4 * nPointsTemplate));
for i = 1:nPointsTemplate
    D(i,(4 * i-3):(4 * i)) = [pointsTemplate(i,:) 1];
end
W = sparse(eye(nPointsTemplate));
N = sparse(zeros(nPointsTemplate, 4 * nPointsTemplate));
for i = 1:nPointsTemplate
    N(i,(4 * i-3):(4 * i)) = [normalsTemplate(i,:) 1];
end

%% Non Rigid ICP
nAlpha = numel(Options.alphaSet);
X = repmat([eye(3); [0 0 0]], nPointsTemplate, 1);
oldX = 10*X;

for i = 1:nAlpha
    
    % Update stiffness
    alpha = Options.alphaSet(i);
    
    while norm(X - oldX) >= Options.epsilon
        
        % Update nodes
        pointsTransformed = D*X;
        
        % Update plot
        set(h, 'Vertices', pointsTransformed);
        drawnow;

        % Determine closest points on target U 
        idTemp = knnsearch(pointsTarget, pointsTransformed);
        U = pointsTarget(idTemp,:);

        % Specify B and C
        B = [...
            alpha .* kron(Options.M, G); 
            W * D;
            ];
        C = [...
            zeros(size(Options.M,1)*size(G,1), 3);
            W * U;
            ];

        % Get X and remember oldX
        oldX = X;
        X = (B' * B) \ (B' * C);
        
    end
end

%% Compute transformed points and normals
pointsTransformed = D*X;
normalsTransformed = N*X;

%% Snap nodes to closest corresponding points on target
% idx = knnsearch(pointsTarget, pointsTransformed);
% pointsTransformed = pointsTarget(idx,:);
%         set(h, 'Vertices', pointsTransformed);
%         drawnow;

end

