function [ vertsTransformed, X ] = nricp( Source, Target, Options )
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
%   Source: structured object with fields - 
%       Source.vertices: V x 3 vertices of template model
%       Source.faces: F x 3 list of connected vertices.
%       Source.normals: (Optional) FV x 3 list of surface normals. Make
%           sure to set Options.normals = 1 if using normals.
% 
%   Target : stuctured object as above for target model.
% 
%   Options : structured object with fields:
%       gamm : real valued, weights differences in the rotational and skew
%           part of the deformation against the translational part.
%       epsilon : real values, tolerence for change in transformation.
%       lambda : If using the bi-directional distance metric this weights
%           the contribution of the target -> source term.
%       alphaSet : decreasing vector of real-valued stiffness parameters. 
%           High stiffness parameters force global transformations whereas 
%           low values allow for local deformations.
%       biDirectional : logical, specifies that a bi-directional distance 
%           is used.
%       useNormals : logical, specifies that surface normals are to be used
%           to project the source onto the target surface. If this term is
%           used then the Source input should contain a normals field.
%       plot : logical, specifies that the transformations should be
%           plotted.
%       rigidInit : logical, specifies that rigid ICP should be performed
%           first before allowing non-rigid and non-global deformations.
%
% Outputs:
%   vertsTransformed : N X 3 vertices of transformed source mesh,
%   X : (4N) X 3 stacked matrix of transformations.

% Set default parameters
if ~isfield(Options, 'gamm')
    Options.gamm = 1;
end
if ~isfield(Options, 'epsilon')
    Options.epsilon = 1e-4;
end
if ~isfield(Options, 'lambda')
    Options.lambda = 1;
end
if ~isfield(Options, 'alphaSet')
    Options.alphaSet = linspace(100, 10, 20);
end
if ~isfield(Options, 'biDirectional')
    Options.biDirectional = 0;
end
if ~isfield(Options, 'useNormals')
    Options.useNormals = 0;
end
if ~isfield(Options, 'plot')
    Options.plot = 0;
end
if ~isfield(Options, 'rigidInit')
    Options.rigidInit = 1;
end
if ~isfield(Options, 'ignoreBoundary')
    Options.ignoreBoundary = 1;
end
if ~isfield(Options, 'normalWeighting')
    Options.normalWeighting = 1;
end

% Optionally plot source and target surfaces
if Options.plot == 1
    clf;
    PlotTarget = rmfield(Target, 'normals');
    p = patch(PlotTarget, 'facecolor', 'b', 'EdgeColor',  'none', ...
              'FaceAlpha', 0.5);
    hold on;
    
    PlotSource = rmfield(Source, 'normals');
    h = patch(PlotSource, 'facecolor', 'r', 'EdgeColor',  'none', ...
        'FaceAlpha', 0.5);
    material dull; light; grid on; xlabel('x'); ylabel('y'); zlabel('z');
    view([60,30]); axis equal; axis manual;
    legend('Target', 'Source', 'Location', 'best')
    drawnow;
end

% Get source vertices 
vertsSource = Source.vertices;
nVertsSource = size(vertsSource, 1);

% Get target vertices
vertsTarget = Target.vertices;

% Optionally get source / target normals
if Options.normalWeighting == 1
    normalsSource = Source.normals;
    normalsTarget = Target.normals;
end

% Get subset of target vertices if Options.biDirectional == 1
if Options.biDirectional == 1
    samplesTarget = sampleVerts(Target, 15);
    nSamplesTarget = size(samplesTarget, 1);
end

% Set matrix G (equation (3) in Amberg et al.) 
G = diag([1 1 1 Options.gamm]);

% Set incidence matrix M 
A = triangulation2adjacency(Source.faces, Source.vertices);
M = adjacency2incidence(A)';

% Precompute kronecker product of M and G
kron_M_G = kron(M, G);

% Set matrix D (equation (8) in Amberg et al.)
I = (1:nVertsSource)';
J = 4*I;
D = sparse([I;I;I;I],[J-3;J-2;J-1;J],[vertsSource(:);ones(nVertsSource,1)],nVertsSource, 4*nVertsSource);

% Set weights vector
wVec = ones(nVertsSource,1);

% Get boundary vertex indices on target surface if required.
if Options.ignoreBoundary == 1
    bdr = find_bound(vertsTarget, Target.faces);
end

% Set target points matrix tarU and target weights matrix tarU
% if Options.biDirectional == 1.
if Options.biDirectional == 1
    tarU = samplesTarget;
    tarW = eye(nSamplesTarget);
end

% Do rigid iterative closest point if Options.rigidInit == 1
if Options.rigidInit == 1
    disp('* Performing rigid ICP...');
    if Options.ignoreBoundary == 0
        bdr = 0;
    end
    [R, t] = icp(vertsTarget', vertsSource', 50, 'Verbose', true, ...
                 'EdgeRejection', logical(Options.ignoreBoundary), ...
                 'Boundary', bdr', 'Matching', 'kDtree');
    X = repmat([R'; t'], nVertsSource, 1);
    vertsTransformed = D*X;
    
    % Update plot
    if Options.plot == 1
        set(h, 'Vertices', vertsTransformed);
        drawnow;
    end
else
    % Otherwise initialize transformation matrix X with identity matrices
    X = repmat([eye(3); [0 0 0]], nVertsSource, 1);
end

% get number of element in the set of stiffness parameters Options.alphaSet
nAlpha = numel(Options.alphaSet);

% Enter outer loop of the non-rigid iterative closest point algorithm. The
% outer loop iterates over stiffness parameters alpha.
disp('* Performing non-rigid ICP...');
for i = 1:nAlpha
    
    % Update stiffness
    alpha = Options.alphaSet(i);
    
    % set oldX to be very different to X so that norm(X - oldX) is large on 
	% first iteration
	oldX = 10*X;
    
    % Enter inner loop. For each stiffness setting alternate between 
    % updating correspondences and getting optimal transformations X. 
    % Break the loop when consecutive transformations are similar.
    while norm(X - oldX) >= Options.epsilon 
        
        % Transform source points by current transformation matrix X
        vertsTransformed = D*X;
        
        % Update plot 
        if Options.plot == 1
            set(h, 'Vertices', full(vertsTransformed));
            drawnow;
        end

        % Determine closest points on target U to transformed source points
        % pointsTransformed.
        targetId = knnsearch(vertsTarget, vertsTransformed);
        U = vertsTarget(targetId,:);
        
        % Optionally give zero weightings to transformations associated
        % with boundary target vertices.
        if Options.ignoreBoundary == 1
            tarBoundary = ismember(targetId, bdr);
            wVec = ~tarBoundary;
        end
        
        % Optionally transform surface normals to compare with target and
        % give zero weight if surface and transformed normals do not have
        % similar angles.
        if Options.normalWeighting == 1
            I = (1:nVertsSource)';
            J = 4*I;
            N = sparse([I;I;I;I],[J-3;J-2;J-1;J],[normalsSource(:);ones(nVertsSource,1)],nVertsSource, 4*nVertsSource);

            normalsTransformed = N*X;
            corNormalsTarget = normalsTarget(targetId,:);
            crossNormals = cross(corNormalsTarget, normalsTransformed);
            crossNormalsNorm = sqrt(sum(crossNormals.^2,2));
            dotNormals = dot(corNormalsTarget, normalsTransformed, 2);
            angle = atan2(crossNormalsNorm, dotNormals);
            wVec = wVec .* (angle<pi/4);
        end
            
        % Update weight matrix
        W = spdiags(wVec, 0, nVertsSource, nVertsSource);

        % Get closest points on source tarD to target samples samplesTarget
        if Options.biDirectional == 1
            transformedId = knnsearch(vertsTransformed, samplesTarget);
            tarD = sparse(nSamplesTarget, 4 * nVertsSource);
            for j = 1:nSamplesTarget
                cor = transformedId(j);
                tarD(j,(4 * cor-3):(4 * cor)) = [vertsSource(cor,:) 1];
            end
        end

        % Specify B and C (See equation (12) from paper)
        A = [...
            alpha .* kron_M_G; 
            W * D;
            ];
        B = [...
            zeros(size(M,1)*size(G,1), 3);
            W * U;
            ];
        
        % Concatentate additional terms if Options.biDirectional == 1.
        if Options.biDirectional == 1
            A = [...
                A;
                Options.lambda .* tarW * tarD
                ];
            B = [...
                B;
                Options.lambda .* tarW * tarU
                ];
        end

        % Get optimal transformation X and remember old transformation oldX
        oldX = X;
        X = (A' * A) \ (A' * B);
    end
end

% Compute transformed points 
vertsTransformed = D*X;

% If Options.useNormals == 1 project along surface normals to target 
% surface, otherwise snap to closest points on target.
if Options.useNormals == 1
    disp('* Projecting transformed points onto target along surface normals...');
    
    % Get template surface normals
    normalsTemplate = Source.normals;
    
    % Transform surface normals with the X matrix
    I = (1:nVertsSource)';
    J = 4*I;
    N = sparse([I;I;I;I],[J-3;J-2;J-1;J],[normalsTemplate(:);ones(nVertsSource,1)],nVertsSource, 4*nVertsSource);

    normalsTransformed = N*X;
    
    % Project normals to target surface
    vertsTransformed = projectNormals(vertsTransformed, Target, ...
                                       normalsTransformed);
else
    % Snap template points to nearest vertices on surface
    targetId = knnsearch(vertsTarget, vertsTransformed);
    corTargets = vertsTarget(targetId,:);
    if Options.ignoreBoundary == 1
        tarBoundary = ismember(targetId, bdr);
        wVec = ~tarBoundary;
    end
    vertsTransformed(wVec,:) = corTargets(wVec,:);
end

% Update plot and remove target mesh
if Options.plot == 1
    set(h, 'Vertices', vertsTransformed);
    drawnow;
    pause(2);
%     delete(p);
end

function [projections] = projectNormals(sourceVertices, Target, normals)
% projectNormals takes a set of vertices and their surface normals and
% projects them to a target surface.
%
% Inputs:
%   sourceVertices: N x 3 vertices of source surface.
%   Target: Structured object with fields - 
%                   Target.vertices: V x 3 vertices of template model
%                   Target.faces: D x 3 list of connected vertices.
%   normals: N x 3 surface normals of sourceVertices.

% Get number of source vertices
nVerticesSource = size(sourceVertices, 1);

% Pre-allocate space for projections
projections = zeros(nVerticesSource, 3);

% Loop over source vertices projecting onto the target surface
for i=1:nVerticesSource
    
    % Get vertex and normal
    vertex = sourceVertices(i,:);
    normal = normals(i,:);
    
    % Define line in direction normal that passes through vertex
    line = createLine3d(vertex, normal(1), normal(2), normal(3));
    
    % Compute the intersection of the line and the source surface
    intersection = intersectLineMesh3d(line, Target.vertices, Target.faces); 
    
    % If multiple intersections choose the one closest to the source vertex
    if ~isempty(intersection)
        [~,I] = min(sqrt(sum((intersection - ...
            repmat(vertex,size(intersection,1),1)).^2, 2)));
        projections(i,:) = intersection(I,:);
    else
        % If no intersections just keep the source vertex position
        projections(i,:) = vertex;
    end
end


function [ samples ] = sampleVerts( Mesh, radius )
% sampleVerts sub samples the vertices of a mesh. Vertices are selected 
% so that no other nodes lie within a pre-determined radius.
% 
% Inputs:
%   Mesh : structured object with fields:
%                   Mesh.vertices: N x 3 vertices of Mesh.
%                   Mesh.faces: M x 3 list of connected vertices.
%   radius : controls the spacing of the vertices.

samples = [];
vertsLeft = Mesh.vertices;
itt = 1;
while size(vertsLeft, 1) > 0
    nVertsLeft = size(vertsLeft, 1);
    
    % pick a sample from remaining points
    vertN = randsample(nVertsLeft, 1);
    vert = vertsLeft(vertN, :);
    
    % Add it to sample set
    samples(itt,:) = vert;
    
    % Remove nearby vertices
    idx = rangesearch(vertsLeft, vert, radius);
    idRemove = idx{1};
    vertsLeft(idRemove, :) = [];
    
    itt = itt + 1;
end


function bound = find_bound(pts, poly)
%  From Iterative Closest Point: 
% http://www.mathworks.com/matlabcentral/fileexchange/27804-iterative-closest-point

% Boundary point determination. Given a set of 3D points and a
% corresponding triangle representation, returns those point indices that
% define the border/edge of the surface.
% Correcting polygon indices and converting datatype 

poly = double(poly);
pts = double(pts);

%Calculating freeboundary points:
TR = triangulation(poly, pts);
FF = freeBoundary(TR);

%Output
bound = FF(:,1);
