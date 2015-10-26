function [ pointsCor, ittErr ] = nricp( Source, Target, Options )
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
nVertsTemplate = size(Source.vertices, 1);
% nPointsTarget = size(pointsTarget, 1);
pointsCor = Source.vertices;
pointsTarget = Target.vertices;
nPointsTarget = size(pointsTarget, 1);


%% Set parameters 
G = diag([1 1 1 Options.gamm]);

%% Plot model
clf;
patch(Target, 'facecolor', 'b', 'EdgeColor', ...
      'none', 'FaceAlpha', 0.5);
hold on;
patch(Source, 'facecolor', 'r', 'EdgeColor', ...
      'none', 'FaceAlpha', 0.5);
view([30,20]);
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
axis equal; 
axis manual;
legend('Target', 'Source', 'Location', 'best')

%% Non Rigid ICP
% s = sprintf('Starting non-rigid deformation... \n\nIteration   Error');
% disp(s);
itt = 1;
nAlpha = numel(Options.alphaSet);
ittErr = zeros(nAlpha, 3);
oldError = 1000;
newError = Inf;

for i = 1:nAlpha
    
    % Update stiffness
    alpha = Options.alphaSet(i);
    s = sprintf('\nalpha = %3.3f\nIteration  Error', alpha);
    disp(s);
    
    while norm(oldError - newError) >= Options.epsilon

        % keep track of old error
        oldError = newError;
    
        % Specify D_temp
        X = zeros(nVertsTemplate, 4 * nVertsTemplate);
        for i = 1:nVertsTemplate
            X(i,(4 * i-3):(4 * i)) = [pointsCor(i,:) 1];
        end

        % Determine closest points on target U 
        [idTemp, distTemp] = knnsearch(Target.vertices, pointsCor);
        U = Target.vertices(idTemp,:);
        
        % Specify W 
        % 1.(weight according to correspondence distance)
%         sim = 1 ./ (distTemp + 1);
%         W = diag(sim ./ max(sim));
        
        % 2. No weighting
        W = eye(nVertsTemplate);
        
        % Set target landmarks tarU
        tarU = Target.vertices;
        
        % Specify tarD
        [idTarget, distTar] = knnsearch(pointsCor, tarU);
        tarX = zeros(nPointsTarget, 4 * nVertsTemplate);
        for i = 1:nPointsTarget
            cor = idTarget(i);
            tarX(i,(4 * cor-3):(4 * cor)) = [pointsCor(cor,:) 1];
        end
        
        % Specify tarW
        % 1.(weight according to correspondence distance)
%         simTar = 1 ./ (distTar + 1);
%         tarW = diag(sim ./ max(simTar));
        
        % 2. No weighting
        tarW = eye(nPointsTarget);

        % Specify B and C
        B = [...
            alpha .* kron(Options.M, G); 
            W * X;
            Options.lambda .* tarW * tarX;
            ];
        C = [...
            zeros(size(alpha .* kron(Options.M, G), 1), 3);
            W * U;
            Options.lambda .* tarW * tarU
            ];

        % Get A
        A = (B' * B) \ (B' * C);

        % Update nodes
        pointsCor = X * A;

        % Error
        errorData = norm(W * (X * A - U), 'fro')^2 + ...
                    norm(tarW * (tarX * A - tarU), 'fro')^2;
%         errorStiffness = norm(kron(M, G) * X, 'fro')^2;
        errorTotal = norm(B * A - C, 'fro')^2;
        newError = errorData;

        % Store iteration and error
        ittErr(itt,:) = [itt, errorData, errorTotal];

        % Print iteration and error_test
        s = sprintf('%2d          %5.3f', itt, errorData);
        disp(s);
        itt = itt + 1;

        % Update plot
        set(h, ...
            'XData', pointsCor(:, 1), ...
            'YData', pointsCor(:, 2), ...
            'ZData', pointsCor(:, 3));
        pause(0.01);
    end
    oldError = 1000 .* newError;
end

%% Snap nodes to closest corresponding points on target
idx = knnsearch(Target.vertices, pointsCor, 'K', 1);
% pointsCor = Target.vertices(idx,:);
set(h, ...
    'XData', pointsCor(:, 1), ...
    'YData', pointsCor(:, 2), ...
    'ZData', pointsCor(:, 3));
pause(0.01);

end

