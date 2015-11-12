function [ landmarks ] = getnodes( Template, radius )
% getNodes gets evenly spaced nodes on mesh for use in an embedded
% deformation graph. Nodes are selected so that no other nodes lie within a
% pre_determined radius.
% 
% Inputs:
%   template : structured object with fields:
%                   template.vertices: N x 3 vertices of template model
%                   template.faces: M x 3 list of connected vertices.
%   radius : controls the spacing of the nodes.
%
% Example:
% model = models(1);
% V = model.vertices;
% N = getnodes(model, 0.05);
% 
% Plot all vertices:
% clf;
% subplot(1,2,1);
% patch(model, 'EdgeColor', 'None', 'FaceColor', 'b');
% l = light('Position', [-1 -1 1]); view([-1 -0.5 0.5]);
% hold on; axis off;
% scatter3(V(:,1), V(:,2), V(:,3), 15, 'r', 'filled');
% axis equal; title('(a)');
%
% Plot nodes:
% subplot(1,2,2);
% patch(model, 'EdgeColor', 'None', 'FaceColor', 'b');
% l = light('Position', [-1 -1 1]); view([-1 -0.5 0.5]);
% hold on; axis off; 
% scatter3(N(:,1), N(:,2), N(:,3), 15, 'r', 'filled'); 
% axis equal; title('(b)');

landmarks = [];
pointsLeft = Template.vertices;
n = size(pointsLeft, 1);
itt = 1;
tol = 0.05;
while size(pointsLeft, 1) > 0
    nPointsLeft = size(pointsLeft, 1);
%     s = sprintf('Iteration: %d  Points left: %d  Landmarks: %d', ...
%                 itt, nPointsLeft, size(landmarks, 1));
%     disp(s);
    pointN = randsample(nPointsLeft, 1);
    point = pointsLeft(pointN, :);
    landmarks(itt,:) = point;
    idx = rangesearch(pointsLeft, point, radius);
    idRemove = idx{1};
    pointsLeft(idRemove, :) = [];
    itt = itt + 1;
end

end

