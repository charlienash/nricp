verts = Target.vertices;
normals = Target.normals;
faceMap = [find(verts(:,2) < 0)'; 1:numel(find(verts(:,2) < 0))]';
newVerts = verts(verts(:,2) < 0, :);
newNormals = normals(verts(:,2) < 0, :);

flatFaces = reshape(Target.faces, numel(Target.faces), 1);
flatFacesNew = zeros(size(flatFaces));
% Flattens by cols
for v=1:numel(flatFaces)
    oldVert = flatFaces(v);
    newVert = faceMap(faceMap(:,1) == oldVert, 2);
    if isempty(newVert)
        flatFacesNew(v) = NaN;
    else
        flatFacesNew(v) = newVert;
    end
    disp(v);
end
newFaces = reshape(flatFacesNew, numel(flatFacesNew)/3, 3);
newFaces = newFaces(~any(isnan(newFaces),2),:); 

TargetMissing.vertices = newVerts;
TargetMissing.faces = newFaces;
TargetMissing.normals = newNormals;

save('data/faceTargetMissing.mat', 'TargetMissing')