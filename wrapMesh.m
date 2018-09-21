function meshes = wrapMesh(meshes)
%% Copyright (c) 2016 DENSEanalysis Contributors
% The codes below are modified from the script "surfacemesh.m".

	verts = reshape(meshes.vertices(meshes.faces.', :), 3, [], 3);

	% Compute the surface normals using two edges
	e1 = squeeze(verts(2,:,:) - verts(1,:,:));
	e2 = squeeze(verts(3,:,:) - verts(1,:,:));
	meshes.normals = normr(cross(e1, e2));

	% Make sure that these normals point towards the centroid
	meshes.centroids = squeeze(mean(verts, 1));
	cdiff = bsxfun(@minus, meshes.centroids, mean(meshes.centroids, 1));

	direction = sign(dot(cdiff, meshes.normals, 2));

	% If the normals point away from the centroid (overall) then
	% flip them all
	if mode(direction) == 1
		meshes.normals = -meshes.normals;
	end
	
	% the struct of "meshes" got to have the fields in the order of 'faces', 'vertices', 'normals', 'centroids'
	meshes = struct('faces', meshes.faces, 'vertices', meshes.vertices, 'normals', meshes.normals, 'centroids', meshes.centroids);
end