function points = intersectVerts(verts, faces, location, normal, ordered)
% Later Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)
% Last Modified: 17:28 August 7, 2018

	% Find intersected points closed to the slice of interest
	plane = [reshape(normal,1,[]) -dot(location,normal)];
	points = sum(bsxfun(@times,[verts ones(size(verts,1),1)],plane),2);
	intersectFaces = points(faces);
	intersectFaces = reshape(intersectFaces,size(faces));
	facesBelow = sum(intersectFaces<0,2) == 1;
	%% fast but duplicate boundary points:
	%{ 
	facesAbove = sum(intersectFaces>0,2) == 1;
	points = [reshape(find(facesBelow),[],1);reshape(find(facesAbove),[],1)];
	points = verts(unique(faces(points)),:);
	 %}
	
	faces = faces(facesBelow,:);
	intersectFaces = intersectFaces(facesBelow,:);
    [~,J] = min(intersectFaces,[],2);
    I = sub2ind(size(faces),repmat(1:size(faces,1),size(faces,2),1)',mod([J J+1 J+2]-1,3)+1);
    lambda = intersectFaces(I(:,2:3))./bsxfun(@minus,intersectFaces(I(:,2:3)),intersectFaces(I(:,1)));
    BC = sparse( ...
      repmat((1:size(faces,1)*2)',1,2), ...
      [repmat(faces(I(:,1)),2,1) reshape(faces(I(:,2:3)),size(faces,1)*2,1)], ...
      [lambda(:) 1-lambda(:)], ...
      size(faces,1)*2,size(verts,1));
	points = BC * verts;
	% scatter3(points(:,1),points(:,2),points(:,3),10,'b','x');
	
	tolerance = 1e-2*sqrt(sum((max(verts)-min(verts)).^2,2));
	[~,ind] = unique(round(points/tolerance),'rows','stable');
	points = points(ind,:);
	% scatter3(ans(:,1),ans(:,2),ans(:,3),10,'b','x');
	
	if exist('ordered', 'var')
		distances = bsxfun(@minus, points.', permute(points, [2 3 1]));
		distances = squeeze(sum(distances.^2, 1));
		
		no = size(points,1);
		idxPt = 1;
		order = [];
		count = 0;
		while count < no
			order = cat(2, order, idxPt);
			distances(idxPt,:) = deal(NaN);
			[~,ind] = min(distances(:,idxPt), [], 1);
			idxPt = ind;
			count = count + 1;
		end
		points = points(order,:);
	end
end