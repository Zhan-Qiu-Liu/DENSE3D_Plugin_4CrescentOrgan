function circumferentialParameterize(self, index)
%% Copyright (c) 2016 DENSEanalysis Contributors
% The codes below are modified from the function "parameterizeEndocardialMesh" in "DENSE3D.m".
% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)
% Last Modified: 19:32 June 27, 2018				

%% Circumferential parameterization
	verts = self.dataObj.EndocardialMesh(index).vertices;
	faces = self.dataObj.EndocardialMesh(index).faces;
	apex = self.dataObj.Apex(index,:);

	% For each isoline find the point that is closest to the insertion
	nLines = 50;

	% descending order:
	isolines = linspace(1, 0, nLines + 2);
	isolines([1 end]) = [];

	points = repmat({zeros(0, 3)}, 1, 2);
	cparam = repmat({zeros(0, 1)}, 1, 2);;


	% Find the point in the mesh that corresponds to the apex
	[~, apexind] = ismember(apex, verts, 'rows');
	
	for k = 1:numel(isolines)
		% Draw an isoline at this particular value
		[isoline, tris, ~, ~, SU] = slice_isolines(verts, ...
			faces, self.dataObj.Parameterization(index).Longitudinal(:), isolines(k), 'Manifold', true);

		% Expects that the base is 1 and the apex is 0
		% Do a quick check here to ensure that this is actually the
		% case
		assert(SU(apexind) < 0.5, ...
			'The apex is supposed to be 0 instead of 1 now')

		meshkeep = SU <= (isolines(k) + eps);

		trikeep = all(ismember(tris, find(meshkeep)), 2);
		tris = tris(trikeep, :);

		B = ordered_outline(tris);

		isoline = isoline(B,:);

		% Take into account that PositionB is the right-most so we need
		% to shift our parameterization by 1/6 later on
		% tmp = self.anteriorInsertion(isoline, index, isolines(k));
		tmp = self.anteriorInsertion(mean(isoline), isolines(k));
		if ~isempty(tmp); insertion = tmp; end
		arrayfun(@assignCparam, 1:2);
	end

	% Shift the Circumferential parameterization by 1/6 to account for PositionB
	% being the EDGE of segment 1 and not the LV insertion
	% XXX This is now handled by anteriorInsertion itself
	%cparam = mod(cparam - (1/6), 1);

	% Try to get Circ locations for each vertex in the initial mesh
	ind = dsearchn(points{1}, verts);
	self.dataObj.Parameterization(index).Circumferential = cparam{1}(ind);
	ind = dsearchn(points{2}, verts);
	self.dataObj.Parameterization(index).reverseCircumferential = cparam{2}(ind);
	
	function assignCparam(idxCparam)% assignDir=@lt; if assignDir(D, 0)
		insertion_distances = sum(bsxfun(@minus, insertion(idxCparam,:), isoline).^2, 2);

		[~, startind] = min(insertion_distances);

		% Shift array circularly to make startind as 1st index
		isoline = circshift(isoline, [1-startind, 0]);
		

		% Now parameterize by arc length
		dists = sqrt(sum(diff(isoline([1:end 1],:), [], 1).^2, 2));
		arclength = [0; cumsum(dists(:))];

		arclength = arclength ./ arclength(end);
		arclength(end) = [];

		% Now check whether the orientation of this sampling is
		% correct.

		% Compute the normal vector to the plane based upon the
		% contour orientation
		tmp = normr(bsxfun(@minus, isoline, mean(isoline, 1)));
		mean_norm = mean(cross(tmp(1:end-1,:), tmp(2:end,:)),1);
		D = point2planeDistance(apex, isoline(1,:), mean_norm);

		% D<0: indicate EndoPts are ordered (ascendingly) from septum to freewall. In order to make sure CircumParam ascend adversely (from freewall to septum), reverse the arclength
		% If the apex was not on the "correct" side, then the
		% orientation of the contour was the opposite of what we
		% expected and it needs to be flipped
		if D < 0; arclength = 1 - arclength; end
		if idxCparam < 2
			ind = dsearchn(isoline, insertion(2,:));
			self.hPickSlice.insertion{index}(k,1:2) = [isolines(k),arclength(ind)];
		else
			% if D > 0; arclength = 1 - arclength; end		
			ind = dsearchn(isoline, insertion(1,:));
			self.hPickSlice.insertion{index}(k,3) = arclength(ind);
		end

		points{idxCparam} = cat(1, points{idxCparam}, isoline);
		cparam{idxCparam} = cat(1, cparam{idxCparam}, arclength(:));
	end
end