function computeStrains(self, points)
%% Copyright (c) of the following section of the codes: 2016 DENSEanalysis Contributors
% The codes below are from the function "computeStrains" in "DENSE3D.m".
	% Determine the "ownership" of each point:[indices,distances ]
	[~, epidist] = dsearchn(self.dataObj.EpicardialMesh.vertices, points);

	nEndo = numel(self.dataObj.EndocardialMesh);

	[indices, endodist] = deal(nan(size(points, 1), nEndo));

	for k = 1:nEndo
		[indices(:,k), endodist(:,k)] = dsearchn(self.dataObj.EndocardialMesh(k).vertices, points);
	end

	if self.hShowMesh.isBiv
		% Only need this in the biventricular case
		% true for the samplePoints in septum
		isSeptum = all(bsxfun(@lt, endodist, epidist), 2);
	else
		isSeptum = false(size(epidist));
	end

	% Find the closest endo to each point
	[~, endo_ownership] = min(endodist, [], 2);

	% Now re-assign the septum such that it doesn't belong to any (septum=0, LVendo=1, RVendo=2)
	% particular ventricle
	endo_ownership(isSeptum) = 0;

	function strains = assignStrains(ind)
		% Mapping: any point within the myocardium->the nearest point on the endocardial surface mesh to.=>determine circumferential and longitudinal position using the circumferential and longitudinal parameterization of the endocardial surfaces
		strains = queryStrains(self.dataObj.EndocardialMesh(k), points(ind,:), self.dataObj.Apex(k,:), self.dataObj.Interpolants);

		inds = indices(ind, k);

		% self.dataObj.strains.Parameterization: Assign normalized parameters corresponding to each element of strain components
		param.Circumferential = self.dataObj.Parameterization(k).Circumferential(inds);
		param.reverseCircumferential = self.dataObj.Parameterization(k).reverseCircumferential(inds);
		param.Longitudinal = self.dataObj.Parameterization(k).Longitudinal(inds);

		strains.Parameterization = param;
		strains.Locations = points(ind,:);	
	end

%% Copyright (c) of the following section of the codes: Zhan-Qiu Liu (lafeir.lew@gmail.com)
% Modified By: Zhan-Qiu Liu (lafeir.lew@gmail.com)
% Last Modified: 19:32 June 27, 2018				
	self.dataObj.Strains = [];
	try rmfield(self.hPickSlice, 'Strains'); end
	for k = 1:nEndo
		% Strains([LVendo; RVendo] X [Septum=1, Freewall=2, Global=3])
		% Use the septal datapoints + the ones that belong to
		% the ventricle in this loop
		self.hPickSlice.Strains(k, 1) = assignStrains(isSeptum);
		self.hPickSlice.Strains(k, 2) = assignStrains(endo_ownership == k);
		
		fld = fieldnames(self.hPickSlice.Strains);
		for ii = 1: numel(fld)
			self.hPickSlice.Strains(k, 3).(fld{ii}) = cat(1, self.hPickSlice.Strains(k, 1).(fld{ii}), self.hPickSlice.Strains(k, 2).(fld{ii}));
		end
		param.Circumferential = cat(1, self.hPickSlice.Strains(k, 3).Parameterization.Circumferential);
		param.reverseCircumferential = cat(1, self.hPickSlice.Strains(k, 3).Parameterization.reverseCircumferential);
		param.Longitudinal = cat(1, self.hPickSlice.Strains(k, 3).Parameterization.Longitudinal);
		self.hPickSlice.Strains(k, 3).Parameterization = param;
		
		self.dataObj.Strains = cat(2, self.dataObj.Strains, self.hPickSlice.Strains(k, 3));
	end
% ii = assignStrains(isSeptum | endo_ownership == k);
%% plot sample points for RV endo:
% self.dataObj.preview('Longitudinal'); hold on;
% scatter3(points(touse,1),points(touse,2),points(touse,3),3,'k','*');
	
end