function mat2ParaView
% mat2ParaView  - transform .dns inFiles of SA and LA
%
%
%	transform .dns inFiles from Version 0.4 2009.06 to Version: 0.5
%
%
% USAGE:
%   mat2ParaView();
%
%
% See also:
%   
%
% Last Modified: 12:04 June 4, 2019
% Copyright (c) Zhan-Qiu Liu (lafeir.lew@gmail.com)
	import plugins.dense3D_plugin_4crescentorgan.*
	
	homepath = chars(java.lang.System.getProperty('user.home'));
	pathDNS = fullfile(homepath,'Dropbox','Analysis_Human')%Directory where DNS stored by default
	pathMAT = fullfile(homepath,'Dropbox','MAT_Human')%Directory where DENSE2D outputs(post-analyzed .MAT) stored by default

	while true
		[uifile, uipath] = uigetfile({'*.mat','DENSE2D(v0.4) outputs(post-analyzed *.MAT)'},['Select MAT file for SA ROI of the input:'], pwd);
		if ~isequal(uipath,0) % && ~isequal(uifile,0)
			breaks
		end
		pathAnalysisMAT = fullfile(uipath,uifile);
	end
	load(pathAnalysisMAT, '-mat') % 'AnalysisInfo','DisplacementInfo', 'ImageInfo', 'SequenceInfo',
	
	%%% Displacement Vectors:
	NodesRef = [DisplacementInfo.X DisplacementInfo.Y];
	NodesRef(:,3) = deal(0);
	flagCINE = true;
	
	%% Displacement Vectors w.r.t. ref. frame:
	disps = sqrt(DisplacementInfo.dX.^2 + DisplacementInfo.dY.^2 + DisplacementInfo.dZ.^2);
	frameStart = find(~isnan(DisplacementInfo.dX(1,:)),1); % 37
	for t = frameStart:DENSEInfo.Number
		Nodes = NodesRef + [DisplacementInfo.dX(:,t) DisplacementInfo.dY(:,t) DisplacementInfo.dZ(:,t)];
		if flagCINE; Nodes = tformfwd(SequenceInfo(1).tform, Nodes); end
		VectorsToVTK(Nodes, disps(:,t), 'Displacements', [DisplacementInfo.dX(:,t),DisplacementInfo.dY(:,t),DisplacementInfo.dZ(:,t)], 'DisplacementVectors', fullfile(uipath,['DENSEdisp',num2str(t)]))
	end

	%% Displacement Vectors w.r.t. a specific frame:
temp(1,:) = abs(mean(DisplacementInfo.dX)); [~, frameMax] = max(temp(1,:));
temp(2,:) = abs(mean(DisplacementInfo.dY)); [~, frameMax(end+1)] = max(temp(2,:));
temp(3,:) = abs(mean(DisplacementInfo.dZ)); [~, frameMax(end+1)] = max(temp(3,:));

frameRef = round(mean(frameMax));
	frameRef = 50;
	frameRef = 54; % 65;
	for t = frameStart:DENSEInfo.Number
		VectorsToVTK(Nodes, [], '', [DisplacementInfo.dX(:,t)-DisplacementInfo.dX(:,frameRef),DisplacementInfo.dY(:,t)-DisplacementInfo.dY(:,frameRef),DisplacementInfo.dZ(:,t)-DisplacementInfo.dZ(:,frameRef)], 'DisplacementVectors', fullfile(uipath,['DENSEdispREL',num2str(t)]))
	end	
	
	clear dispsVect;
	if isempty(frameRef)
		dispsVect(:,1,:) = DisplacementInfo.dX; dispsVect(:,2,:) = DisplacementInfo.dY; dispsVect(:,3,:) = DisplacementInfo.dZ;
	else
		dispsVect(:,1,:) = DisplacementInfo.dX - repmat(DisplacementInfo.dX(:,frameRef), 1, DENSEInfo.Number);
		dispsVect(:,2,:) = DisplacementInfo.dY - repmat(DisplacementInfo.dY(:,frameRef), 1, DENSEInfo.Number);
		dispsVect(:,3,:) = DisplacementInfo.dZ - repmat(DisplacementInfo.dZ(:,frameRef), 1, DENSEInfo.Number);
	end
	function
		disps = reshape(sqrt(sum(dispsVect.^2,2)), [], DENSEInfo.Number);
		tmp = abs(mean(disps)); [~, frameMin] = min(tmp);
		if frameMin ~= frameRef; msgbox(''); end
		[~, frameMax] = max(tmp)
temp(1,:) = abs(mean(dispsVect(:,1,:))); [~, frameMax(end+1)] = max(temp(1,:));
temp(2,:) = abs(mean(dispsVect(:,2,:))); [~, frameMax(end+1)] = max(temp(2,:));
temp(3,:) = abs(mean(dispsVect(:,3,:))); [~, frameMax(end+1)] = max(temp(3,:));
		
		for t = frameStart:DENSEInfo.Number
			Nodes = NodesRef + [DisplacementInfo.dX(:,t) DisplacementInfo.dY(:,t) DisplacementInfo.dZ(:,t)];
% Nodes = ROIInfo.Contour{t}; Nodes(:,3) = deal(0);
			% 2d to 3d:
			if flagCINE; Nodes = tformfwd(SequenceInfo(1).tform, Nodes); end
			VectorsToVTK(Nodes, disps(:,t), 'Displacements', dispsVect(:,:,t), 'DisplacementVectors', fullfile(uipath,['DENSEdispREL',num2str(t)]))
		end
	end
	
	%%% CMR images:
	col = double(SequenceInfo(1, 1).Columns); row = double(SequenceInfo(1, 1).Rows);
	% col = size(ImageInfo.Mag(:,:,1),1); row = size(ImageInfo.Mag(:,:,1),2);
	locs = [];
%% Test Allignment of Segmetation with CMR Mag Images:
% figure; imagesc(ImageInfo.Mag(:,:,1)); hold on; plot(DisplacementInfo.X, DisplacementInfo.Y, 'o'); 
	clear tmp; XX=1:col; YY=1:row;
	%% rotate 90 degree:
	for ii = XX % YY
		tmp(:,2) = (YY)'; tmp(:,1) = deal(ii);
		locs = [locs; tmp];
	end
	%{ 
	if DENSEInfo.NegFlag(2); XX=col:-1:1; else; XX=1:col; end
	if DENSEInfo.NegFlag(1); YY=row:-1:1; else; YY=1:row; end
	for ii = YY
		% locs = [locs; [repmat(ii,col,1), (1:col)']];
		tmp = (XX)'; tmp(:,2) = deal(ii);
		locs = [locs; tmp];
		% arrayfun(@(ii)[repmat(ii,col,1), (1:col)'], 1:row, 'UniformOutput', false)
	end
	 %}
	locs(:,3) = deal(0);
	
	%% DENSE magnitude images:
	for t = frameStart:DENSEInfo.Number
		mag = ImageInfo.Mag(:,:,t);
		VectorsToVTK(locs, mag(:), 'DENSEmag', [], '', fullfile(uipath,['DENSEmag',num2str(t)]))
	end

	%% cine SSFP images:
	% RCS to PCS:
	locs = tformfwd(SequenceInfo(1).tform, locs);
	
frameStart=24; uipath=pwd;
% SequenceInfo(1).ImagePositionPatient
	DirPx = SequenceInfo(1).tform.tdata.T(3,1:3);
	% DirPx = SequenceInfo(1).ImageOrientationPatient(4:end); 
	% DirPx = SequenceInfo(1).ImageOrientationPatient(1:3); 
	DirPx = cross(SequenceInfo(1).ImageOrientationPatient(1:3),SequenceInfo(1).ImageOrientationPatient(4:end));
	DirPx = reshape(DirPx, [], 3); DirPx = repmat(DirPx,col*row,1);
	for t = frameStart:DENSEInfo.Number
		mag = ImageInfo.Mag(:,:,t);
		VectorsToVTK(locs, mag(:), 'DENSEmag', DirPx, 'PixelDirection', fullfile(uipath,['cine',num2str(t)]))
	end
tformfwd(SequenceInfo(1).tform, [1.328125 1.328125 0])
1.9922 0.6641
	
	%%% Nodes Path:
	Nodes = [DisplacementInfo.X DisplacementInfo.Y];
	Nodes(:,3) = deal(0);
	VectorsToVTK(Nodes, [], '', [], '', 'DENSEnodes36')
	for t = frameStart:DENSEInfo.Number
		VectorsToVTK(Nodes + [DisplacementInfo.dX(:,t) DisplacementInfo.dY(:,t) DisplacementInfo.dZ(:,t)], [], '', [], '', fullfile(uipath,['DENSEnodes',num2str(t)]))
	end

end

VectorsToVTK(Nodes, [], '', [], '', 'DENSEnodes0)
for t = 1:33
	VectorsToVTK(Nodes + [DisplacementInfo.dX(:,t+36) DisplacementInfo.dY(:,t+36) DisplacementInfo.dZ(:,t+36)], [], '', [], '', fullfile(uipath,['DENSEnodes',num2str(t)]))
end

for t = 0:32
	VectorsToVTK(Nodes, [], '', [], '', fullfile(uipath,['DENSEnodes',num2str(t)]))
	Nodes = Nodes + [DisplacementInfo.dX(:,t+frameStart) DisplacementInfo.dY(:,t+frameStart) DisplacementInfo.dZ(:,t+frameStart)];
end

for t = frameStart:DENSEInfo.Number
	VectorsToVTK(Nodes, [], '', [], '', fullfile(uipath,['DENSEnodes',num2str(t-1)]))
	Nodes = Nodes + [DisplacementInfo.dX(:,t) DisplacementInfo.dY(:,t) DisplacementInfo.dZ(:,t)];
end

scatter3(Nodes(:,1),Nodes(:,2),Nodes(:,3),10,'r','x');

%%% mapping contours: voxel to patient coordinate system
            for k = 1:numel(ROIInfo.RestingContour)
                tmp = roi.RestingContour{k};
                tmp(:,3) = 0;
                RestingContour3D{k} = tformfwd(seq.tform, tmp);
            end

            for k = 1:numel(ROIInfo.Contour)
                tmp = roi.Contour{k};
                tmp(:,3) = 0;
                Contour3D{k} = tformfwd(seq.tform, tmp);
            end
