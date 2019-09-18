function result4SingleSubj(varargin)
% Last Modified: 1:06 PM Thursday, October 15, 2015
% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)

import plugins.DENSE3D_Plugin_4CrescentOrgan.*

%%%	Goto the directory in Matlab where Inputs are
%set current Dir.
files = DENSE2DobjectFileName()
button = questdlg('Do you happy with the auto-selected files in printout?','Warning','Yes');
if isempty(files) || ~strcmpi(button,'Yes')
	clear files folders;
	ii = 0;
	while true
		[uifile, uipath, uipopup] = uigetfile({'*.mat', 'Select DENSE2D Outputs (*.mat)'},'Open',pwd);
		if ~uipopup %isequal(uifile,0) || isequal(uipath,0)
			break
		else
			ii = ii + 1;
			files{ii} = fullfile(uipath,uifile);
			% tmp = fileparts(uipath);
			% [~,folders{ii}] = fileparts(tmp);
		end
	end
	% msgbox('Please place all .mat generated for each slice with DENSE2D in the current directory of Matlab','Error'); return
end

%% Standardize the Data loaded from DENSE2Dobject:
% fields = {'ImageInfo','AnalysisInfo','SequenceInfo','TransmuralStrainInfo','ROIInfo'};
fields = {'AnalysisInfo','SequenceInfo','ROIInfo'};
data = cellfun(@(x)load(x, fields{:}, '-mat'), files);
nFiles = numel(data);
% nFiles = size(files,2);
% recognize Input Files type: SA or LA?
indSA = []; indLA = [];
for ii = 1:nFiles
	if ~isempty(strfind(data(ii).SequenceInfo(1,1).DENSEanalysisName, 'SA'))
		indSA = [indSA, ii];
	elseif ~isempty(strfind(data(ii).SequenceInfo(1,1).DENSEanalysisName, 'LA'))
		indLA = [indLA, ii];
	else
		error(sprintf('%s:invalidInput',mfilename),'%s',...
			'Invalid Input #',sprintf('%d',ii),': ',sprintf('%s',data(ii).SequenceInfo(1,1).DENSEanalysisName.Orientation),'. Cannot recognize the type for SA or LA.');
		return
	end
end
nSlicesSA = numel(indSA);

locs = cellfun(@(x)x(1).SliceLocation, {data.SequenceInfo});
locs = locs(indSA);
% SA(top:base -> bottom:apex) -> LA
[~, indSA] = sort(locs, 'descend');
data = horzcat(data(indSA),data(indLA));


%% Load DENSE3Dobject:
uipath = pwd;	
[uifile,uipath] = uigetfile(...
	{'*.mat','Open a DENSE3D workspace(*.mat)'},...
	'Select MAT file',uipath);
if isequal(uipath,0) || isequal(uifile,0)
	while true
		tmp = input('1-by-3 DENSE3D Inputs of Positive Integers(Meshing Parameters Matrix)--"[#ofCircElemt,#ofRadialElemt,#ofLongiElemt]"[press enter for "[36,3,12]" by default]:');
		if isempty(tmp)
			tmp = [36,3,12];
			% Regions = 36; Layers = 3; Interp = 2;
		end
		% #ofLongiElemt = [(nSlicesSA-1)*Interp+nSlicesSA] - 1
		Regions = tmp(1); Layers = tmp(2); Interp = (tmp(3)+1-nSlicesSA)/(nSlicesSA-1);
		if Interp >= 0 || mod(Interp,1)==0
			Interp
			break
		end
	end
	dbstop if error
	[~,tmp] = fileparts(pwd);	
	filename = strcat('DENSE3Dobject','Mesh',num2str(Regions),num2str(Layers),num2str(Interp),'_',tmp);
	DENSE3Dobject = DENSE3D(files(indSA),'Layers',Layers,'Interp',Interp,'Regions',Regions,'Gif',filename)%,'Output','Type','Linear','Constant',0,'Smooth',0)
	save(filename,'DENSE3Dobject');
	% save(fullfile(pwd,folders{ii},['DENSE3Dobject',Insertion,SegDistro,num2str(Regions),num2str(Layers),num2str(Interp)),'DENSE3Dobject']);
else
	filename = fullfile(uipath,uifile);	
	% load file
	load(filename,'-mat');
	if ~isprop(DENSE3Dobject,'Mesh') && ~isprop(DENSE3Dobject,'Data')
		msgbox('"DENSE3Dobject.mat" file does not contain the necessary information');
		return
	end
	[~,filename] = fileparts(uifile);
	Interp = str2num(filename(end));
end


%% Define the frames of interest:
nFrames = size(data(1).ROIInfo.Contour,1)+1;
frame1 = 1;
while frame1 < 2 || mod(frame1,1)~=0 % ~isinteger(frame1)
	frame1 = input(strcat(['Input a No. of the initial frame where cavity shown up (frame#1 ~ frame#',num2str(nFrames),') [With RestingContour, frame#4(actually 3) by default when enter is pressed]:  ']));
	if isempty(frame1)
		frame1 = 4;
	end
end

%%% curve of Cavity Volume: cone approximation
flag_vol = questdlg('Do you wish to create a Cavity Volume curve using Cone Approximation?','Warning','Yes');
frameI = 0;
if strcmpi(flag_vol,'Yes')
	if sum(indLA)
		%% Select Frames of Interest:	
		while frameI < 1 || frameI > frameJ || frameJ < frame1 || mod(frameI,1)~=0 || mod(frameJ,1)~=0
			frameI = input(strcat(['Input a No. of the initial frame of interest (frame#1 ~ frame#',num2str(nFrames),') [frame#1(RestingContour) by default when enter is pressed]:  ']));
			if isempty(frameI)
				frameI = 1;
				% msgbox 'An input is needed';
				% return;
			end
			frameJ = input(strcat(['Input a No. of the last frame of interest (frame#',num2str(frameI),' ~ frame#',num2str(nFrames),') [With RestingContour, ES-frame#11(actually 10) by default when enter is pressed]:  ']));
			if isempty(frameJ)
				frameJ = 11;
				% msgbox 'An input is needed';
				% return;
			end
		end
	else	
		%% Load LA contours from .dns & auto-assign Frames of Interest:	
		tmp = fullfile(char(java.lang.System.getProperty('user.home')),'Dropbox','Analysis_Human');%Directory where DNS stored by default
		[uifile,uipath] = uigetfile(...
			{'*.dns','Open a DENSE2D workspace(*.dns)'},...
			'Select DNS file',tmp);
		if isequal(uipath,0) || isequal(uifile,0)
			% Load by dafault:
			[~,filename] = fileparts(pwd);	
			filename = fullfile(tmp,strcat(filename,'.dns'));
			load(filename,'-mat');
		else
			filename = fullfile(uipath,uifile);	
			% load file
			load(filename,'-mat');
		end
		indLA = strcmpi({roi.Type},'LA');
		ind = find(indLA);
		indFrames = 1;
		for ii = 1:sum(indLA)
			indFrames = indFrames.*cellfun(@(x)~isempty(x), {roi(ind(ii)).Position{:,1}});
			Frames = find(indFrames);
			% indFrames = find(cellfun(@(x)~isempty(x), {roi(ind(ii)).Position{:,1}}));
			for jj = 1: numel(Frames)
				samplePt = size(roi(ind(ii)).Position{Frames(jj),1},1);			
				sampleVal = roi(ind(ii)).Position{Frames(jj),1}(:,1);
				data(nFiles+ii).ROIInfo.Contour{Frames(jj),1}(:,1) = spline(1:samplePt, sampleVal, 1:samplePt/100:samplePt);
				sampleVal = roi(ind(ii)).Position{Frames(jj),1}(:,2);
				data(nFiles+ii).ROIInfo.Contour{Frames(jj),1}(:,2) = spline(1:samplePt, sampleVal, 1:samplePt/100:samplePt);
				
				samplePt = size(roi(ind(ii)).Position{Frames(jj),2},1);			
				sampleVal = roi(ind(ii)).Position{Frames(jj),2}(:,1);
				data(nFiles+ii).ROIInfo.Contour{Frames(jj),2}(:,1) = spline(1:samplePt, sampleVal, 1:samplePt/90:samplePt);
				sampleVal = roi(ind(ii)).Position{Frames(jj),2}(:,2);
				data(nFiles+ii).ROIInfo.Contour{Frames(jj),2}(:,2) = spline(1:samplePt, sampleVal, 1:samplePt/90:samplePt);
			end
			data(nFiles+ii).SequenceInfo(1, 1)= seq(roi(ind(ii)).SeqIndex(1));
		end
		
		if isempty(Frames)
			msgbox('CANNOT find a frame with LA contours!','Error');
			return		
		end
		
		% Validate Interpolation:
		figure
		plot(roi(ind(ii)).Position{Frames(jj),1}(:,1),roi(ind(ii)).Position{Frames(jj),1}(:,2),'LineWidth',1,'color','b','Marker','x','MarkerSize',5);  
		hold on;
		plot(data(nFiles+ii).ROIInfo.Contour{Frames(jj),1}(:,1),data(nFiles+ii).ROIInfo.Contour{Frames(jj),1}(:,2),'LineWidth',1,'color', 'r','Marker', '+', 'MarkerSize',5); 
		title('Validate Interpolation!');

		% Re-define constants:
		nFiles = numel(data);
		% With RestingContour: frame#4 is actually frame#3 
		Frames = Frames + 1;	
	end
else
	while frameI<1 || (frameI<frame1&&frameI~=1) || mod(frameI,1)~=0
		frameI = input(strcat(['Input a No. for the Frame of Interest (frame#1 ~ frame#',num2str(nFrames),') [With RestingContour, ES-frame#11(actually 10) by default when enter is pressed]:  ']));
		if isempty(frameI)
			frameI = 11;
			% msgbox 'An input is needed';
			% return;
		end
		% frameI = input(,'s');
		% frameI = sscanf(frame, '%i');
	end
	frameJ = frameI;
end

if exist('frameJ','var')
	if frameI < frame1
		Frames = [1, frame1 : frameJ];
	else
		Frames = [frameI : frameJ];
	end
end

%% Manually Fix dysfunction of SA & LA:
%% Continue to calculate the rest of frames for volume curve:
% Frames = [12 : frameJ];

% Compute Transformation Matrix: make sure the loaded .mat file is for SA
IOP = data(1).SequenceInfo(1, 1).ImageOrientationPatient;
%norm(IOP(1:3))
%norm(IOP(4:6))
T = [IOP(1:3),IOP(4:6),cross(IOP(1:3),IOP(4:6))];
% T*T'
sliceThickness = data(1).SequenceInfo(1, 1).SliceThickness;
for frame = Frames
	for slice=1:nFiles
		%load voxel coordinates
		IPP = data(slice).SequenceInfo(1, 1).ImagePositionPatient;
		IOP = data(slice).SequenceInfo(1, 1).ImageOrientationPatient';
		PS = data(slice).SequenceInfo(1, 1).PixelSpacing';
		if frame == 1
			epi = data(slice).ROIInfo.RestingContour{1, 1};
			endo = data(slice).ROIInfo.RestingContour{1, 2};		
		else
			epi = data(slice).ROIInfo.Contour{frame-1, 1};
			endo = data(slice).ROIInfo.Contour{frame-1, 2};
		end

		%transfer into world coordinates
		[X,Y,Z] = im2world(epi(:,1), epi(:,2), IPP, IOP, PS);
		world3D{frame}{slice,1} = [X(:), Y(:), Z(:)];	
		[X,Y,Z] = im2world(endo(:,1), endo(:,2), IPP, IOP, PS);
		world3D{frame}{slice,2} = [X(:), Y(:), Z(:)];

		%transfer into imaging coordinates
		for nPoint=1:size(world3D{frame}{slice,1},1)
			imaging3D{frame}{slice,1}(nPoint,:) = (inv(T) * world3D{frame}{slice,1}(nPoint,:)')';
		end
		for nPoint=1:size(world3D{frame}{slice,2},1)
			imaging3D{frame}{slice,2}(nPoint,:) = (inv(T) * world3D{frame}{slice,2}(nPoint,:)')';
		end
		
		[~,tmp] = size(data(slice).ROIInfo.Contour);
		if tmp > 2
			if frame == 1
				endoRV = data(slice).ROIInfo.RestingContour{1, 3};
			else
				endoRV = data(slice).ROIInfo.Contour{frame-1, 3};
			end		
			[X,Y,Z] = im2world(endoRV(:,1), endoRV(:,2), IPP, IOP, PS);
			world3D{frame}{slice,3} = [X(:), Y(:), Z(:)];
			for nPoint=1:size(world3D{frame}{slice,3},1)
				imaging3D{frame}{slice,3}(nPoint,:) = (inv(T) * world3D{frame}{slice,3}(nPoint,:)')';
			end
		end
	end


	%%% Aligning contours with coordinates created by DENSE3D program:
	%{ 
	IOP = [1 0 0 0 1 0];
	for slice = 1:nSlicesSA
		IPP = [0 0 DENSE3Dobject.Data(slice).SequenceInfo(1, 1).SliceLocation];
		PS = DENSE3Dobject.Data(slice).SequenceInfo(1, 1).PixelSpacing;
		
		endo = DENSE3Dobject.Data(slice).ROIInfo.Contour{frame, 1};
		epi = DENSE3Dobject.Data(slice).ROIInfo.Contour{frame, 1};

		[X,Y,Z] = im2world(endo(:,1), endo(:,2), IPP, IOP, PS);
		dense{frame}{slice,2} = [X(:), Y(:), Z(:)];
		[X,Y,Z] = im2world(epi(:,1), epi(:,2), IPP, IOP, PS);
		dense{frame}{slice,1} = [X(:), Y(:), Z(:)];
	end
	%%% Get translation{frame} here
	translation{frame} = [0 0 0];
	for ii = 1:nSlicesSA
		for jj = 1:2
			translation{frame} = mean(dense{frame}{ii,jj})-mean(imaging3D{frame}{ii,jj}) + translation{frame};
		end
	end
	translation{frame} = (translation{frame} / (ii*jj))

	%%% Plot
	% figure;
	% hold on;
	for slice=1:nFiles
		%translate imaging coordinates into dense coordinates
		for nPoint=1:size(imaging3D{frame}{slice,1},1)
			imaging3D{frame}{slice,1}(nPoint,:) = imaging3D{frame}{slice,1}(nPoint,:) + translation{frame};
		end
		for nPoint=1:size(imaging3D{frame}{slice,2},1)
			imaging3D{frame}{slice,2}(nPoint,:) = imaging3D{frame}{slice,2}(nPoint,:) + translation{frame};
		end

		%% debug:
		% if slice <= nSlicesSA
			% LineColor = 'b';
		% else
			% LineColor = 'r';
		% end		
		% plot3(imaging3D{frame}{slice,1}(:,1),imaging3D{frame}{slice,1}(:,2),imaging3D{frame}{slice,1}(:,3),'LineWidth',1.5,'color', LineColor);
		% plot3(imaging3D{frame}{slice,2}(:,1),imaging3D{frame}{slice,2}(:,2),imaging3D{frame}{slice,2}(:,3),'LineWidth',1.5,'color', LineColor);
		%% 'LineStyle' = '-'; 'Marker' = 'none';
	end
	% vertices = DENSE3Dobject.Mesh.strains.vtrj(:,:,frame);
	% scatter3(vertices(:,1),vertices(:,2),vertices(:,3),3,'k','*');
	% %'LineStyle' = 'none'; 'Marker' = '*';
	% hold off;
	% %axes('XAxisLocation','origin');
	% axis equal;
	% xlabel('x');
	% ylabel('y');
	% zlabel('z');
	% grid on;

	% translation{frame} @ cell array:
	%{ 
	for k = 1 : 2
		slice = nFiles + k;
		for ii = 1 : 2
			tmp = numel(imaging3D{frame}{slice, ii});
			if tmp == 0;continue;end
			for jj = 1 : tmp
				imaging3D{frame}{slice, ii}{jj} = imaging3D{frame}{slice, ii}{jj} + repmat(translation{frame},size(imaging3D{frame}{slice, ii}{jj},1),1);
			end
		end
	end
	 %}

	button = questdlg('Do you wish to continue aligning contours with meshing of strain calculation?','Warning','No');
	if strcmpi(button,'No') || strcmpi(button,'Cancel')
		button = questdlg('Are you sure to quit?','Warning');
		if strcmpi(button,'Yes')
			return
		end
	end
	%}



	%%% Aligning contours with meshing of strain calculation:
	%% Assumption: last slices of contour and mesh
	% 2 elemt bet each 2 slices
	if ~isfield(DENSE3Dobject.Mesh,'Interp')
		DENSE3Dobject.Mesh.Interp = 1;
	end
	nSlice_vertices = nSlicesSA + (nSlicesSA-1)*DENSE3Dobject.Mesh.Interp;
	interval = size(DENSE3Dobject.Mesh.strains.vtrj(:,:,frame),1)/nSlice_vertices;
	%%% Approximated value of Atrioventricular plane displacement
	AVPD{frame} = imaging3D{frame}{1,1}(1,3) - mean(DENSE3Dobject.Mesh.strains.vtrj(1:interval,3,frame));

			
	% center of mesh
	center_mesh{frame} = mean(DENSE3Dobject.Mesh.strains.vtrj(:,:,frame));
	tmp = mean(DENSE3Dobject.Mesh.strains.vtrj(end-interval+1:end,:,frame));
	center_mesh{frame}(3) = tmp(3);

	% center of contour
	for slice = 1:nSlicesSA
		center_contour{frame}(slice,:) = (mean(imaging3D{frame}{slice,1})+mean(imaging3D{frame}{slice,2}))/2;
		% mesh3D {slice,1} = DENSE3Dobject.Mesh.strains.vtrj(interval*(nSlice_vertices-1)+1:interval*nSlice_vertices,:,frame);
	end
	center_contour{frame} = [mean(center_contour{frame}(:,1:2)),center_contour{frame}(nSlicesSA,3)];

	%% translation{frame}
	translation{frame} = center_mesh{frame} - center_contour{frame};

	%%% Plot
	% figure;
	% hold on;
	for slice=1:nFiles
		%translate imaging coordinates into dense coordinates
		% for nPoint=1:size(imaging3D{frame}{slice,1},1)
			% imaging3D{frame}{slice,1}(nPoint,:) = imaging3D{frame}{slice,1}(nPoint,:) + translation{frame};
		% end
		imaging3D{frame}{slice,1} = imaging3D{frame}{slice,1} + repmat(translation{frame},size(imaging3D{frame}{slice,1},1),1);
		% for nPoint=1:size(imaging3D{frame}{slice,2},1)
			% imaging3D{frame}{slice,2}(nPoint,:) = imaging3D{frame}{slice,2}(nPoint,:) + translation{frame};
		% end
		imaging3D{frame}{slice,2} = imaging3D{frame}{slice,2} + repmat(translation{frame},size(imaging3D{frame}{slice,2},1),1);

		%% debug:
		% if slice <= nSlicesSA
			% LineColor = 'b';
		% else
			% LineColor = 'r';
		% end		
		% plot3(imaging3D{frame}{slice,1}(:,1),imaging3D{frame}{slice,1}(:,2),imaging3D{frame}{slice,1}(:,3),'LineWidth',1.5,'color', LineColor);
		% plot3(imaging3D{frame}{slice,2}(:,1),imaging3D{frame}{slice,2}(:,2),imaging3D{frame}{slice,2}(:,3),'LineWidth',1.5,'color', LineColor);
	end
	% vertices = DENSE3Dobject.Mesh.strains.vtrj(:,:,frame);
	% scatter3(vertices(:,1),vertices(:,2),vertices(:,3),3,'k','*');
	%% Only plot one slice of vertices:
	% slice = 1;
	% scatter3(mesh3D{slice,1}(:,1),mesh3D{slice,1}(:,2),mesh3D{slice,1}(:,3),3,'k','*');
	% hold off;
	% axis equal;
	% xlabel('x');
	% ylabel('y');
	% zlabel('z');
	% grid on;
end


%% Plot to verify Transformation:
for frame = 1 % Frames
	figure;
	hold on;
	for slice=1:nFiles
		% debug:
		if slice <= nSlicesSA
			LineColor = 'b';
		else
			LineColor = 'r';
		end		
		plot3(imaging3D{frame}{slice,1}(:,1),imaging3D{frame}{slice,1}(:,2),imaging3D{frame}{slice,1}(:,3),'LineWidth',1.5,'color', LineColor);
		plot3(imaging3D{frame}{slice,2}(:,1),imaging3D{frame}{slice,2}(:,2),imaging3D{frame}{slice,2}(:,3),'LineWidth',1.5,'color', LineColor);
		[~,tmp] = size(imaging3D{frame});
		if tmp > 2
			plot3(imaging3D{frame}{slice,3}(:,1),imaging3D{frame}{slice,3}(:,2),imaging3D{frame}{slice,3}(:,3),'LineWidth',1.5,'color', 'k');
		end
	end
	% vertices = DENSE3Dobject.Mesh.strains.vtrj(:,:,frame);
	% scatter3(vertices(:,1),vertices(:,2),vertices(:,3),3,'k','*');
	% for slice=1:4 %GoodSlices{frame}
		% plot3(imaging3D{frame}{nFiles+2,1}{slice}(:,1),imaging3D{frame}{nFiles+2,1}{slice}(:,2),imaging3D{frame}{nFiles+2,1}{slice}(:,3),'LineWidth',1.5,'color', 'k');
		% plot3(imaging3D{frame}{nFiles+2,2}{slice}(:,1),imaging3D{frame}{nFiles+2,2}{slice}(:,2),imaging3D{frame}{nFiles+2,2}{slice}(:,3),'LineWidth',1.5,'color', 'k');
	% end
	% for slice=1:numel(imaging3D{frame}{nFiles+3,1})
		% plot3(imaging3D{frame}{nFiles+3,1}{slice}(:,1),imaging3D{frame}{nFiles+3,1}{slice}(:,2),imaging3D{frame}{nFiles+3,1}{slice}(:,3),'LineWidth',1.5,'color', 'k');
		% plot3(imaging3D{frame}{nFiles+3,2}{slice}(:,1),imaging3D{frame}{nFiles+3,2}{slice}(:,2),imaging3D{frame}{nFiles+3,2}{slice}(:,3),'LineWidth',1.5,'color', 'k');
	% end
	hold off;
	axis equal;
	xlabel('x');
	ylabel('y');
	zlabel('z');
	grid on;
	title('Validate Transformation!');
	% hgexport(gcf,'AllCoutours',hgexport('factorystyle'), 'Format', 'jpeg');
end

%%% Create Apex
if nSlicesSA == nFiles
	msgbox('No .mat files for LA slice detected!','Error');	
	return
end
flag_shape = input('The flag_shape of section(z-axis as normal) for Apex [press "1" for 2-D circle(no fitting); "2" for 2-D ellipse; "3" for 2-D ellipse(exact fitting); "4" for 3-D circle; any key for 3-D ellipse]: ');
if isempty(flag_shape) || (flag_shape~=1 && flag_shape~=2 && flag_shape~=3 && flag_shape~=4)
	flag_shape = 5;
end
nPointsInterp = input('The No. of points on a ellipse slice of apex(12 by default when enter is pressed): ');
if isempty(nPointsInterp)
    nPointsInterp = 12;
end

%epi-
if strcmpi(flag_vol,'Yes')
	button = 'Yes';
else
	button = questdlg('Do you wish to create Apex for Epicardium?','Warning','Yes');
end
% Frames1=?:frameJ
if strcmpi(button,'Yes')
	col = 1; %1st column of imaging3D{frame}
	for frame = Frames %Frames1
		flag_plot = true;
		% flag_vol = false;
		[imaging3D{frame},output{frame,col},~,nSlicesLA{frame,col},FitParams{frame,col},merit{frame,col},GoodSlices{frame,col}]=buildApex(imaging3D{frame},col,frame,nFiles,sliceThickness,nSlicesSA, nPointsInterp, flag_shape, flag_plot, 'No');
		% flag_vol='No': only endocardium is needed to compute cavity volumn
	end
end

%endo-
if strcmpi(flag_vol,'Yes')
	button = 'Yes';
else
	button = questdlg('Do you wish to create Apex for Endocardium?','Warning','Yes');
end
% Frames1=?:frameJ
if strcmpi(button,'Yes')
	col = 2; %2nd column of imaging3D{frame}
	for frame = Frames %Frames1
		flag_plot = true;
		% flag_vol = false;
		[imaging3D{frame},output{frame,col},volume{frame},nSlicesLA{frame,col},FitParams{frame,col},merit{frame,col},GoodSlices{frame,col}]=buildApex(imaging3D{frame},col,frame,nFiles,sliceThickness,nSlicesSA, nPointsInterp, flag_shape, flag_plot, flag_vol);
	end
end



%%% Output:
button = questdlg('Do you wanna save the result?','Yes');
if strcmpi(button,'Yes')
	uipath = pwd;
	if exist('frameJ','var')
		if frameI == frameJ 
			defaultpath = fullfile(uipath,['Volume_Frame',num2str(frameI)])
		else
			defaultpath = fullfile(uipath,['Volume_Frame',num2str(frameI),'-',num2str(frameJ)])
		end
	else
		defaultpath = fullfile(uipath,['Volume_Frame',sprintf(['_','%d'], Frames)])
	end
	while true
		[uifile, uipath, uipopup] = uiputfile(...
			{'*.mat', 'Save a DENSE3D workspace (*.mat)'},...
			'Save',defaultpath);
		if uipopup
			break
		else
			button = questdlg('Do you wanna quit? You might lose data!','Warning','No');
			if strcmpi(button,'Yes')
				return
			end
		end
	end

	%% Append the correct file extension if needed
	filename = fullfile(uipath,uifile);
	[~,~,ext] = fileparts(filename);
	if ~strcmp(ext, '.mat')
		filename = strcat(filename, '.mat');
	end
			
	% save(filename,'imaging3D','world3D','dense','output');
	save(filename,'imaging3D');
	if exist('output', 'var')
		save(filename,'output','-append');
	end
	if exist('volume', 'var')
		% interpolation for cavity-no-shown frames:
		if ~isempty(volume{1}) && frame1>2 && ~isempty(volume{frame1})
			tmp = interp1([1,frame1], [volume{1},volume{frame1}], 1:frame1);
			for ii = 2:numel(tmp)-1
				volume{ii} = tmp(ii);
			end
		end
		save(filename,'volume','-append');
	end
	
	% Table of the whole volume curve:
	if exist('frameJ','var')
		xlswrite(strcat('Volume_Frame',num2str(frameI),'-',num2str(frameJ)),volume');%[num2cell(Frames);volume]);
	else
		xlswrite(strcat('Volume_Frame',sprintf(['_','%d'], Frames)),volume');%[num2cell(Frames);volume]);
	end	
	plot(Frames,cell2mat(volume),'LineWidth',1,'color','b','Marker','x','MarkerSize',5); 
	title({'Volume vs Frame','Frame #1 is the resting contour','!!! Actual Frame # = Frame # in the figure - 1'});
	[ESV,ES] = min([volume{:}]);
	ES = Frames(ES)
	tmp = input('The Early-Diastolic frame #:');
	if ~isempty(tmp)
		ES = tmp;
	end
	
	[~,filename]=fileparts(pwd);
	ind = round(nSlicesSA/2);
	if mod(nSlicesSA,2) == 0
		if mod(Interp,2) == 0
			MidSA = ind + Interp/2*(1/(Interp+1));
		else
			msgbox('When #ofSlices is even and ROI is mid-ventricle: DENSE3D Input "Interp" must be even!','Error');
			return
		end
	else
		MidSA = ind;
		% MidSA = [ind-1/(Interp+1), ind];
	end
	ind = abs(MidSA-DENSE3Dobject.Mesh.absind)<1e-2;
	ind2 = (2 == DENSE3Dobject.Mesh.layerid);
	index = ind&ind2;
	MidSAelemt=DENSE3Dobject.Mesh.strains.ptrj(index,:,ES);
	
	% '1st point on Apical SA Slice:',num2str(imaging3D{1,ES}{1,1}(1,:))];['1st point on Basal SA Slice:',num2str(imaging3D{1,ES}{nSlicesSA,1}(1,:))];
	Report={'This is an auto generated report:';['Dataset',filename];'!!! Actual Frame # = Reported Frame #  - 1 (Frame #1 is the resting contour)';['ES frame #',num2str(ES)];['ESV=',num2str(ESV)];['ED frame #',num2str(Frames(end))];['EDV=',num2str(volume{end})];['1st point on Basal SA Slice:',num2str(DENSE3Dobject.Mesh.strains.ptrj(1,:,ES))];['Last point on Apical SA Slice:',num2str(DENSE3Dobject.Mesh.strains.ptrj(end,:,ES))];'#ofElement bet. Basal & Api. SA in Long. dir. = ? (This No. is used for DENSE3D Mesh Adjustment)';['Mid-ventricular Anterior RV Insert Point:',num2str(MidSAelemt(1,:))];['Mid-ventricular Inferior RV Insert Point:',num2str(MidSAelemt(data(round(MidSA)).AnalysisInfo.SegDistribution{2},:))]}
	
end

end