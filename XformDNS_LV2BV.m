function XformDNS_LV2BV(button,config,self)
% XformDNS_LV2BV  - transform .dns inFiles of SA and LA
%
%
%	transform .dns inFiles from Version 0.4 2009.06 to Version: 0.5
%
%
% USAGE:
%   XformDNS_LV2BV();
%
%
% See also:
%   
%
% Last Modified: 17:38 June 14, 2017
% Copyright (c) Zhanqiu Liu (lafeir.lew@gmail.com)

	% dbstop if error
	
	%% DIRECTORY LOAD: default inputs
	homepath = char(java.lang.System.getProperty('user.home'));
	pathDNS = fullfile(homepath,'Dropbox\Analysis')%Directory where DNS stored by default
	pathMAT = fullfile(homepath,'Dropbox\MAT')%Directory where DENSE2D outputs(post-analyzed .MAT) stored by default

	if ~exist('button','var')
		button = input('Press any key for using a User Interface to select Datasets of interest.\nOtherwise leave it blank!\n','s');
	end
	
	if exist('config','var') && exist(config,'file')
		inFiles{1} = config;
		[~, fname{1}, ext] = fileparts(config);
		files{1} = [fname{1}, ext];
	else
		if isempty(button)
			tmp = dir(pathDNS); files = {tmp.name}; idx = strwcmpi(files,'*.dns'); files = files(idx);	
			tmp = arrayfun(@num2str, 1:numel(files), 'UniformOutput', false);
			tmp = strcat(tmp, '.', files);
			idx = input(strcat([sprintf('%s\n',tmp{:}),'Input a MATRIX of the numbers of Datasets of interest[use Space/Space/Semicolon(;) to seperate]:']));
			idx = reshape(idx,[],1); idx = round(abs(idx)); idx(idx==0) = [];
			% Default input datasets:
			if isempty(idx)
				idx = [9 19 8 18 5 15 6 16 7 17];
				files(idx)
			end
			files = files(idx);
			inFiles = fullfile(pathDNS,files);
			[~, fname] = cellfun(@fileparts,files,'UniformOutput',false);
		else
			clear inFiles files fname;
			ii = 0;
			while true
				[uifile, uipath, uipopup] = uigetfile({'*.dns','DENSEanalysis(v0.4) workspace(*.dns)'},'Select a to-be-transfered DNS file once a time',pathDNS);
				if uipopup %~isequal(uifile,0) || ~isequal(uipath,0)
					ii = ii + 1;
					inFiles{ii} = fullfile(uipath,uifile);
					files{ii} = uifile;
					[~, fname{ii}] = fileparts(files{ii});
					% if ii > 1 && strcmpi(inFiles{ii}, inFiles{ii-1})
						% break
					% end
				else
					break
				end
			end
			if ~exist('files','var'); return; end;
		end
	end
	nFiles = numel(files);

	if isempty(button)
		tmp = {'[Use the Data Cursor tool to confirm the directions of X-axis and Y-axis]','From X-axis to Y-axis: ','Press ENTER for ClockWise (rightward X-axis and downward Y-axis)','ANY INPUTS for CounterClockWise (rightward X-axis and upward Y-axis)','Your observation:'};%or input "CW"
		clcDir = input(strcat([sprintf('%s\n',tmp{:})]),'s');
		if isempty(clcDir)
			clcDir = -1;
			% dir_clc = 'CW';
		else
			clcDir = 1;		
			% dir_clc = 'CCW';		
		end		
	else
		clcDir = questdlg({'[Use the Data Cursor tool to confirm the directions of X-axis and Y-axis]','From X-axis to Y-axis: ','ClockWise? (rightward X-axis and downward Y-axis)'},'Clock Direction of Image Coordinate System');
		if strcmpi(clcDir,'Yes')
			clcDir = -1;		
		else
			clcDir = 1;		
		end
	end	
	
	if ~exist('self','var') || isempty(self.aside)
		self.aside = struct('ant',3,'pos',0);
	end
	%% Setting of distance of Endo-RV-insertPT away from Epi-RV-insertPT in SA view:
	while true
		if isempty(button)
			tmp = {'Input a INTEGER from -60 to 120:','+1 for setting Anterior-Endo-RV-insertPT 17 degree away from Anterior-Epi-RV-insertPT in the CCW direction in the Image Viewer','[Press ENTER for the default value]:'};%{['Input a INTEGER:','+1 for setting Anterior-Endo-RV-insertPT 17 degree away from Anterior-Epi-RV-insertPT in the ',clcDir,' direction'],'Sign Convention:','Input ANY STRINGS EXCEPT "yes"'}
			self.aside.ant = input(strcat([sprintf('%s\n',tmp{:})]));
			if isempty(self.aside.ant)
				self.aside.ant = 3;
			end
		else
			tmp = inputdlg('Input a INTEGER from -60 to 120: +1 for setting Anterior-Endo-RV-insertPT 17 degree away from Anterior-Epi-RV-insertPT in the CCW direction in the Image Viewer','Observe the outFiles and Adjust: Positive=CCW!',1,{num2str(self.aside.ant)});
			%% Allowed return will destory the value backup of self.aside
			if isempty(tmp)
				return
			else
				tmp = sscanf(sprintf('%s*', tmp{:}), '%f*');
				% self.aside.ant = str2num(self.aside.ant{:});%slow
			end
		end
		if mod(tmp,1)==0 && abs(tmp)<60 % self.aside.ant >= 0 || isinteger(self.aside.ant)
			self.aside.ant = tmp;
			break
		end		
	end
	while true
		if isempty(button)
			tmp = {'Input a INTEGER from -60 to 120:','+1 for setting Posterior-Endo-RV-insertPT 17 degree away from Posterior-Epi-RV-insertPT in the CCW direction','[Press ENTER for the default value]:'};
			self.aside.pos = input(strcat([sprintf('%s\n',tmp{:})]));
			if isempty(self.aside.pos)
				self.aside.pos = 0;
			end
		else
			tmp = inputdlg('Input a INTEGER from -60 to 120: +1 for setting Posterior-Endo-RV-insertPT 17 degree away from Posterior-Epi-RV-insertPT in the CCW direction','Observe the outFiles and Adjust: Positive=CCW!',1,{num2str(self.aside.pos)});
			if isempty(tmp)
				return
			else
				tmp = sscanf(sprintf('%s*', tmp{:}), '%f*');
			end
		end
		if mod(tmp,1)==0 && abs(tmp)<60
			self.aside.pos = tmp;
			break
		end		
	end
	
	if ~isempty(button)
		hwait = waitbartimer();
		cleanupObj = onCleanup(@(x)delete(hwait(isvalid(hwait))));
		hwait.String = 'Converting Contours...';
		hwait.WindowStyle = 'modal';
		hwait.AllowClose = false;
		hwait.Visible = 'on';
		hwait.start;
		drawnow
	end

	%% Setting of distance of Endo-RV-insertPT away from Epi-RV-insertPT in LA view:
	self.aside.LA = 1;
	
	% Define some constants:
	old = {'curve','line','SAFull', 'LAFull'};
	new = {'curved','poly','sadual','ladual'};
	% ROITYpe:
	% rois.CurvedRegion; rois.PolylineRegion; rois.LVShortAxis; rois.LVLongAxis; rois.ClosedContour; rois.OpenContour; plugins.dense3D_plugin.SADualContour; plugins.dense3D_plugin.LADualContour; 
	nText = numel(old);
	%% Setting of #ofPoints of SA RV-endo:
	% 12 pts for SA RV-endo by default(#1 & #7 are corner pts):
	nRVendoPt = 12;% got to be even
	RVendoCorner = repmat(false, [nRVendoPt,1]); RVendoCorner(1) = true; RVendoCorner(nRVendoPt/2+1) = true;
	
	for ii = 1 : nFiles
		%% test:
		% t1 = load(fullfile(homepath,'Dropbox','Analysis','20150628.dns'),'-mat');
		load(inFiles{ii},'-mat');
		if ~all(cellfun(@(x)exist(x,'var'),{'seq','img','dns','roi'}))
			msgbox('".dns" file does not contain the necessary information');
		end
		
		%load DENSE2D outputs(post-analyzed .MAT)
		patientName = regexpi(fname{ii}, '\_', 'split'); patientName = patientName{1};
		pathPostMAT = fullfile(pathMAT,patientName);
		if exist(pathPostMAT,'dir')
			MagIndex = reshape([dns.MagIndex],3,[]);MagIndex = MagIndex(1,:);
			tmp = dir(pathPostMAT);
			tmp = {tmp.name};
			%% DENSEanalysis(v0.5) cannot save variable "status"
			pathDNSori = fullfile(pathDNS,[patientName,'.dns']);
			if (~exist('status','var') || isempty(status.SOI)) && exist(pathDNSori,'file')
				try
					load(pathDNSori,'status','-mat');
				catch
					if isempty(button)
						warning('NO "status" information found!(variable "status" not in the workspace file .dns)');
					else
						msgbox('NO "status" information found!(variable "status" not in the workspace file .dns)');
					end
				end
			end
			if exist('status','var') && ~isempty(status.SOI)
				tmp2 = dns(status.SOI(1)).Name;
				tmp1 = find(strncmpi(tmp, tmp2, numel(tmp2)));	
				pathAnalysisMAT = fullfile(pathPostMAT,tmp{tmp1(1)});
			else
				tmp1 = strncmpi(tmp, 'auto.', 5);
				tmp2 = find(tmp1);
				pathAnalysisMAT = fullfile(pathPostMAT,tmp{tmp2(ceil(sum(tmp1)/2))});
				% msgbox('Variable "status" NOT exist!','s');
				% return
			end
			load(pathAnalysisMAT, 'AnalysisInfo', '-mat');
			if ~isnan(AnalysisInfo.SegDistribution{1, 2})
				SegPos = AnalysisInfo.SegDistribution{1, 2};
			else
				SegPos = 15
			end
		else
			msgbox(['CANNOT FOUND the directory where DENSE2D outputs stored for the input',files(ii),'A update is required for the Codes!']);
			keyboard; %continue				
			%{ 
			clckPOS = input('The Longitudinal Ventricular Location of Insertion Points of interest[press enter for "Mid", "1" for Base, "2" for "Apex", "3" for "Whole", "4" for "Average"; or Input a n-by-2 Matrix]: ');
			if isempty(clckPOS)
				clckPOS = 'Mid';
			end							
			SeptumRatio = 
			AnalysisInfo.PositionB =
			SegPos/60*360
			 %}
		end

		nROI = numel(roi);

		% Clear points of SA RV-endo:
		RVendo = zeros(nRVendoPt,2);
		
		idx = find(strcmpi({roi.Type},'SA'));
		% idx = strcmpi({roi.Type},'SA');
		if ~isempty(idx)
		for jj = idx
		% for jj = find(idx)
			roiName = [dns(MagIndex == roi(jj).SeqIndex(1)).Name,'_',roi(jj).Name];
			pathAnalysisMAT = fullfile(pathPostMAT,[roiName,'.mat']);
			if ~exist(pathAnalysisMAT,'file')
				[uifile, uipath] = uigetfile({'*.mat','DENSE2D(v0.4) outputs(post-analyzed *.MAT)'},['Select MAT file for SA ROI:',roiName,' of the input:',files{ii}],pathPostMAT);
				if isequal(uipath,0) || isequal(uifile,0)
					continue
				end
				pathAnalysisMAT = fullfile(uipath,uifile);					
			end
			load(pathAnalysisMAT, 'AnalysisInfo','DisplacementInfo', '-mat')

			idx1 = find(strcmpi({roi.Name},'Short Axis Dual Ventricle'));
			idx2 = cellfun(@(x)isequal(x,roi(jj).SeqIndex), {roi(idx1).SeqIndex});%,'UniformOutput', false)
			switch sum(idx2)
				case 1
					ROI = idx1(idx2);
					% warning(['Old auto-built ROI detected: Plz CHANGE ROI NAME ',roiName,' if you have modified this SA ROI for the input',files(ii)]);
				case 0
					nROI = nROI + 1;
					ROI = nROI;
				otherwise
					msgbox(['Duplicate converted ROIs detected:',num2str(sum(idx1)),' converted ROIs found for the SA ROI:',roiName,' of the input',files(ii)]);
					continue
			end
			
			roi(ROI).Name = 'Short Axis Dual Ventricle';
			roi(ROI).Type = 'sadual';
			roi(ROI).UID = dicomuid;
			roi(ROI).SeqIndex = roi(jj).SeqIndex;
			nFrames = size(roi(jj).Position,1);
			roi(ROI).IsClosed = repmat({true}, [nFrames,3]);
			roi(ROI).IsCurved = repmat({true}, [nFrames,3]);
			% roi(jj).IsCurved = repmat({true}, [nFrames,2]);
			roi(jj).IsCurved = cellfun(@logical,roi(jj).IsCurved,'UniformOutput',false);
			roi(ROI).IsCorner = repmat({false}, [nFrames,3]);
			% roi(ROI).IsCorner = repmat({false,false,RVendoCorner}, [nFrames,1]);
			% roi(jj).IsCorner = repmat({false}, [nFrames,2]);
			roi(jj).IsCorner = cellfun(@logical,roi(jj).IsCorner,'UniformOutput',false);
			
			roi(ROI).Position = roi(jj).Position;
			% roi(ROI).Position{:, 2} = roi(jj).Position{:, 2};
			idx1 = cellfun(@isempty, roi(jj).Position);
			if ~isequal(idx1(:,1),idx1(:,2))
				msgbox(['Missing LV-endo or LV-epi countours for the SA ROI:',roiName,' of the input',files(ii)]);			
			end
			idx1 = find(~idx1(:,1))'; % idx1(:,2) = [];
			
			if isempty(idx1)
				if ROI < numel(roi)
					msgbox(['No contours found in the SA ROI:',roiName,' of the input',files(ii),'But it has been converted before!']);
					return
				else
					roi(ROI) = [];
					nROI = nROI - 1;
				end				
				continue;
			else
				for kk = idx1
					roi(ROI).IsCorner{kk,3} = RVendoCorner;
					% roi(ROI).IsCorner{kk,3} = false;
					
					%% Project InsertPt:
					dx = fnvalmod(DisplacementInfo.spldx,[AnalysisInfo.PositionB([2 1]),kk]');
					dy = fnvalmod(DisplacementInfo.spldy,[AnalysisInfo.PositionB([2 1]),kk]');
					InsertPt = AnalysisInfo.PositionB + [dx,dy];
					% NO PIXEL# FOUND!:find(abs(DisplacementInfo.dX(:,kk)-dx) < 1e-04)
					% InsertPt = AnalysisInfo.PositionB+[DisplacementInfo.dX(:,kk),DisplacementInfo.dY(:,kk)];
					
					%% Interpolation
					% Need a closed circle:
					LVepi = [roi(jj).Position{kk, 1};roi(jj).Position{kk, 1}(1,:)];
					nSamplePt = size(LVepi,1);
					% AnalysisInfo.Nmodel = 60: not stored!
					LVepiInterp(:,1) = spline(1:nSamplePt, LVepi(:,1), 1:(nSamplePt-1)/59:nSamplePt);
					LVepiInterp(:,2) = spline(1:nSamplePt, LVepi(:,2), 1:(nSamplePt-1)/59:nSamplePt);
					%% find nearest point:
					tmp = bsxfun(@minus,LVepiInterp,InsertPt);
					[~,idx_AntInsPt] = min(sum(tmp.^2,2));
					
					%% Insert 2 RVinsertPts into LVepi
					%{
					% find nearest point: work for over-flat LVepi- ellipse
					if clcDir == -1
						LVepiInterp_ReOrd = LVepiInterp([idx_AntInsPt:-1:1,end:-1:idx_AntInsPt+1],:);
					else
						LVepiInterp_ReOrd = LVepiInterp([idx_AntInsPt:end,1:idx_AntInsPt-1],:);
					end
					for ll = [1,SegPos]
						% roi(jj).Position{kk,1}: NOT save the iteration result 
						LVepi = roi(ROI).Position{kk,1};
						tmp = bsxfun(@minus,LVepi,LVepiInterp_ReOrd(ll,:)); 
						[~,idx2] = sort(sum(tmp.^2,2),'ascend');
						if abs(idx2(1)-idx2(2)) > 1
							warning(['Insertion points NOT SAVED in the LVepi since random order of LVendo- points stored at the frame#',num2str(kk),'for the SA ROI:',roiName,' of the input',files{ii}]);%msgbox
							continue
						end
						idx3 = sort(idx2(1:2),'ascend');
						v1 = LVepi(idx2(2),:)-LVepi(idx2(1),:); v2 = LVepiInterp_ReOrd(ll,:)-LVepi(idx2(1),:);
						% angle bet 2 vectors = mod(-180/pi * angle, 360);
						if atan2(abs(det([v1;v2])),dot(v1,v2)) < pi/2
							roi(ROI).Position{kk,1} = [LVepi(1:idx3(1),:);LVepiInterp_ReOrd(ll,:);LVepi(idx3(2):end,:)];
						elseif idx2(1) < idx2(2)
							roi(ROI).Position{kk,1} = [LVepi(1:idx2(1)-1,:);LVepiInterp_ReOrd(ll,:);LVepi(idx2(1):end,:)];						
						else % idx2(1) > idx2(2)
							roi(ROI).Position{kk,1} = [LVepi(1:idx2(2),:);LVepiInterp_ReOrd(ll,:);LVepi(idx2(2)+1:end,:)];
						end
					end
					%}

					%% Insert 2 Pts near the 2 RVinsertPts into LVepi
					tmp = circshift(LVepiInterp,-(idx_AntInsPt-1+clcDir*SegPos));
					% tmp(1,:): LVepiInterp(idx_AntInsPt+clcDir*SegPos,:)
					for ll = {InsertPt, tmp(1,:)}
						LVepi = roi(ROI).Position{kk,1};
						tmp = bsxfun(@minus,LVepi,ll{:});
						[~,tmp] = sort(sum(tmp.^2,2),'ascend');
						idx_NearInsPt = min(tmp([1 2]))+1;
						idx_NearInterpInsPt = dsearchn(LVepiInterp,LVepi(tmp([1 2]),:));
						
						nPt = length(LVepi)+1;
						tmp = [1:nPt]; tmp(idx_NearInsPt) = [];
						roi(ROI).Position{kk,1} = nan(nPt,2);
						roi(ROI).Position{kk,1}(tmp,:) = LVepi;
						roi(ROI).Position{kk,1}(idx_NearInsPt,:) = LVepiInterp(round(mean(idx_NearInterpInsPt)),:);
					end
					
					% Added point & RVinsertPt could be the different side:
					%{ 
					tmp = bsxfun(@minus,LVepi,InsertPt);
					[~,idx_NearAntInsPt] = min(sum(tmp.^2,2));
					tmp = bsxfun(@minus,LVepiInterp,LVepi(idx_NearAntInsPt,:));
					[~,idx_NearInterpAntInsPt1] = min(sum(tmp.^2,2));
					% tmp(1,:): LVepi(idx_NearAntInsPt+clcDir,:)
					tmp = circshift(LVepi,-(idx_NearAntInsPt-1-1));
					tmp = bsxfun(@minus,LVepiInterp,tmp(1,:));
					[~,idx_NearInterpAntInsPt2] = min(sum(tmp.^2,2));
					
					tmp = circshift(LVepiInterp,-(idx_AntInsPt-1+clcDir*SegPos));
					tmp = bsxfun(@minus,LVepi,tmp(1,:));
					[~,idx_NearPosInsPt] = min(sum(tmp.^2,2));
					tmp = bsxfun(@minus,LVepiInterp,LVepi(idx_NearPosInsPt,:));
					[~,idx_NearInterpPosInsPt1] = min(sum(tmp.^2,2));
					% tmp(1,:): LVepi(idx_NearPosInsPt-clcDir,:)
					tmp = circshift(LVepi,-(idx_NearPosInsPt-1+1));
					tmp = bsxfun(@minus,LVepiInterp,tmp(1,:));
					[~,idx_NearInterpPosInsPt2] = min(sum(tmp.^2,2));
					
					if clcDir == -1
						idx_NearInsPt = [idx_NearAntInsPt+1,idx_NearPosInsPt+1];
					else
						idx_NearInsPt = [idx_NearAntInsPt+2,idx_NearPosInsPt];
					end
					
					nPt = length(LVepi)+length(idx_NearInsPt);
					tmp = [1:nPt];
					tmp(idx_NearInsPt) = [];
					roi(ROI).Position{kk,1} = nan(nPt,2);
					roi(ROI).Position{kk,1}(tmp,:) = LVepi;
					roi(ROI).Position{kk,1}(idx_NearInsPt,:) = [LVepiInterp(round((idx_NearInterpAntInsPt1+idx_NearInterpAntInsPt2)/2),:);LVepiInterp(round((idx_NearInterpPosInsPt1+idx_NearInterpPosInsPt2)/2),:)];
					 %}
					
					%% Generate RVendo:
					% Shift array circularly but fail to reverse the order
					%{ 
					idx3 = circshift(1:60,idx_AntInsPt-clcDir*self.aside.ant);
					LVepiInterp_ReOrd = LVepiInterp(idx3);
					 %}
					idx3 = idx_AntInsPt+clcDir*self.aside.ant;
					if idx3 < 1
						idx3 = 60 + idx3;
					elseif idx3 > 60
						idx3 = idx3 - 60;					
					end
					if clcDir == -1
					% if strcmpi(dir_clc,'CW')
					%% LV points stored CW under rightward X-axis and downward Y-axis: 
						LVepiInterp_ReOrd = LVepiInterp([idx3:-1:1,end:-1:idx3+1],:);
					else
					%% LV points stored CCW under rightward X-axis and upward Y-axis(common coord. sys.):
						LVepiInterp_ReOrd = LVepiInterp([idx3:end,1:idx3-1],:);
					end
					RVendo(1,:) = LVepiInterp_ReOrd(1,:);
					SegPos_new = SegPos - self.aside.ant + self.aside.pos;
					if SegPos_new < 1
						msgbox(['Wrong selections: Posterior-aside-distance ',num2str(self.aside.pos),'is too large in magnitude compared with ','Anterior-aside-distance ',num2str(self.aside.ant)]);
						return				
					end
					RVendo(nRVendoPt/2+1,:) = LVepiInterp_ReOrd(SegPos_new,:);
					interval = (SegPos_new-1)/(nRVendoPt/2);
					% interval = round((SegPos_new-1)/(nRVendoPt/2));
					%{ 
					interval = (SegPos_new-1)/(nRVendoPt/2);
					if mod(interval,1) > 0.6
						interval = ceil(interval);					
					else
						interval = floor(interval);
					end
					 %}
					LVendo = roi(jj).Position{kk, 2};
					center = mean(LVendo);
					% center = mean([mean(LVepi);mean(LVendo)])=AnalysisInfo.PositionA+[dx,dy]?
					% tmp = bsxfun(@minus,LVendo,center); rad = mean(sqrt(sum(tmp.^2,2)));
					%% Project Direction:
                    v1 = (RVendo(nRVendoPt/2+1,:)+RVendo(1,:))/2-center;
					% v1 = RVendo(nRVendoPt/2+1,:)-RVendo(1,:); v1 = [v1(2),-v1(1)];
					v1 = v1/norm(v1);
					v2 = (RVendo(nRVendoPt/2+1,:)-RVendo(1,:));
					for ll = 1:nRVendoPt/2-1
						RVendo(1+ll,:) = LVepiInterp_ReOrd(round(1+interval*ll),:);
						% RVendo(end+1-ll,:) = RVendo(1+ll,:) + v1*rad;
						%% aside-distant of RV freewall Pt = dist of the Pt on the arc to the corresponding chord
						tmp = dot(v2, RVendo(1+ll,:) - RVendo(1,:)) / dot(v2, v2);
						intersectPt = RVendo(1,:) + tmp * v2;
						dist = norm(RVendo(1+ll,:)-intersectPt);
						RVendo(end+1-ll,:) = RVendo(1+ll,:) + v1*dist;
					end
					roi(ROI).Position{kk,3} = RVendo;
					
					% if kk == 3; breakpoint; end;
					
%{ 			
% Validate RVendo generation:
figure
% plot(LVepi(:,1),LVepi(:,2),'LineWidth',1,'color','b','Marker','x','MarkerSize',5);
plot(roi(ROI).Position{kk,1} (:,1),roi(ROI).Position{kk,1} (:,2),'LineWidth',1,'color', 'b','Marker', 'x', 'MarkerSize',5); 
hold on;
plot(LVendo(:,1),LVendo(:,2),'LineWidth',1,'color', 'r','Marker', '+', 'MarkerSize',5); 
title('Validate RVendo generation!');
plot(RVendo(:,1),RVendo(:,2),'LineWidth',1,'color', 'k','Marker', '+', 'MarkerSize',5); 
set(gca,'Ydir','reverse');
axis equal;
plot(LVepiInterp(1,1),LVepiInterp(1,2),'LineWidth',1,'color', 'k','Marker', 'o', 'MarkerSize',10);
plot(LVepiInterp(4,1),LVepiInterp(4,2),'LineWidth',1,'color', 'k','Marker', 's', 'MarkerSize',10);
plot(LVepiInterp(SegPos,1),LVepiInterp(SegPos,2),'LineWidth',1,'color', 'k','Marker', 'd', 'MarkerSize',10);
% plot(InsertPt(1,1),InsertPt(1,2),'LineWidth',1,'color', 'k','Marker', 'o', 'MarkerSize',10);
% plot(RVendo(1,1),RVendo(1,2),'LineWidth',1,'color', 'k','Marker', 'o', 'MarkerSize',5);
% plot(RVendo(nRVendoPt/2+1,1),RVendo(nRVendoPt/2+1,2),'LineWidth',1,'color', 'k','Marker', 'o', 'MarkerSize',5);		
% plot(AnalysisInfo.PositionB(1),AnalysisInfo.PositionB(2),'color', 'k','Marker', 'o', 'MarkerSize',10);
% plot(ROIInfo.Contour{kk, 1}(1,1),ROIInfo.Contour{kk, 1}(1,2),'LineWidth',1,'color', 'k','Marker', 'x', 'MarkerSize',10);
figure
plot(LVepiInterp(:,1),LVepiInterp(:,2),'LineWidth',1,'color', 'k');
axis equal;
 %}
				end
			end
		end
		end
		
		idx = find(strcmpi({roi.Type},'LA'));
		idx1 = strcmpi({roi(idx).Name},'Auto Generated LA'); idx(idx1) = [];
		if ~isempty(idx)
		for jj = idx
			roiName = [dns(MagIndex == roi(jj).SeqIndex(1)).Name,'_',roi(jj).Name];
			idx1 = find(strcmpi({roi.Name},'Long Axis Dual Ventricle'));
			idx2 = cellfun(@(x)isequal(x,roi(jj).SeqIndex), {roi(idx1).SeqIndex});
			switch sum(idx2)
				case 1
					ROI = idx1(idx2);
				case 0
					nROI = nROI + 1;
					ROI = nROI;
				otherwise
					msgbox(['Duplicate converted ROIs detected:',num2str(sum(idx1)),' converted ROIs found for the LA ROI:',roiName,' of the input',files(ii)]);
					continue
			end
			
			roi(ROI).Name = 'Long Axis Dual Ventricle';
			roi(ROI).Type = 'ladual';
			roi(ROI).UID = dicomuid;
			roi(ROI).SeqIndex = roi(jj).SeqIndex;
			nFrames = size(roi(jj).Position,1);
			roi(ROI).IsClosed = repmat({false}, [nFrames,3]);
			roi(jj).IsCurved = cellfun(@logical,roi(jj).IsCurved,'UniformOutput',false);
			% roi(ROI).IsCurved = [roi(jj).IsCurved,roi(jj).IsCurved(:,2)];
			roi(ROI).IsCurved = roi(jj).IsCurved;
			roi(jj).IsCorner = cellfun(@logical,roi(jj).IsCorner,'UniformOutput',false);
			roi(ROI).IsCorner = roi(jj).IsCorner;
			
			roi(ROI).Position = roi(jj).Position;			
			idx1 = cellfun(@isempty, roi(jj).Position);
			if ~isequal(idx1(:,1),idx1(:,2))
				msgbox(['Missing LV-endo or LV-epi countours on the LA ROI:',roiName,' for the input',files(ii)]);			
			end
			idx1 = find(~idx1(:,1))';
			if isempty(idx1)
				if ROI < numel(roi)
					msgbox(['No contours found on the LA ROI:',roiName,' for the input',files(ii),'But it has been converted before!']);
					return
				else
					roi(ROI) = [];
					nROI = nROI - 1;
				end				
				continue;
			else
				for kk = idx1			
					LVepi = roi(jj).Position{kk, 1};
					LVendo = roi(jj).Position{kk, 2};
					EpiApex = size(LVepi,1); EpiApex = (EpiApex+1)/2;
					EndoApex = size(LVendo,1); EndoApex = (EndoApex+1)/2;
					if floor(EpiApex)~=EpiApex || floor(EndoApex)~=EndoApex
					%% NOT FAIL for inf: 
					% if mod(EpiApex,1)~=0 || mod(EndoApex,2)~=0
						msgbox(['CANNOT FOUND the APEX because of the even number of LV points on LA ROI ',roiName,' for the input ',files(ii)],'A update is required for the Codes!');
						continue
					end
					nRVendoPt = 2*(EpiApex-self.aside.LA)-1;%size(LVepi,1)-2*self.aside.LA;
					roi(ROI).IsCurved{kk, 3} = repmat(true, [nRVendoPt,1]);
					roi(ROI).IsCorner{kk, 3} = repmat(false, [nRVendoPt,1]);
					% maximum curvature here: err for max(inf)?
					roi(ROI).IsCorner{kk, 3}(EpiApex-self.aside.LA) = true;
					
					if LVepi(1,1)-LVepi(EpiApex,1) < LVepi(EpiApex,1)-LVepi(EpiApex,1)
						LVepi = LVepi(end:-1:1,:);
						LVendo = LVendo(end:-1:1,:);
					end
					v2 = LVepi(EpiApex,:)-LVendo(EndoApex,:);
					RVendo = [LVepi(1:EpiApex-self.aside.LA,:);zeros(EpiApex-1-self.aside.LA,2)];
					for ll = 1:EpiApex-1-self.aside.LA
						tmp = dot(v2, LVendo(ll,:) - LVendo(EndoApex,:)) / dot(v2, v2);
						intersectPt = LVendo(EndoApex,:) + tmp * v2;
						dist = norm(LVendo(ll,:)-intersectPt);
						v1 = LVepi(ll,:) - LVendo(ll,:); v1 = v1/norm(v1);
						RVendo(end+1-ll,:) = LVepi(ll,:) + v1*dist;
						thk = norm(LVepi(ll,:) - LVendo(ll,:));
						roi(ROI).Position{kk,1}(ll,:) = RVendo(end+1-ll,:) + v1*thk/3;
					end
					ll = EpiApex-self.aside.LA;
					v1 = LVepi(ll,:) - LVendo(ll,:); v1 = v1/norm(v1);
					thk = norm(LVepi(ll,:) - LVendo(ll,:));
					roi(ROI).Position{kk,1}(ll,:) = RVendo(ll,:)+v1*thk/3;
					
					roi(ROI).Position{kk,3} = RVendo;
%{ 			
% Validate RVendo generation:
figure
plot(LVepi(:,1),LVepi(:,2),'LineWidth',1,'color','b','Marker','x','MarkerSize',5);
hold on;
plot(LVendo(:,1),LVendo(:,2),'LineWidth',1,'color', 'r','Marker', '+', 'MarkerSize',5); 
title('Validate RVendo generation!');
plot(RVendo(:,1),RVendo(:,2),'LineWidth',1,'color', 'r','Marker', '+', 'MarkerSize',5); 
set(gca,'Ydir','reverse');
axis equal;

figure
plot(roi(ROI).Position{kk,1}(:,1),roi(ROI).Position{kk,1}(:,2),'LineWidth',1,'color','b','Marker','x','MarkerSize',5);
hold on;
plot(roi(jj).Position{kk, 2}(:,1),roi(jj).Position{kk, 2}(:,2),'LineWidth',1,'color', 'r','Marker', '+', 'MarkerSize',5); 
title('Validate RVendo generation!');
plot(roi(ROI).Position{kk,3}(:,1),roi(ROI).Position{kk,3}(:,2),'LineWidth',1,'color', 'r','Marker', '+', 'MarkerSize',5); 
set(gca,'Ydir','reverse');
axis equal;
 %}
					end
			end
		end
		end
		
		% Make compatible:
		for kk = 1 : nText
			idx = find(strcmpi({roi.Type},old{kk}));
			if ~isempty(idx)
				for jj = idx
					roi(jj).Type = new{kk};%Err if new(kk) will make roi.Type as cell instead of char
				end
			end
		end
		
		%% Save outFiles:
		if strwcmpi(fname{ii},'*_v05')
			outFiles{ii} = inFiles{ii}
		else
			outFiles{ii} = fullfile(pathDNS,[fname{ii},'_v05.dns']);
		end
		
		if nFiles == 1
			while true
				[uifile, uipath] = uiputfile({'*.dns', 'Save a DENSEanalysis(v0.5) workspace (*.dns)'},'Save DNS file', outFiles{ii});
				if isequal(uipath,0) || isequal(uifile,0)
					tmp = questdlg('Do you wanna quit? You might lose data!','Warning');
					if strcmpi(tmp,'Yes')
						return
					else
						continue
					end
				else
					%% Append the correct file extension if needed
					[~,fname] = fileparts(uifile);
					outFiles{ii} = [uipath,fname,'.dns'];
					break
				end
			end
		end
			
		if exist('status','var')
			save(outFiles{ii}, 'seq','img','dns','roi','status');
		else
			save(outFiles{ii}, 'seq','img','dns','roi');	
		end
		
	end
	
end

function Debug
	%% Incorrect Anterior-Epi-RV-insertPT: average the correct ones
	k = 1;
	tmp{k} = AnalysisInfo.PositionB;
	k = k+1;
	AnalysisInfo.PositionB = mean(vertcat(tmp{:}));
end

