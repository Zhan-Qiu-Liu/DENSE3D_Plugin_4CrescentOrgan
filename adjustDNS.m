% adjustDNS  - adjust .dns files of SA and LA
%
%
%	adjust .dns files of LA to SA 
%
%
% USAGE:
%   adjustDNS();
%
%
% See also:
%   
%
% Last Modified: 11:39 PM Monday, November 2, 2015
% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)

% function adjustDNS(varargin)

	%% DIRECTORY LOAD
	% default inputs
	uipath = pwd;
	%{ 
	% locate a SA file to load
	[uifile,uipath] = uigetfile(...
		{'*.dns','Copy DNS from a DENSE workspace(*.dns)'},...
		'Select DNS file',uipath);
	if isequal(uipath,0) || isequal(uifile,0)
		return
	end
	%% test:
	% uifile = '20150823_6SA.dns';
	filename = fullfile(uipath,uifile);
	% load file
	tmp1 = load(filename,'-mat');
	if ~all(isfield(tmp1,{'seq','img','dns','roi'}))
		msgbox('".dns" file does not contain the necessary information');
		return
	end

	% locate a LA file to load
	[uifile,uipath] = uigetfile(...
		{'*.dns','Copy DNS to a DENSE workspace(*.dns)'},...
		'Select DNS file',uipath);
	if isequal(uipath,0) || isequal(uifile,0)
		return
	end
	%% test:
	% uifile = '20150823_LA.dns';
	filename = fullfile(uipath,uifile);
	% load file
	tmp1 = load(filename,'-mat');
	if ~all(isfield(tmp1,{'seq','img','dns','roi'}))
		msgbox('".dns" file does not contain the necessary information');
	end
	 %}
	 
	
for kk = 1:10	
	% tmp1 = load('vol1_scan1.dns','-mat'); tmp1 = load('vol1_scan2.dns','-mat');
	inFile = ['vol',num2str(kk),'_SectionA.dns'];
	outFile = ['vol',num2str(kk),'_SectionA_redo.dns'];
	tmp1 = load(inFile,'-mat');
%{ 
tmp = dir(uipath); tmp = {tmp.name};
temp = strncmpi(tmp, '.dns', -4);
for ii = 1:size(tmp,2)
	if temp2(ii)
		jj = jj + 1;
		files(jj) = fullfile(path,tmp(ii));%load names
	end
end
%}

	%% adjustDNS  - adjust .dns files of LA to SA 
	nseq1 = numel(tmp1.seq);
	ndns1 = numel(tmp1.dns);
	nroi1 = numel(tmp1.roi);
	
	for ii = 1:nroi1
		tmp1.roi(ii).UID = dicomuid;
		
		%% Add noises to contours:
		[nFr, nContours] = size(tmp1.roi(ii).Position);
		for jj = 1:nContours
			for fr = 1:nFr
				pts = tmp1.roi(ii).Position{fr, jj};
				center = mean(pts);
				tmp = bsxfun(@minus,pts,center); rad = mean(sqrt(sum(tmp.^2,2)));
				[Nx, Ny] = size(pts);
%% Find the best Power : mm
%{ 
for mm = 1:10
	figure; hold on; title(num2str(mm));
	for ll = 1:100
		tmp = rad/100*mm*randn(Nx, Ny) + pts; plot(tmp(:,1),tmp(:,2),'r');
	end
	plot(pts(:,1),pts(:,2),'kx');
end
 %}
				mm = 2;
				tmp1.roi(ii).Position{fr, jj} = rad/100*mm*randn(Nx, Ny) + pts;
%% Plot difference:
% figure; plot(pts(:,1),pts(:,2),'k'); hold on; plot(tmp1.roi(ii).Position{fr, jj}(:,1),tmp1.roi(ii).Position{fr, jj}(:,2),'b'); title(num2str(mm));
			end
		end
	end

	%{ 
	for k = 1:nseq2
		tmp1.seq(nseq1+k) = tmp1.seq(k);
		tmp1.img(nseq1+k) = tmp1.img(k);
	end

	for k = 1:ndns2
		tmp1.dns(ndns1+k) = tmp1.dns(k);
		tmp1.dns(ndns1+k).MagIndex = tmp1.dns(k).MagIndex + nseq1;
		tmp1.dns(ndns1+k).PhaIndex = tmp1.dns(k).PhaIndex + [nseq1,nseq1,nseq1];	
		tmp1.dns(ndns1+k).Name = strcat('auto.',sprintf('%d',ndns1+k));	
	end

	% if nroi2 ~= 0
		for k = 1:nroi2
			tmp1.roi(nroi1+k) = tmp1.roi(k);
			tmp1.roi(nroi1+k).SeqIndex = tmp1.roi(k).SeqIndex + [nseq1,nseq1,nseq1,nseq1];
		end
	% end

	uifile = 0;
	while isequal(uifile,0)
		[uifile, uipath] = uiputfile(...
			{'*.dns', 'Save a DENSE workspace (*.dns)'},...
			[], uipath);
		if isequal(uipath,0) || isequal(uifile,0)
			button = questdlg('Do you wanna quit? You might lose data!','Warning','No');
			if strcmpi(button,'Yes')
				return
			else
				uipath = pwd;
			end
		end
	end

	%% Append the correct file extension if needed
	% Method #1: Can replace an existing one
	filename = fullfile(uipath,uifile);
	% Method #2: Never replace an existing one
	[~,~,ext] = fileparts(filename);
	if ~strcmpi(ext, '.dns')
		filename = strcat(filename, '.dns');
	end
	%}
	
	save(outFile, '-struct', 'tmp1');
	% save(filename,'tmp1.seq','tmp1.img','tmp1.dns','tmp1.roi');
end
