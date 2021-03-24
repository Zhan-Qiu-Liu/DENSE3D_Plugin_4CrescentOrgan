% mergeDNS  - merge .dns files of SA and LA
%
%
%	merge .dns files of LA to SA 
%
%
% USAGE:
%   mergeDNS();
%
%
% See also:
%   
%
% Last Modified: 11:39 PM Monday, November 2, 2015
% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)

function mergeDNS(varargin)

	%% DIRECTORY LOAD
	% default inputs
	uipath = pwd;
		
	% locate a SA file to load
	[uifile,uipath] = uigetfile(...
		{'*.dns','Open a DENSE workspace of SA(*.dns)'},...
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
		{'*.dns','Open a DENSE workspace of LA(*.dns)'},...
		'Select DNS file',uipath);
	if isequal(uipath,0) || isequal(uifile,0)
		return
	end
	%% test:
	% uifile = '20150823_LA.dns';
	filename = fullfile(uipath,uifile);
	% load file
	tmp2 = load(filename,'-mat');
	if ~all(isfield(tmp2,{'seq','img','dns','roi'}))
		msgbox('".dns" file does not contain the necessary information');
	end

	%% mergeDNS  - merge .dns files of LA to SA 
	nseq1 = numel(tmp1.seq);
	ndns1 = numel(tmp1.dns);
	nroi1 = numel(tmp1.roi);
	nseq2 = numel(tmp2.seq);
	ndns2 = numel(tmp2.dns);
	nroi2 = numel(tmp2.roi);

	for k = 1:nseq2
		tmp1.seq(nseq1+k) = tmp2.seq(k);
		tmp1.img(nseq1+k) = tmp2.img(k);
	end

	for k = 1:ndns2
		tmp1.dns(ndns1+k) = tmp2.dns(k);
		tmp1.dns(ndns1+k).MagIndex = tmp2.dns(k).MagIndex + nseq1;
		tmp1.dns(ndns1+k).PhaIndex = tmp2.dns(k).PhaIndex + [nseq1,nseq1,nseq1];	
		tmp1.dns(ndns1+k).Name = strcat('auto.',sprintf('%d',ndns1+k));	
	end

	% if nroi2 ~= 0
		for k = 1:nroi2
			tmp1.roi(nroi1+k) = tmp2.roi(k);
			tmp1.roi(nroi1+k).SeqIndex = tmp2.roi(k).SeqIndex + [nseq1,nseq1,nseq1,nseq1];
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
		
	save(filename, '-struct', 'tmp1');
	% save(filename,'tmp1.seq','tmp1.img','tmp1.dns','tmp1.roi');

end