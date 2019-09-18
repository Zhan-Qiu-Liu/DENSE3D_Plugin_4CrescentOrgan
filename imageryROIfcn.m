function imageryROIfcn(self,handles,flag)
% Copyright (c) 2016 DENSEanalysis Contributors
% The codes below are modified from the function "switchstate" in "DENSEanalysis.m".
% Last Modified: 11:10 September 13, 2019
% Modified By: Zhan-Qiu Liu (lafeir.lew@gmail.com)

	import plugins.DENSE3D_Plugin_4CrescentOrgan.*
	obj = handles.hdata;
	if flag; keyboard; end
	
    if isempty(obj.seq) || isempty(obj.dns)
        return
    end
	
	fidx = numel(self.roi)+1;
	%% if self.roi is a struct:
	% if isempty(self.roi); self.roi = struct('fidx', 1); end
	
	self.roi{fidx}.ridx = handles.hdense.ROIIndex;
    if isempty(self.roi{fidx}.ridx); return; end
	
switch handles.LastTab
case 2
	while true
		tmp = arrayfun(@num2str, 1:numel(handles.hdense.dnsbuf), 'UniformOutput', false);
		tmp = strcat(tmp, '.', reshape(handles.hdense.dnsbuf,1,[]));
		tmp = inputdlg(strcat([sprintf('%s\n',tmp{:}),'Input a MATRIX of the numbers of seq images used to auto-generate ROIs[use Space/Colon/Semicolon(;) to seperate]:']));
		if isempty(tmp)
			self.roi{fidx}.didx = handles.hdense.DENSEIndex;
			break
		end
		tmp = str2num(tmp{:}); tmp = round(abs(tmp)); tmp(tmp==0) = [];
		if ~isempty(tmp)
			self.roi{fidx}.didx = reshape(tmp,1,[]);
			break
		end
	end
% case 1
	% self.roi{fidx}.didx = handles.hdicom.SequenceIndex;
otherwise
	return
end


%% Initialize Values
self.uniform.tags = {'Mag','Xpha','Ypha','Zpha'};
self.uniform.Marker = {'x','+','*','o','^','s','v','d','>','p','<','+','*','o','s','d','x'};
self.uniform.Color = {'k','b','r','c','y','g','b','r','k','g','m','c','b','r','k','g','m','c'};
if ~isfield(self.uniform,'CN'); self.uniform.CN=[]; end
if isempty(self.uniform.CN)
	EncFreq = zeros(0,0);
	for ii = 1:numel(obj.seq)
		try
			EncFreq = [EncFreq obj.seq(ii).DENSEdata.EncFreq];
		catch
			continue
		end
	end
	self.uniform.CN{round(max(EncFreq)*100)}=[];
end
name = [obj.seq(1).PatientName.FamilyName,'-DENSE#',num2str(self.roi{fidx}.didx),'-ROI#',handles.hdense.roibuf{self.roi{fidx}.ridx}]

ind=0;
self.uniform.CoN=[]; self.uniform.magMean=[];self.uniform.magStd=[];
hcurve={[],[]};% count=0;
hfig=figure('units','normalized','outerposition',[0 0 1 1],'Name', [name,'-insID#',num2str(fidx)]);
% hplot = axes; hold(hplot,'on');title('Single Subject');
hplot(1)=subplot(1,2,1); hold(hplot(1),'on'); grid(hplot(1),'on');
title(hplot(1), 'Single Subject with View Sharing');xlabel(hplot(1),'Cardiac Phase Number'); ylabel(hplot(1),'CN');
hplot(2)=subplot(1,2,2); hold(hplot(2),'on'); grid(hplot(2),'on');
title(hplot(2), 'Single Subject without View Sharing'); xlabel(hplot(2),'Cardiac Phase Number'); ylabel(hplot(2),'CN');
% axis(hplot,'equal');
% set(hplot,'XMinorTick','off','YMinorTick','off','XMinorGrid','off','YMinorGrid','off');
% 'XLim',[1,roundn(api.Nfr,1)]

self.uniform.N=[];
hCurve={[],[]};
hFig=figure('units','normalized','outerposition',[0 0 1 1],'Name', [name,'-insID#',num2str(fidx)]);
% hPlot = axes; hold(hPlot,'on');title('Single Subject');
hPlot(1)=subplot(1,2,1); hold(hPlot(1),'on'); grid(hPlot(1),'on');
title(hPlot(1), 'Single Subject with View Sharing');xlabel(hPlot(1),'Cardiac Phase Number'); ylabel(hPlot(1),'Mean of Magnitude');
hPlot(2)=subplot(1,2,2); hold(hPlot(2),'on'); grid(hPlot(2),'on');
title(hPlot(2), 'Single Subject without View Sharing'); xlabel(hPlot(2),'Cardiac Phase Number'); ylabel(hPlot(2),'N');
% set(hPlot,'YLim',[0,1]); axis(hPlot,'equal'); 

EncFreq = zeros(0,0); 
for didx = self.roi{fidx}.didx
	% image data
	api = imagedata(obj,didx);
	cndata = contourDataREPL(obj,self.roi{fidx}.ridx);
	fld = fieldnames(cndata);
	for ii = 1: numel(fld)
		api.(fld{ii}) = cndata.(fld{ii});
	end
	api.SeedFrame = handles.hdense.Frame;
	
    % stage 2 check
    xyzvalid = cellfun(@(tag)~isempty(api.(tag)),self.uniform.tags);
    idx = find(xyzvalid,1,'first');
    if ~any(xyzvalid)
        warning('All phase fields are empty!');
		keyboard
    end
	
%% Comment this if images in SOI are not in same size:
% if ~isfield(self.roi{fidx},'Mask') || isempty(self.roi{fidx}.Mask)
	% image size / number of frames
	api.Isz = size(api.(self.uniform.tags{idx})(:,:,1));
	api.Nfr = size(api.(self.uniform.tags{idx}),3);
	% image space
	x = 1:api.Isz(2);
	y = 1:api.Isz(1);
	[X,Y] = meshgrid(x,y);

	%%% input self.roi{fidx}.Mask or Contour/MaskFcn check
	% check contour data
	C = api.Contour;
	if isempty(C) || ~iscell(C) || ~ismatrix(C) || size(C,1) ~= api.Nfr
		warning('Invalid Contour data.');
		keyboard
	end

	tfvalid = cellfun(@(c)(~isempty(c) && isnumeric(c) && ...
		ismatrix(c) && size(c,2)==2), C);
	tfvalid = all(tfvalid,2);

	% MaskFcn
	maskfcn = api.MaskFcn;
	% maskfcn = @maskGeneral;
	if isempty(maskfcn) || ~isa(maskfcn,'function_handle')
		warning('Invalid MaskFcn.');
		keyboard
	end

	try
		self.roi{fidx}.Mask = false([api.Isz,api.Nfr]);
		for fr = 1:api.Nfr
			if tfvalid(fr)
				self.roi{fidx}.Mask(:,:,fr) = maskfcn(X,Y,C(fr,:));
			end
		end
	catch ERR
		ME = MException('The MaskFcn failed.');
		ME.addCause(ERR);
		throw(ME);
	end

	% check Mask
	checkfcn = @(m)islogical(m) && any(ndims(m)==[2 3]) && ...
		all(size(m(:,:,1))==api.Isz) && size(m,3)==api.Nfr;
	if ~checkfcn(self.roi{fidx}.Mask)
		warning('Invalid Mask');
		keyboard
	end

	% all valid frames
	tfvalid = squeeze(any(any(self.roi{fidx}.Mask)));
	frvalid = find(tfvalid);
	frvalid = frvalid([1 end])'
	if ~all(tfvalid(frvalid(1):frvalid(end)))
		warning('The ROI is not defined on a continuous frame range.');
		keyboard
	end
	if abs(frvalid(1)-api.ValidFrames(1))>0.1 || abs(frvalid(2)-api.ValidFrames(2))>0.1
		warning('Mismatch for valid frames? Valid frames DENSEanalysis calculated:');
		api.ValidFrames
		keyboard
	end
% end
	
    %%%% display imagery
	% if isfield(self.roi{fidx},'fig')
	try
		figure(self.roi{fidx}.fig);
	% else	
	catch
		self.roi{fidx}.fig = figure(...
		'Name', ['ROI Imagery | created by Zhan-Qiu:#',[name,'-insID#',num2str(fidx)]], ...
		'NumberTitle',   'off', ...
		'color', self.BackgroundColor, ...%'none'
		'units', 'normalized',...
		'outerposition',[0 0 1 1],...
		'renderer', 'opengl',...
		'CloseRequestFcn',@(s,e)deleteFcn(self,fidx));
		% 'KeyPressFcn', @windowkeypress, ...% CANNOT pass "self" to windowkeypress
		% 'MenuBar',       'none', ...
		% 'Toolbar',       'none');
		self.roi{fidx}.delete = @deleteFcn;

		self.roi{fidx}.hax(1) = subplot(2,2,1);  
		self.roi{fidx}.hax(2) = subplot(2,2,2);
		self.roi{fidx}.hax(3) = subplot(2,2,3);
		self.roi{fidx}.hax(4) = subplot(2,2,4);
		
		%% playbar & listeners
		self.roi{fidx}.hplaybar = playbar(self.roi{fidx}.fig);
		% place the playbar at the bottom/self.roi{fidx}.center of the display panel
		self.roi{fidx}.hplaybar.Units = 'normalized';
		self.roi{fidx}.hplaybar.Position = [0.40 0.01 0.15 0.04];
		self.roi{fidx}.hplaybar.Min = frvalid(1);
		self.roi{fidx}.hplaybar.Max = frvalid(2);
		self.roi{fidx}.hplaybar.Enable = 'on';
		self.roi{fidx}.hplaybar.Value = frvalid(1);
		% Most Anonymous Functions have two input arguments by default:
		self.roi{fidx}.hlisten_playbar = addlistener(self.roi{fidx}.hplaybar,'NewValue',@(s,e,x,y)playbackFcn(self,obj,api,fidx));%@(varargin)playbackFcn(varargin{1}.Value)
		self.roi{fidx}.hlisten_playbar.Enabled = true;
		
		%% display first magnitude image
		if ~isempty(api.Mag)
			self.roi{fidx}.him(1)  = imshow(rand(api.Isz),[0 1],'init','fit','parent',self.roi{fidx}.hax(1));
			
			title(self.roi{fidx}.hax(1), self.uniform.tags{1});
			hold(self.roi{fidx}.hax(1),'on');
		end
		
		%% Outline ROI in phase 
		for k = 2:numel(self.uniform.tags)
			if ~isempty(api.(self.uniform.tags{k}))
				idx = numel(self.roi{fidx}.him)+1;
				% 'CData': a matrix of [Isz(1) X Isz(2)X 3] for TrueColor contains element in range 0~1
				self.roi{fidx}.him(idx)=image('parent',self.roi{fidx}.hax(idx),'cdata',rand([api.Isz,3]),'hittest','off');
				title(self.roi{fidx}.hax(idx), self.uniform.tags{k});
				hold(self.roi{fidx}.hax(idx),'on');		
			end
		end
		%{ 
		if ~isempty(api.Xpha)
			rgb = im2rgb(api.Xpha(:,:,fr),tf);
			self.roi{fidx}.him(2)=image('parent',self.roi{fidx}.hax(2),'cdata',rgb,'hittest','off');
			title(self.roi{fidx}.hax(2), self.uniform.tags{2});
			hold(self.roi{fidx}.hax(2),'on');
		end	
		if ~isempty(api.Ypha)
			rgb = im2rgb(api.Ypha(:,:,fr),tf);
			self.roi{fidx}.him(3)=image('parent',self.roi{fidx}.hax(3),'cdata',rgb,'hittest','off');
			title(self.roi{fidx}.hax(3), self.uniform.tags{3});
			hold(self.roi{fidx}.hax(3),'on');
		end
		if ~isempty(api.Zpha)
			rgb = im2rgb(api.Zpha(:,:,fr),tf);
			self.roi{fidx}.him(4)=image('parent',self.roi{fidx}.hax(4),'cdata',rgb,'hittest','off');
			title(self.roi{fidx}.hax(4), self.uniform.tags{4});
			hold(self.roi{fidx}.hax(4),'on');
		end
		 %}
		 
		if isempty(api.Zpha)
			self.roi{fidx}.flag.mag2 = true;
			self.roi{fidx}.him(4) = image('parent',self.roi{fidx}.hax(4),'cdata',rand([api.Isz,3]),'hittest','off');
			set(self.roi{fidx}.hax(4),'xlim',[0 api.Isz(2)]+0.5,'ylim',[0 api.Isz(1)]+0.5);
		else
			self.roi{fidx}.flag.mag2 = false;
		end		
		axis(self.roi{fidx}.hax,'equal');
		
		%% ROI contours
		self.roi{fidx}.hcline = cline;
		self.roi{fidx}.hcline.UndoEnable = true;	
		self.roi{fidx}.himcardiac = arrayfun(@(hax)imcardiac(self.roi{fidx}.hcline,hax),self.roi{fidx}.hax,'uniformoutput',0);
		self.roi{fidx}.himcardiac = [self.roi{fidx}.himcardiac{:}];
	%   ContextOpenClosed.........allow open/closed contour
	%   ContextSmoothCorner.......allow smooth/corner positions
	%   ContextStraightCurved.....allow straight/curved line segments
	%   ContextAdd................allow point addition
	%   ContextDelete.............allow point deletion
	%   IndependentDrag....enable/disable cline children independent drag
		[self.roi{fidx}.himcardiac.ContextOpenClosed] = deal('on');
		% display points as dots:
		[self.roi{fidx}.himcardiac.Enable] = deal('off');
		[self.roi{fidx}.himcardiac.Enable] = deal('on');
		[self.roi{fidx}.himcardiac.Visible] = deal('on');
		% imcline is general version of imcardiac:
		%{
		self.roi{fidx}.himcline = arrayfun(@(hax)imcline(self.roi{fidx}.hcline,hax),self.roi{fidx}.hax,'uniformoutput',0);
		self.roi{fidx}.himcline   = [self.roi{fidx}.himcline{:}];
		[self.roi{fidx}.himcline.Enable]   = deal('off');
		[self.roi{fidx}.himcline.ContextOpenClosed] = deal('on');
		[self.roi{fidx}.himcline.Visible]   = deal('on');
		%}

%% CANNOT pass the updated "api" to playbackFcn
		self.roi{fidx}.hplaybar.Value = api.SeedFrame;
		% Show or hide Plotted objects
		figure(self.roi{fidx}.fig);
		% plotbrowser(self.roi{fidx}.fig,'off');
		% plotedit(self.roi{fidx}.fig, 'off');
	end
	

	%%% CoV:	
% test
% for fr = 1:api.Nfr
	% Mag = obj.img{api.MagIndex(1)}(:,:,fr);
	% mx(fr) = max(max(Mag));
	% mi(fr) = min(min(Mag));
% end
	for fr = 1:api.Nfr
		tf = self.roi{fidx}.Mask(:,:,fr);
		Mag = double(obj.img{api.MagIndex(1)}(:,:,fr));
%% For DISPLAYING: im2double rescales the output from integer data types to the range [0, 1]
% Mag = im2double(Mag);
		% No rescales

		% normalized magnitude values w.r.t. the max value at each frame
		%{ 
		Mag = api.Mag(:,:,fr);
		if sum(min(Mag)<0) || sum(max(Mag)>1)
			warning('Mag contains element out of range 0~1');
			keyboard
		Mag = api.Mag(:,:,fr);
		end
		 %}
% set(self.roi{fidx}.him(1),'cdata',tf.*Mag/double(mx(fr)));
		Mag = Mag(tf);
		% Mag(~tf)=deal(NaN); idx=find(tf); Mag=Mag(idx) % same as above
		CN(fr) = std(Mag)/mean(Mag)*100;
		magMean(fr)=mean(Mag); magStd(fr)=std(Mag);
		CN(fr) = std(Mag)/mean(Mag)*100;
		mx=max(Mag);mi=min(Mag);
		N(fr)=(mx-mi)/(mx+mi)*100;
	end
	self.uniform.CoN=[self.uniform.CoN;CN];
	self.uniform.magMean=[self.uniform.magMean;magMean];self.uniform.magStd=[self.uniform.magStd;magStd];
	ind=ind+1;
	hcurve{1}(end+1) = plot(hplot(1),CN,'LineWidth',1,'color', self.uniform.Color{ind},'Marker', self.uniform.Marker{ind},'MarkerSize',5);
	hcurve{2}(end+1) = plot(hplot(2),[1:2:api.Nfr],CN([1:2:end]),'LineWidth',1,'color', self.uniform.Color{ind},'Marker', self.uniform.Marker{ind},'MarkerSize',5);
	
	self.uniform.N=[self.uniform.N;N];
	hCurve{1}(end+1) = plot(hPlot(1),magMean,'LineWidth',1,'color', self.uniform.Color{ind},'Marker', self.uniform.Marker{ind},'MarkerSize',5);
	hCurve{2}(end+1) = plot(hPlot(2),[1:2:api.Nfr],magMean([1:2:end]),'LineWidth',1,'color', self.uniform.Color{ind},'Marker', self.uniform.Marker{ind},'MarkerSize',5);
	
	idxSeq = api.PhaIndex(1);
	encFreq = obj.seq(idxSeq).DENSEdata.EncFreq;
	EncFreq = [EncFreq encFreq];
	ke = round(encFreq*100);
	self.uniform.CN{ke}(end+1,:) = CN;
	% self.uniform.CN{ke} = [self.uniform.CN{ke};CN];
	
	idx = numel(self.uniform.seq)+1;
	self.uniform.seq(idx).CN = CN;
	self.uniform.seq(idx).EncFreq = encFreq;
	self.uniform.seq(idx).SliceThickness=obj.seq(idxSeq).SliceThickness;
	self.uniform.seq(idx).PixelSpacing=obj.seq(idxSeq).PixelSpacing;
	self.uniform.seq(idx).EchoTime=obj.seq(idxSeq).EchoTime;
	self.uniform.seq(idx).CardiacNumberOfImages=obj.seq(idxSeq).CardiacNumberOfImages;
	
	
	% h = msgbox({'Compare a series of fittings of generated RV endo meshes to the actual RV endo contours','Write down the best fitting parameters','Click any button when you are ready!'},'Check!');
	% waitfor(h);
	legend(hcurve{1}, arrayfun(@num2str,EncFreq,'UniformOutput', false));
	legend(hcurve{2}, arrayfun(@num2str,EncFreq,'UniformOutput', false));
	
	% legend(hCurve{1}, arrayfun(@num2str,EncFreq,'UniformOutput', false));
	legend(hCurve{2}, arrayfun(@num2str,EncFreq,'UniformOutput', false));
end
	set(hplot,'XTick',[1:api.Nfr]);
	set(hPlot,'XTick',[1:api.Nfr]);
	
	self.uniform.MagMean=mean(self.uniform.magMean); self.uniform.MagStd=mean(self.uniform.magStd);
	
	drawnow
	keyboard
	
%% Save:
	filename=fullfile(char(java.lang.System.getProperty('user.home')),'Dropbox','Analysis_Human',name);
	hgexport(hfig,[filename,'.jpg'],hgexport('factorystyle'), 'Format', 'jpeg');
	% hgsave(hfig,tmp);
	hgexport(hFig,[filename,'-mean.jpg'],hgexport('factorystyle'), 'Format', 'jpeg');
	
	out=self.uniform;
	save([filename,'.mat'],'out');
	return
	
	self.uniform.CN=[];
	CN = vertcat(self.uniform.seq.CN);
	
			ave.(fields{jj}){slice,layer}(region,:) = tmp2;
			stderr.(fields{jj}){slice,layer}(region,:) = tmp3;
			hb = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{region+3},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
			% Adjust error bar width: only work in 2013b
			errorbar_tick(hb, 150);%140-(-1)^layer*20
			hold on;	
	
	%% Interactively draw Contour
	% gather data, based on object type
	switch lower(obj.Type)
		case 'sa'
			[epi,endo] = getSA(hax);
			pos   = {epi,endo};
			iscls = {true};
			iscrv = {true};
			iscrn = {false};
		case 'safull'
			[epi,rendo,lendo] = getSAFull(hax);
			pos   = {epi,rendo,lendo};
			iscls = {true};
			iscrv = {true};
			iscrn = {false};
		case 'lafull'
			[epi,rendo,lendo] = getLAFull(hax);
			pos   = {epi,rendo,lendo};
			iscls = {false};
			iscrv = {true};
			iscrn = {false};
		case 'la'
			[epi,endo] = getLA(hax);
			pos   = {epi,endo};
			iscls = {false};
			iscrv = {true};
			iscrn = {false};
		otherwise
			iscrv = ~strcmpi(obj.Type,'line');
			h = getcline(hax,'IsCurved',iscrv,...
				'color','b','markersize',10,'linewidth',2);
			pos   = h.Position;
			iscls = h.IsClosed;
			iscrv = h.IsCurved;
			iscrn = h.IsCorner;
			delete(h);
	end
		
    self.roi{fidx} = struct(...
        'XYZValid',         xyzvalid,...
        'mmperpixel',       pxsz,...
        'Xwrap',            Xwrap,...
        'Ywrap',            Ywrap,...
        'Zwrap',            Zwrap,...
        'Xunwrap',          Xunwrap,...
        'Yunwrap',          Yunwrap,...
        'Zunwrap',          Zunwrap,...
        'Multipliers',      [xfac,yfac,zfac],...
        'spldx',            spldx,...
        'spldy',            spldy,...
        'spldz',            spldz,...
        'irng',             yrng,...
        'jrng',             xrng,...
        'frrng',            frrng,...
        'ResampleMethod',   method,...
        'ResampleDistance', dresamp,...
        'SpatialSmoothing', pspace,...
        'TemporalOrder',    Npoly,...
        'Xseed',            xseed,...
        'Yseed',            yseed,...
        'Zseed',            zseed,...
        'Contour',          {C},...
        'RestingContour',   {C0},...
        'MaskFcn',          maskfcn);
end

function RGB = im2rgb(im,tf)
%% MASKED RGB IMAGE
    % check for empty
    if isempty(im), RGB = []; return; end

    % "black" color where (tf==0)
    blk = [0.25 0.25 0.50];

    % "white" color where (tf==0)
    wht = [0.50 0.50 1.00];

    % normalize image
    mn = -pi;%min(im(:));
    mx = pi;%max(im(:));
    im = (im - mn)/(mx-mn);

    % generate RGB image
    R = tf.*im + ~tf.*((im*(wht(1)-blk(1))) + blk(1));
    G = tf.*im + ~tf.*((im*(wht(2)-blk(2))) + blk(2));
    B = tf.*im + ~tf.*((im*(wht(3)-blk(3))) + blk(3));
    RGB = cat(3,R,G,B);
end

%% MASK FUNCTIONS
%   In order to recover a black and white region-of-interest self.roi{fidx}.Mask from the
%   contour definition, the user must enter a MASKFCN in the api input.
%   This function must accept a matrix of X/Y pixel values, defined as:
%       api.Isz = size(api.Mag(:,:,1);
%       [X,Y] = meshgrid(1:api.Isz(2),1:api.Isz(1));
%   as well as a contour input "C" (cell array of contours). The function
%   must output a black-and-white logical self.roi{fidx}.Mask defining the region of
%   interest. Thus, the function for the frame "fr" is of the form:
%       function BW = maskfcn(X,Y,C(fr,:))
%
%   For example, consider the case when the contour C is a single cell
%   containing an [Nptx2] contour.  A valid MASKFCN could be:
%       maskfcn = @(X,Y,C)inpolygon(X,Y,C{1}(:,1),C{1}(:,2));
function tf = maskSA(X,Y,C)
    [inep,onep] = inpolygon(X,Y,C{1}(:,1),C{1}(:,2));
    [inen,onen] = inpolygon(X,Y,C{2}(:,1),C{2}(:,2));
    tf = (inep & ~inen) | onep | onen;
end
function tf = maskLA(X,Y,C)
    C = cat(1,C{:});
    tf = inpolygon(X,Y,C(:,1),C(:,2));
end
function tf = maskSAFull(X,Y,C)
    [inep,onep] = inpolygon(X,Y,C{1}(:,1),C{1}(:,2));
    [inRen,onRen] = inpolygon(X,Y,C{2}(:,1),C{2}(:,2));
    [inLen,onLen] = inpolygon(X,Y,C{3}(:,1),C{3}(:,2));
    tf = (inep & ~inRen & ~inLen) | onep | onRen | onLen;
end
function tf = maskLAFull(X,Y,C)
    C = cat(1,C{1},C{3},C{2});
    tf = inpolygon(X,Y,C(:,1),C(:,2));
end
function tf = maskGeneral(X,Y,C)
    tf = false(size(X));
    for n = 1:numel(C)
        tf = tf | inpolygon(X,Y,C{n}(:,1),C{n}(:,2));
    end
end


function deleteFcn(self,fidx)
	try delete(self.roi{fidx}.fig(ishandle(self.roi{fidx}.fig))); end
	self.roi{fidx} = [];
	
	% try rmfield(self.viewerObj.cache,'regionalStrain'); end
	
	try delete(findall(0,'type','figure','-and','tag','WaitbarTimer')); end
end

function playbackFcn(self,obj,api,fidx)
%% Copyright (c) of the following section of the codes: Zhan-Qiu Liu (lafeir.lew@gmail.com)
% Modified By: Zhan-Qiu Liu (lafeir.lew@gmail.com)
% Last Modified: 19:32 June 27, 2018
% Debug Meshes in the current frame:
	%% Image Frames:
	self.roi{fidx}.frame = self.roi{fidx}.hplaybar.Value;
	%% Referential self.roi{fidx}.Frame + Image Frames:
	% self.roi{fidx}.Frame = self.roi{fidx}.frame + 1;
    fr = self.roi{fidx}.frame;
	
	% update cLine to new positions
    self.roi{fidx}.hcline.reset(...
        'position',     obj.roi(self.roi{fidx}.ridx).Position(fr,:),...
        'isclosed',     obj.roi(self.roi{fidx}.ridx).IsClosed(fr,:),...
        'iscurved',     obj.roi(self.roi{fidx}.ridx).IsCurved(fr,:),...
        'iscorner',     obj.roi(self.roi{fidx}.ridx).IsCorner(fr,:));

	%% Show ROI only in magnitude
    tf = self.roi{fidx}.Mask(:,:,fr);
	set(self.roi{fidx}.him(1),'cdata',tf.*api.Mag(:,:,fr));
	% set(self.roi{fidx}.him(1),'Visible','on','HitTest','off');
	
	if self.roi{fidx}.flag.mag2
		set(self.roi{fidx}.him(4),'cdata',api.Mag(:,:,[fr fr fr]));
	end
	% test:
	% hl(1) = line('Parent',self.roi{fidx}.hax(1),'XData',api.Contour{1,fr}(:,1),'YData',api.Contour{1,fr}(:,2),'color','w');
	
    %% Outline ROI in phase 
	for k = 2:numel(self.roi{fidx}.him)
        rgb = im2rgb(api.(self.uniform.tags{k})(:,:,fr),tf);
		set(self.roi{fidx}.him(k),'cdata',rgb)
	end
	
	drawnow
end

%% CANNOT pass "self" to windowkeypress
function windowkeypress(~,evnt,self)
% Last Modified: 8:37 PM Friday, November 20, 2015
% Modified By: Zhan-Qiu Liu (lafeir.lew@gmail.com)
	if isempty(evnt.Modifier)
		key = evnt.Key;
	else
		key = [evnt.Modifier{1},evnt.Key];
		% modifiers = sort(evnt.Modifier);
		% key = strcat(sprintf('%s-', modifiers{:}), evnt.Key);
	end
	
	% avoid hplaybar being clear:
	playbar = self.roi{fidx}.hplaybar;
	
	switch key
	% toolbar
		case 'equal'
			zoom(self.roi{fidx}.hax, 2);
		case 'hyphen'
			zoom(self.roi{fidx}.hax, 0.5);
		case 'altp'
		% Toggle for panning images
			% pan(self.roi{fidx}.fig);
			% msgbox ok;
			pan(self.roi{fidx}.fig);
			hmanager = uigetmodemanager(self.roi{fidx}.fig);
			set(hmanager.WindowListenerHandles,'enable','off');
			set(self.roi{fidx}.fig,'WindowKeyReleaseFcn',@(s,e)pan(self.roi{fidx}.fig));
		case 'altd'
		% Toggle for data cursor
			datacursormode(self.roi{fidx}.fig);
			hmanager = uigetmodemanager(self.roi{fidx}.fig);
			set(hmanager.WindowListenerHandles,'enable','off');
			set(self.roi{fidx}.fig,'WindowKeyReleaseFcn',@(s,e)pan(self.roi{fidx}.fig));
		case 'altr'
		% Toggle for 3D rotation tool
			rotate3d(self.roi{fidx}.fig);
			hmanager = uigetmodemanager(self.roi{fidx}.fig);
			set(hmanager.WindowListenerHandles,'enable','off');
			set(self.roi{fidx}.fig,'WindowKeyReleaseFcn',@(s,e)rotate3d(self.roi{fidx}.fig));
		case 'altb'
		% Toggle for brush tool
			brush(self.roi{fidx}.fig);
			hmanager = uigetmodemanager(self.roi{fidx}.fig);
			set(hmanager.WindowListenerHandles,'enable','off');
			set(self.roi{fidx}.fig,'WindowKeyReleaseFcn',@(s,e)brush(self.roi{fidx}.fig));
		% 'command' for Macintosh computers
		case {'controls','command-s'}
			filemenufcn(self.roi{fidx}.fig,'FileSave');
		case {'controlp','command-p'}
			printdlg(self.roi{fidx}.fig);
	%Viewer
		case {'d', 'rightarrow'}%'n', 
			if ~isempty(playbar.Value) && ~playbar.IsPlaying
				playbar.Value = mod(playbar.Value, playbar.Max) + 1;
			end
		case {'a', 'leftarrow'}%'b', 
			if ~isempty(playbar.Value) && ~playbar.IsPlaying
				playbar.Value = mod((playbar.Value - 2), playbar.Max) + 1;
			end
		% case {'w', 'uparrow'}
		% case {'s', 'downarrow'}
		case {'1','2','3','4','5','6','7','8','9','0'}%'numpad1','numpad2','numpad3','numpad4','numpad5','numpad6','numpad7','numpad8','numpad9','numpad0',
		% Change the active frame
			if ~isempty(playbar.Value) && ~playbar.IsPlaying
				sel = str2double(key(end));
				if isequal(sel,0),sel = double(10);end
				if sel <= playbar.Max
					playbar.Value = sel;
				end
			end
		case 'space'
		% Toggle movie playing
			if isempty(playbar.Value)
				return
			end
			if playbar.IsPlaying
				playbar.stop()
			else
				playbar.play()
			end
	end
end			
