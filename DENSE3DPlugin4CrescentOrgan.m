classdef DENSE3DPlugin4CrescentOrgan < plugins.DENSEanalysisPlugin
    % DENSE3DPlugin4CrescentOrgan - A DENSEanalysis plugin
    %
    %   A plugin for analyzing 3D DENSE data in RV
    %
	% 
	% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)
	% Last Modified: 11:54 July 19, 2017
    
	properties (Constant)
        BackgroundColor = [.78 .93 .8];%eye-protected
		% FOI = 'SeedPoints.m';%  {'surfacemesh.m','Generate_Mesh.m','Fit_Surface.m','Psi.m','Adjust_Mesh_Mex.*','*.mat'};
    end
	
    properties
        handles
		Handles
		
		%% 2D DENSE:
		aside = [];
		dns = struct([]);
		% dns = struct('SegPos_idx',[],'SegPos',[]);
        status = struct('SOI',[],'nSA',[],'nLA',[]);
        % SPLINE information
        spl = repmat(struct,[0 1]);
		straindata = [];

		%% 3D DENSE:
		hShowMesh = [];
		hPickSlice = [];
		viewerObj
		dataObj
		configObj
	end
	
    methods

        function self = DENSE3DPlugin4CrescentOrgan(varargin)

			import plugins.DENSE3D_Plugin_4CrescentOrgan.*
			% load plugin config: plugin.json file
            self@plugins.DENSEanalysisPlugin(varargin{:});
						
            handles = guidata(self.hfig(1));
			% viewer = get(handles.hsidebar.CurrentPanel, 'UserData');

			%%VIEWER---------------------------------------------
			%% Add menu items
			parent = findall(handles.hfig, 'tag', 'menu_file');
			parent = get(parent, 'Parent');
			self.handles.menu_append = uimenu('Parent', parent, 'Label', 'Plugin_DENSE3D4CrescentOrgan');
			uimenu('Parent', self.handles.menu_append, 'Label', 'Contour Interpolation', 'Callback', @(s,e)linearInterp(handles.hdata,handles.hdense.hroi), 'Accelerator', 'D');
			uimenu('Parent', self.handles.menu_append, 'Label', 'Auto-build:SA RVendo(LVendo+epi required)', 'Callback', @(s,e)XformDNS_LV2BV(true,fullfile(get(handles.config,'locations.matpath',userdir()),get(handles.config, 'locations.matfile',userdir())),self),'Accelerator','B');%'DENSEanalysis workspace v0.4->v0.5'
			uimenu('Parent', self.handles.menu_append, 'Label', 'Run Analysis & Export DENSE3D Inputs', 'Callback', @(s,e)menu_runanalysis_REPL(self), 'Accelerator', 'A');%'Tag','menu_runanalysis',
			uimenu('Parent', self.handles.menu_append, 'Label', 'Check 2D computed Contours', 'Callback', @(s,e)chk2DcomputedContours(self));%'Tag','menu_chk2Dcontours',
			uimenu('Parent', self.handles.menu_append, 'Label', 'Initialize DENSE3DPlugin4CrescentOrgan', 'Callback', @(s,e)self.initGUI());
			uimenu('Parent', self.handles.menu_append, 'Label', 'Delete the Pre-defined Slice of Interest', 'Callback', @(s,e)self.deleteSOI());
			uimenu('Parent', self.handles.menu_append, 'Label', 'Export GIF Movie', 'Callback', @(s,e)self.saveGIF());

			%% Remap for all click events
            % set(findobj(handles.hfig, 'tag', 'menu_runanalysis'), 'Callback', @(s,e)menu_runanalysis_REPL(self));
            set(findobj(handles.hfig, 'tag', 'menu_new'), 'Callback', @(s,e)loadFcnREPL(self,true));
            set(findobj(handles.hfig, 'tag', 'tool_new'), 'ClickedCallback', @(s,e)loadFcnREPL(self,true));
            set(findobj(handles.hfig, 'tag', 'menu_open'), 'Callback', @(s,e)loadFcnREPL(self,false));
            set(findobj(handles.hfig, 'tag', 'tool_open'), 'ClickedCallback', @(s,e)loadFcnREPL(self,false));
            set(findobj(handles.hfig, 'tag', 'menu_save'), 'Callback', @(s,e,x)saveFcnREPL(self,false));
            set(findobj(handles.hfig, 'tag', 'tool_save'), 'ClickedCallback', @(s,e,x)saveFcnREPL(self,false));
			set(findobj(handles.hfig, 'tag', 'menu_saveas'), 'Callback', @(s,e,x)saveFcnREPL(self,true));
			
			%% Add menu items
			%{ 
			% load without selection GUI
			parent = findall(handles.hfig, 'tag', 'menu_analysis');
			self.handles.menu_chk2Dcontours = uimenu('Parent', parent, 'Label', 'Check 2D computed Contours', 'Tag','menu_chk2Dcontours','Callback', @(s,e)chk2DcomputedContours(self));
			self.handles.menu_setSOI = uimenu('Parent', parent, 'Label', 'Set SOI', 'Tag','menu_setSOI','Callback', @(s,e)setSOI(self,handles));
			
			parent = get(findall(handles.hfig, 'tag', 'menu_file'), 'Parent');
			self.handles.menu_append = uimenu('Parent', parent, 'Label', 'Plugin_DENSE3D_LV');
			% @(s)XformDNS_LV2BV(true): Too many input arguments.
			uimenu('Parent', self.handles.menu_append, 'Label', 'DENSEanalysis workspace v0.4->v0.5', 'Callback', @(s,e)XformDNS_LV2BV(true));
			% uimenu('Parent', self.handles.menu_append, 'Label', 'Summary of Hotkeys', 'Callback', @HotkeysSummary);
			% uimenu('Parent', self.handles.menu_append, 'Label', 'Known Bugs and Issues', 'Callback', @ReadMe);
			 %}
			
			%% Add toolbar items
			%{ 
			cdata = load(fullfile(self.InstallDir,'icons.mat'),'-mat');
			htoolbar = get(findall(handles.hfig, 'tag', 'tool_open'));
			self.handles.tool_saveas = uipushtool(...
			'Parent',htoolbar,...
			'Separator',        'on', ...
			'ClickedCallback',@(s,e,x)saveFcnREPL(self,handles,true),...
			'CData',cdata.tool_saveas,...
			'TooltipString','Save Workspace As(Ctrl+Shift+S)',...
			'BusyAction','cancel',...
			'Interruptible','off',...
			'Tag','tool_saveas');
			
			% link some menu enable to tool enable
			setappdata(self.handles.tool_saveas,'linkMenuToolEnable',linkprop([self.handles.tool_saveas,self.handles.menu_saveas],'Enable'))
			 %}

			return
			
			%{ 
			% Need to update? Check version!
			tmp = dir(self.InstallDir); files = {tmp.name}; files = files([tmp.isdir] == 0);  idx = strwcmpi(files,'*.version'); [~,files] = fileparts(files{idx});
            Configuration(fullfile(self.InstallDir,[self.Version,'.version']));
			if strwcmpi(files,self.Version)
			 %}
			
			self.handles.tool_DataCursor = uitoolfactory(htoolbar,'Exploration.DataCursor');
			set(self.handles.tool_DataCursor,'Separator','on','TooltipString','Data Cursor(Alt+D)');
			self.handles.tool_reload = uipushtool(...
			'Parent',htoolbar,...
			'Separator',        'on', ...
			'ClickedCallback',@(s,e)reloadFcn(self),...
			'CData',cdata.tool_reload,...
			'TooltipString','Reload Workspace(Ctrl+R)',...
			'Tag','tool_reload');
			%{ 
			% plugins.DENSEanalysisPlugin is a valid base class, not DENSEanalysis:
			% NO update passed to roitool: self.roiidx is empty forever
			uipushtool(...
			'Parent',htoolbar,...
			'Separator',        'on', ...
			'ClickedCallback',@(s,e)viewer.hroi.cut(),...
			'CData',cdata.tool_roi_cut,...
			'TooltipString','Cut ROI(Ctrl+X)',...
			'Tag','tool_roi_cut');
			uipushtool(...
			'Parent',htoolbar,...
			'ClickedCallback',@(s,e)copy(viewer),...
			'CData',cdata.tool_roi_copy,...
			'TooltipString','Copy ROI(Ctrl+C)',...
			'Tag','tool_roi_copy');
			uipushtool(...
			'Parent',htoolbar,...
			'ClickedCallback',@(s,e)viewer.hroi.paste(),...
			'CData',cdata.tool_roi_paste,...
			'TooltipString','Paste ROI(Ctrl+V)',...
			'Tag','tool_roi_paste');
			 %}


			%% Bind more keyboard events
			accelerators = struct(...
				'menu_exportmat',        'M', ...
				'menu_exportexcel',       'X', ...
				'menu_exportroi',        'R', ...
				'menu_test', 'T');
			func = @(x,y)set(findobj(handles.hfig, 'tag', x), 'Accelerator', y);
			cellfun(func, fieldnames(accelerators), struct2cell(accelerators));
			
        end

        function h = uimenu(varargin)
		%% Make sure PluginMenu.m DONNOT create a menu entry for this plugin:
            % Enable checkAvailability:
			% h = gobjects(1,1);
            % Disable checkAvailability:
            h = [];		
        end

        function chk2DcomputedContours(self)
			if isempty(self.straindata)
				self.menu_runanalysis_REPL;
			else
				import plugins.DENSE3D_Plugin_4CrescentOrgan.*
				opts = MeshControl(self.straindata);
				drawnow
				if isempty(opts); return; end
				tags = fieldnames(opts);
				for ti = 1:numel(tags)
					self.straindata.(tags{ti}) = opts.(tags{ti});
				end				
			end
		end

		%{ 
        function setSOI(self)
					
			% tmp = [];
			% while isempty(tmp)
				% tmp = inputdlg('Slices of Interest(enter Space/Comma(,)/Semicolon(;)-separated numbers):');
				% tmp = str2num(tmp{:});
			% end
			tmp = {NaN};		
			while any(isnan(tmp{:}))
				tmp = inputdlg(strcat('Current Slices of Interest(SOI): [',num2str(self.status.SOI),']. [enter Comma(,)/Semicolon(;)-separated slice numbers to reset] (EG: if slicce name is "xyz: auto.1", then slice number is "1")'),'Set SOI');
				if isempty(tmp)
					tmp = {self.status.SOI};
				else
					tmp = regexpi(tmp,',|;','split'); tmp = cellfun(@str2double,tmp,'UniformOutput',false);
				end
			end
			% if isempty(tmp) || any(isnan(tmp{:}))
				% h = errordlg(errstr,'','modal');
				% waitfor(h);
			% else
				self.status.SOI = tmp{1};
			% end
		end
		 %}

        function self = initGUI(self)
		
			if isfield(self.Handles,'hShowMesh')
				self.delete(true);
			end
			
			handles = guidata(self.hfig(1));
			% handles = guidata(self.hShowMesh.fig);
						
			str = 'hmanager';
			try
				POI = cellfun(@(x)strcmpi(x,'3D DENSE Analysis'),{handles.(str).Plugins.Name});
			catch
				str = 'hpluginmenu';
				POI = cellfun(@(x)strcmpi(x,'3D DENSE Analysis'),{handles.(str).Plugins.Name});
			end
			tmp = sum(POI);
			if  tmp ~= 1
				errordlg({'The required Plugin ''3D DENSE Analysis'' CANNOT be solely found:',[num2str(tmp),' Plugins with the same name!'],'The latest compatible version of the Plugin ''3D DENSE Analysis'' is 0.1.0','Try to replace "handles.(str).Plugins" in codes'},'Debug Mode','modal');
				keyboard;
			end
			% Create a new DENSE3DPlugin:
			%{ 
            classes = handles.(str).findAllPlugins('plugins.DENSEanalysisPlugin');
			% pluginType = 'plugins.DENSEanalysisPlugin';
			newplugins = feval(classes(1, 1).Name);
			 %}
			try
				self.viewerObj = handles.(str).Plugins(POI).Handles.hviewer;
				if ~isa(self.viewerObj, 'DataViewer');err;end
			catch
				errordlg({'The Object ''self@plugins.dense3D_plugin.DENSE3Dviewer'' CANNOT be found:','Try to replace "handles.(str).Plugins(POI).hviewer" in codes'},'Debug Mode','modal');
				keyboard;
			end
			
			try
				self.dataObj = self.viewerObj.Data;
			catch
				self.dataObj = handles.(str).Plugins(POI).hdense;
			end
			if ~isa(self.dataObj, 'hgsetget')
				errordlg({'The Object ''self@plugins.dense3D_plugin.DENSE3D'' CANNOT be found:','Try to replace "handles.(str).Plugins(POI).hdense" in codes'},'Debug Mode','modal');
				keyboard;
			end
			
			self.configObj = handles.(str).Plugins(POI).Config;
			
            if isempty(self.dataObj.Data)
				buttonEnable = 'off';
				% errordlg('Please add Files into the 3D Workspaces!','','modal');
				% return
			else
				buttonEnable = 'on';
			end
			

			%% CANNOT add "private" folders.
			% addpath(fullfile(handles.(str).Plugins(POI).InstallDir,'private'));
			
			%% Copy folder:
			if ~exist(fullfile(self.InstallDir,'private'),'dir')
				clear mex;
				try
					copyfile(fullfile(handles.(str).Plugins(POI).InstallDir,'private'),fullfile(self.InstallDir,'private'),'f');
					% copyfile(fullfile(self.InstallDir,self.FOI),fullfile(self.InstallDir,'private'));
				catch ERR
					errordlg(ERR.message,'','modal');
					keyboard;
				end
			end
			%% Copy files:
			% All files:
			%{ 
			tmp = dir(fullfile(handles.(str).Plugins(POI).InstallDir,'private'));
			tmp = {tmp(~[tmp.isdir]).name};
			path_from = cellfun(@(x)fullfile(handles.(str).Plugins(POI).InstallDir,'private',x),self.FOI,'UniformOutput',false);%tmp
			cellfun(@(x)copyfile(x,fullfile(self.InstallDir,'private')),path_from);
			 %}
			
			%% It will replace the original DataViewer: DENSE3Dviewer
			%{ 
		function self = DENSE3Dviewer(data, varargin)
            fakedata = DENSEdata;
			% get protected properties: hparent_display, hdisplay, exportaxes, etc
            self.hviewer = self@handles.(str).Plugins(1).Handles.hviewer;
            self.hviewer = self@DataViewer(struct([]), fakedata, handles.(str).Plugins(POI).Handles.hpanel,handles.(str).Plugins(POI).Handles.hresults);
			self.Data = handles.(str).Plugins(POI).hdense;
            self.hlistener = addlistener(handles.(str).Plugins(POI).hdense, 'NewData', @(s,e)self.newDataCallback());
			 %}
			
			%% Initialize configuration of SeedPoints
			% only numbers as variable name: NOT SUPPORT!
			%{ 
			self.handles.seedPt = Configuration(fullfile(self.InstallDir,'SeedPoints.json'));
			% self.handles.seedPt = Configuration(fullfile(handles.(str).Plugins(3, 1).InstallDir,'SeedPoints.json'));
			set(self.handles.seedPt, '20150628.center', [2.51346206948696	7.08576724017031	-3.12295877567370])
			 %}
	
			%% Add toolbar items
			hbuttongroup = self.viewerObj.Handles.hbuttongroup;
			options = get(hbuttongroup,'Children');
						
			% bottoms = linspace(0, 1, numel(options) + 1);
			% pos = get(options(end),'Position');
			height = 1/(numel(options) + 1);
			self.Handles.hShowMesh = uicontrol( ...
						'Parent', hbuttongroup, ...
						'Style', 'radiobutton', ...
						'Units', 'normalized', ...
						'Enable', buttonEnable, ...
						'Value', 0, ...
						'String', 'Refine 3D Mesh', ...
						'Position', [0.05 -height 0.9 height]);
			set(self.Handles.hShowMesh, 'Callback', @(s,e)showMesh(self), 'Value', 0);
			self.Handles.hPickSlice = uicontrol( ...
						'Parent', hbuttongroup, ...
						'Style', 'radiobutton', ...
						'Units', 'normalized', ...
						'Enable', buttonEnable, ...
						'Value', 0, ...
						'String', 'Slice of Interest for Strain Analysis', ...
						'Position', [0.05 -2*height 0.9 height]);
						% 'String', 'Mid-ventricular Polar Strains', ...
			set(self.Handles.hPickSlice, 'Callback', @(s,e)pickSlice(self), 'Value', 0);
			
			% newDataCallback: enable the buttons
			% buttons = self.viewerObj.options;
			% buttons(end+1).String = 'Slice of Interest for Strain Analysis';
			% buttons(end).Fcn = @(s,e)pickSlice(self);
			self.viewerObj.options(end+1).Handle = self.Handles.hShowMesh;
			self.viewerObj.options(end+1).Handle = self.Handles.hPickSlice;
			% self.viewerObj.options = buttons;
			
			%% Remap for all click events
			%{ 
			buttonsNm = arrayfun(@(x)get(x,'String'),options,'UniformOutput',false);
			BOI = cellfun(@(x)strcmpi(x,'3D Mesh'),buttonsNm);
			tmp = sum(BOI);
			if  tmp ~= 1
				errordlg({['The required Plugin ''3D DENSE Analysis'' CANNOT be solely found:'],[num2str(tmp),' Plugins with the same name!'],'The latest compatible version of the Plugin ''3D DENSE Analysis'' is 0.1.0'},'','modal');
				return
			end
			if strcmpi(get(options(BOI),'Enable'),'off')
				errordlg('Please add Files into the 3D Workspaces!','','modal');
				return			
			end

			self.Handles.hShowMesh = options(BOI);
			set(options(BOI),'Callback',@(s,e)showMesh(handles));
			 %}
			
			% tmp = cellfun(@(x)regexpi(x,'\w*3d*'),handles.hsidebar.TabNames,'UniformOutput',false);
			handles.hsidebar.ActiveTab = handles.hsidebar.NumberOfTabs;
        end

		function deleteSOI(self)
			try deleteHandles(self.hPickSlice); end
			try delete(self.hPickSlice.hPts); end
			try rmfield(self.viewerObj.cache,'regionalStrain'); end
			self.dataObj.Parameterization = [];
			self.dataObj.Strains = [];			
			self.hPickSlice = [];
		end
		
        function delete(self,reset)
            % Make sure all UI components are removed
			self.deleteSOI();
			try self.hShowMesh.delete(); end
			deleteHandles(self.Handles);
			
			% clean MEX Files still open in Matlab
			clear mex;
			% Files solely for the plugin CANNOT be deleted: "SeedPoints"
			try
				rmdir(fullfile(self.InstallDir,'private'),'s');
			catch ERR
				warning(ERR.message);
			end
			
			if ~exist('reset', 'var')
				% Call superclass destructor
				delete@plugins.DENSEanalysisPlugin(self);
			end
        end
		
		function validate(varargin)
            % validate - Check if the plugin can run.
            %
            %   Performs validation to ensure that the state of the program
            %   is correct to be able to run the plugin.
            %
			%   When the user clicks on your "plugin" menu item within the DENSEanalysis GUI and your validate method doesn't produce an error,
			%   the DENSEdata object will automatically be passed to your plugin's run method.
            %
            % USAGE:
            %   DENSE3DPluginLV.validate(data)
            %
            % INPUTS:
            %   data:   Object, DENSEdata object containing all series, images, ROI, and spline information from the DENSEanalysis program.
        end

        function run(varargin)
		%% CANNOT DELETE: otherwise it will become a abstract class
            % run - Method executed when user selects the plugin
            %
            % USAGE:
            %   DENSE3DPlugin4CrescentOrgan.run(data)
            %
            % INPUTS:
            %   data:   Object, DENSEdata object containing all underlying
            %           data from the DENSEanalysis program.
        end

		% function getParentHandles(self)
            % self.hParent = guidata(self.hfig(1));
        % end
		
		function loadFcnREPL(self,flag_new)
		% Copyright (c) 2016 DENSEanalysis Contributors
		% Last Modified: 13:12 July 28, 2017
		% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)
			import plugins.DENSE3D_Plugin_4CrescentOrgan.*
			
            handles = guidata(self.hfig(1));
			% Get directory path by default:
			if flag_new
				type = 'dicom';
				startpath = get(handles.config, 'locations.dicomfolder', userdir());
			else
				type = 'dns';
				startpath = get(handles.config, 'locations.matpath', userdir());
			end

			% try to load new data
			try
				[uipath,uifile] = load(handles.hdata,type,startpath);
			catch ERR
				uipath = [];
				errstr = ERR.message;
				h = errordlg(errstr,'','modal');
				ERR.getReport()
				waitfor(h);
			end
			if isempty(uipath), return; end

			% save path to configure
			if flag_new
				set(handles.config, 'locations.dicomfolder', uipath)
				set(handles.config, 'locations.matfile', '')
				f = 'new';			
			else
				set(handles.config, 'locations.matpath', uipath)
				set(handles.config, 'locations.matfile', uifile)
				[~,f,~] = fileparts(uifile);
				% set(handles.config, 'locations.dnsname', f)
			end

			guidata(handles.hfig,handles);
			% figure name
			set(handles.hfig,'Name',['DENSEanalysis: ' f]);
			% update figure
			resetFcnREPL(handles.hfig);
			

			%% Additional fields
			fields = {'SOI','nSA','nLA'};
			idx = 1:numel(fields);
			if flag_new
				self.dns = handles.hdata.dns;
				
				% create dns for non-DENSE sequences:
				idxSeq = 1:numel(handles.hdata.seq);
				tmp = unique([[handles.hdata.dns.MagIndex],[handles.hdata.dns.PhaIndex]]);
				tmp(isnan(tmp)) = [];
				idxSeq(tmp) = [];
				newDNS = numel(handles.hdata.dns);
				for ii = idxSeq
					newDNS = newDNS+1;
					self.dns(newDNS).UID = dicomuid;
					self.dns(newDNS).MagIndex = repmat(ii,1,numel(self.dns(1).MagIndex));
					self.dns(newDNS).PhaIndex = repmat(ii,1,numel(self.dns(1).PhaIndex));					
					self.dns(newDNS).Number = handles.hdata.seq(ii).CardiacNumberOfImages;
					self.dns(newDNS).PixelSpacing = handles.hdata.seq(ii).PixelSpacing;
					self.dns(newDNS).Type = self.dns(1).Type;
					self.dns(newDNS).Scale = self.dns(1).Scale;
					self.dns(newDNS).EncFreq = self.dns(1).EncFreq;
					self.dns(newDNS).SwapFlag = self.dns(1).SwapFlag;
					self.dns(newDNS).NegFlag = self.dns(1).NegFlag;
				end
				
				% change dns names:
				for ii = 1:numel(self.dns)
					self.dns(ii).Name = handles.hdata.seq(self.dns(ii).MagIndex(1)).DENSEanalysisName;
				end
				
				% reload to make changes take effect:				
				file = fullfile(pwd,'tmp.cache');
				seq = handles.hdata.seq; img = handles.hdata.img; roi = handles.hdata.roi;
				dns = self.dns;
				save(file,'seq','img','dns','roi');
				load(handles.hdata,'dns',file);
				delete(file);
			else
				tmp = {'dns','status'};
				tmp = load(fullfile(uipath,uifile), tmp{:}, '-mat');
				self.dns = tmp.dns;
				
				if isfield(tmp,'status')
					self.status = tmp.status;
					idx = find(~isfield(self.status,fields(:)))';
				end
			end
			
			% load data beyond DENSEdata:
			for ii = idx; self.status.(fields{ii}) = []; end
			
			if ~isfield(self.dns, 'SegPos')
				[self.dns.SegPos] = deal([]);	
				% [self.dns(1:self.status.nSA).SegPos] = deal([]);	
			end
			if ~isfield(self.dns, 'meshCtrl')
				[self.dns.meshCtrl] = deal([]);	
			end
		end

		function reloadFcn(self)
			import plugins.DENSE3D_Plugin_4CrescentOrgan.*
			
            handles = guidata(self.hfig(1));
			if isempty(handles.hdata)
				return
			end
			
			% proper startpath
			startpath = fullfile(get(handles.config, 'locations.matpath', userdir()),get(handles.config, 'locations.matfile', userdir()));

			% try to load new data
			try
				[uipath,uifile] = load(handles.hdata,'dns',startpath);
			catch ERR
				uipath = [];
				errstr = ERR.message;
				h = errordlg(errstr,'','modal');
				ERR.getReport()
				waitfor(h);
			end
			if isempty(uipath), return; end
			
			% update figure
			resetFcnREPL(handles.hfig);
			
			% load beyond DENSEdata:
			fields = {'dns','status'};
			tmp = load(fullfile(uipath,uifile), fields{:}, '-mat');
			self.dns = tmp.dns;
			fields = {'SOI','nSA','nLA'};
			if isfield(tmp,'status')
				self.status = tmp.status;
				idx = find(~isfield(self.status,fields(:)))';
			else
				idx = 1:numel(fields);
			end
			for ii = idx; self.status.(fields{ii}) = []; end

		end

		function saveFcnREPL(self,flag_saveas)
		% Copyright (c) 2016 DENSEanalysis Contributors
		% Last Modified: 13:12 July 28, 2017
		% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)

            handles = guidata(self.hfig(1));
			defaultpath = get(handles.config, 'locations.matpath', '');
			defaultfile = get(handles.config, 'locations.matfile', '');
			file = save(handles.hdata, defaultpath, defaultfile, flag_saveas);
			if isempty(file)
				return
			else
				[p,f,e] = fileparts(file);
				set(handles.config, 'locations.matpath', p);
				set(handles.config, 'locations.matfile', [f, e]);
				set(handles.hfig, 'Name', ['DENSEanalysis: ' f]);
			end
			
			%% save to file beyond DENSEdata:
			tmp = load(file, '-mat'); seq = tmp.seq; img = tmp.img; roi = tmp.roi;
			% if isempty(self.dns)
				% dns = tmp.dns;
			% else
				dns = self.dns;
			% end
			% if isempty(self.status.SOI)
				% save(file,'seq','img','dns','roi');
			% else
				status = self.status;
				save(file,'seq','img','dns','roi','status');
			% end
			
		end

		function menu_runanalysis_REPL(self)
			import plugins.DENSE3D_Plugin_4CrescentOrgan.*

            handles = guidata(self.hfig(1));
			didx = handles.hdense.DENSEIndex;
			ridx = handles.hdense.ROIIndex;
			frame = handles.hdense.Frame;
			% self.hdata = handles.hdata;
			
			% self.spl = handles.hdata.analysis(didx,ridx,frame);
			analysisFcnREPL(self,handles.hdata,didx,ridx,frame);
			if isempty(self.spl), return; end
			%{ 
			if isempty(self.status.nSA)
				for k = 1:numel(handles.hdata.dns)
					tmp = handles.hdata.seq(handles.hdata.dns(k).MagIndex(1)).DENSEanalysisName;
					if ~isempty(strfind(tmp, 'SA'))
						self.status.nSA = self.status.nSA + 1;
					elseif ~isempty(strfind(tmp, 'LA'))
						self.status.nLA = self.status.nLA + 1;
					else
						msgbox(['Invalid Input #',sprintf('%d',k),': ',sprintf('%s',tmp),'. Cannot recognize the type for SA or LA.']);
					end
				end
			end
			 %}

			%{ 
			if mod(self.status.nSA,2) == 0
				if didx == self.status.nSA/2 || didx == self.status.nSA/2+1 
					self.dns(didx).SegPos = handles.hanalysis.straindata.PositionB;
				if ~isempty(self.dns(self.status.nSA/2).SegPos) && ~isempty(self.dns(self.status.nSA/2+1).SegPos)
					handles.hanalysis.straindata.PositionB(1,:) = mean([self.dns(self.status.nSA/2).SegPos(1,:);self.dns(self.status.nSA/2+1).SegPos(1,:)],1);		
					handles.hanalysis.straindata.PositionB(2,:) = mean([self.dns(self.status.nSA/2).SegPos(2,:);self.dns(self.status.nSA/2+1).SegPos(2,:)],1);
					self.dns([1:self.status.nSA/2-1,self.status.nSA/2+2:self.status.nSA]).SegPos = handles.hanalysis.straindata.PositionB;
				end
				end
			elseif didx == round(self.status.nSA/2)
				self.dns(didx).SegPos = handles.hanalysis.straindata.PositionB;
				if didx <= self.status.nSA
					self.dns([1:didx-1,didx+1:self.status.nSA]).SegPos = handles.hanalysis.straindata.PositionB;
				end
			end
			 %}
			%{  
			if isempty(self.status.SOI)
				setSOI(self);
			end
			if any(didx==self.status.SOI)
			 %}
				exportpath = get(handles.config, 'export.mat.location', '');
				if ~isa(handles.hanalysis, 'DataViewer')
					return;
				end	
				% file = handles.hanalysis.exportMat(exportpath);

				% determine strain object
				data = spl2strainFcnREPL(self,didx);%,handles
				if isempty(data)
					return;
				else
					self.straindata = data;
					% save to *.dns
					self.dns(didx).SegPos = data.PositionB;
					self.dns(didx).meshCtrl = data.meshCtrl;		
					
					file = exportMatREPL(self,handles,exportpath);
					% save to file beyond DENSEdata:
					%{ 
					defaultpath = get(handles.config, 'export.mat.location', '');
					handles.hanalysis.exportMat(self, defaultpath);
					defaultpath = get(handles.config, 'export.mat.location', '');
					defaultfile = get(handles.config, 'export.mat.file', '');
					file = fullfile(defaultpath, defaultfile);
					tmp = load(file, '-mat');
					tmp.AnalysisInfo.Nseg = self.straindata.Nseg;
					 %}
				end
				if isempty(file)
				% if isempty(handles.hanalysis.straindata)
					return;
				end
				handles.config.export.mat.location = fileparts(file);

				%{ 
				% self.dns(didx).SegPos = handles.hanalysis.straindata.PositionB;
				self.dns(didx).SegPos = self.straindata.PositionB;
				% self.dns.SegPos_idx(end+1) = didx;
				% self.dns.SegPos(end+1,:) = handles.hanalysis.straindata.PositionB;
				tmp = {self.dns(self.status.SOI).SegPos};
				if ~any(cellfun(@isempty,tmp))
					% tmp1 = []; tmp2 =[];
					% for k = 1:numel(self.status.SOI)
						% tmp1 = [tmp1;tmp{k}(1,:)];
						% tmp2 = [tmp2;tmp{k}(2,:)];
					% end
					% handles.hanalysis.straindata.PositionB(1,:) = mean(tmp1,1);
					% handles.hanalysis.straindata.PositionB(2,:) = mean(tmp2,1);
					self.straindata.PositionB = mean(vertcat(tmp{:}),1);
					% handles.hanalysis.straindata.PositionB = mean(vertcat(tmp{:}),1);
					ind = 1:self.status.nSA; ind(self.status.SOI) = [];
					[self.dns(ind).SegPos] = deal(self.straindata.PositionB);
					% [self.dns(ind).SegPos] = deal(handles.hanalysis.straindata.PositionB);
				else
					msgbox('Next, build the mesh for the second Slice of Interest.');			
				end
			elseif isempty(self.dns(didx).SegPos)
			% elseif ~sum(self.dns.SegPos_idx == didx)
				msgbox('This slice is NOT the Slices of Interest (SOI). First, build the mesh for the SOI. Later, re-build the mesh for this slice again!');
			end
			 %}
			
			% handles.hsidebar.ActiveTab = 3;
		end
		
		function showMesh(self)
		% work in STATIC WORKSPACE: variables in the main function retain
			if isfield(self.hShowMesh,'fig')
				figure(self.hShowMesh.fig);
				return
				%% debug:delete previous figures
                % delete(self.hShowMesh.fig);
				% delete([1:20](ishandle(1:20)));
			end
			
            % if isempty(POI) || isempty(self.dataObj.Data)
				% errordlg('Please add Files into the 3D Workspaces!','','modal');
				% set(self.Handles.hShowMesh,'Value',0);
                % return
            % end
			
			self.hShowMesh.ax = NaN;
			% htitle = NaN;

			self.hShowMesh.meshes = cell(0,1);
			self.hShowMesh.hmeshes = [];

			tmp = [self.dataObj.Data.ROIInfo];
			types = {'SAFull', 'sadual', 'LAFull', 'ladual'};
			self.hShowMesh.isBiv = any(ismember(types, {tmp.ROIType}));
			self.hShowMesh.rois = cat(1, tmp.RestingContour3D);
			self.hShowMesh.hROIs = [];

			self.hShowMesh.rotmat = eye(3);
			self.hShowMesh.center = [0 0 0];
			
			if isempty(self.dataObj.EpicardialMesh)
				% Perform analysis
				self.generateMeshes();
			else
				button = questdlg({'Previously Defined Meshes Detected!','Choose what you want to do with it'},'Warning','Delete','Display','Delete');
				if strcmpi(button,'Delete')
					self.dataObj.reset();
					self.generateMeshes();
				end
			end
			if isempty(self.dataObj.EndocardialMesh); return; end
						
			self.hShowMesh.fig = figure(...
			'Name', 'Refine 3D Mesh | created by Liu', ...
			'NumberTitle',   'off', ...
			'color', self.BackgroundColor, ...%'none'
			'units', 'normalized',...
			'outerposition',[0 0 1 1],...
			'renderer', 'opengl',...
			'WindowStyle', 'docked',...
			'CloseRequestFcn',@(s,e)deleteFcn(), ...
			'KeyPressFcn', @windowkeypress, ...
			'MenuBar',       'none', ...
			'Toolbar',       'none');
			% 'IntegerHandle',   'off', ...
			% 'Position',      pos, ...
			% 'DockControls',      'off', ...% Not for 'docked' WindowStyle
			% 'WindowKeyPressFcn', @windowkeypress, ...
						

			% light color in black BG:
			% color = self.viewerObj.dispclr;
			
			%% Toolbar setup:				
			tb = uitoolbar(self.hShowMesh.fig);
			% uitoolfactory(tb,'Standard.EditPlot');
			uitoolfactory(tb,'Standard.SaveFigure');
			uitoolfactory(tb,'Standard.PrintFigure');
			% Separator in front of a tool
			tmp = uitoolfactory(tb,'Exploration.ZoomIn');
			set(tmp,'Separator','on','TooltipString','Zoom In(=)');
			tmp = uitoolfactory(tb,'Exploration.ZoomOut');
			set(tmp,'TooltipString','Zoom Out(-)');
			tmp = uitoolfactory(tb,'Exploration.Rotate');
			set(tmp,'TooltipString','Rotate(Alt+R)');
			tmp = uitoolfactory(tb,'Exploration.Pan');
			set(tmp,'TooltipString','Pan(Alt+P)');
			tmp = uitoolfactory(tb,'Exploration.DataCursor');
			set(tmp,'TooltipString','Data Cursor(Alt+D)');
			tmp = uitoolfactory(tb,'Exploration.Brushing');
			set(tmp,'TooltipString','Highlight(Alt+B)');
			
			self.hShowMesh.ax = axes('Parent', self.hShowMesh.fig,...
			'Units', 'normalized', ...
			'FontWeight', 'bold', ...
			'Position', [0 0 1 1], ...
			'Visible', 'off');
			% 'xtick', [], ...
			% 'ytick', [], ...
			% 'Clipping', 'off', ...
			% 'box',  'on', ...
			% 'color', 'none', ...
			% 'xcolor', color, ...
			% 'zcolor', color, ...
			% 'ycolor', color, ...
			
			%% Create multiple images windows:
			%{ 
			self.hShowMesh.ax(1) = axes('Units', 'normalized', 'color', 'k', 'Position', [0 0.5 .5 .5]);
			self.hShowMesh.ax(2) = axes('Units', 'normalized', 'color', 'k', 'Position', [0.5 0.5 0.5 0.5]);
			self.hShowMesh.ax(3) = axes('Units', 'normalized', 'color', 'k', 'Position', [0 0 .5 .5]);
			self.hShowMesh.ax(4) = axes('Units', 'normalized', 'color', 'k', 'Position', [0.5 0 .5 .5]);
			 %}

			% htitle = title(self.hShowMesh.ax, '3D Mesh', 'Color', color);
			% hlegend = legend(self.hShowMesh.ax, {'3D Mesh'}, 'AutoUpdate', 'off', 'TextColor', color, 'EdgeColor', color);
			% title(hlegend,'3D Mesh')

			self.hShowMesh.center = mean(self.dataObj.EpicardialMesh.vertices, 1);
			self.hShowMesh.rotmat = self.dataObj.rotationMatrix();

			self.hShowMesh.meshes = {[self.dataObj.EpicardialMesh, self.dataObj.EndocardialMesh]};
			self.hShowMesh.meshes{1} = self.dataObj.rotateMesh(self.hShowMesh.meshes{1}, self.hShowMesh.rotmat, self.hShowMesh.center);

			self.hShowMesh.ROIs{1} = cellfun(@(x)bsxfun(@plus, bsxfun(@minus, x, self.hShowMesh.center) * self.hShowMesh.rotmat, self.hShowMesh.center), self.hShowMesh.rois,'UniformOutput', false);
			%{ 
			self.hShowMesh.rois = [self.dataObj.Data.ROIInfo];
			self.hShowMesh.rois = cat(1, self.hShowMesh.rois.RestingContour3D);
			for k = 1:numel(self.hShowMesh.meshes{1})
				for slice=1:size(self.hShowMesh.rois,1)
					 self.hShowMesh.rois{slice,k} = bsxfun(@plus, bsxfun(@minus, self.hShowMesh.rois{slice,k}, self.hShowMesh.center) * self.hShowMesh.rotmat, self.hShowMesh.center);
				end
			end
			self.hShowMesh.ROIs{1} = self.hShowMesh.rois;
			 %}
			
			names = {'Epi','LVendo', 'RVendo'};
			for k = 1:numel(self.hShowMesh.meshes{1})
				self.hShowMesh.hmeshes(k) = patch(self.hShowMesh.meshes{1}(k), 'EdgeColor', 'k', 'FaceColor', 'w', 'DisplayName', [names{k},' Mesh'], 'FaceAlpha', 1, 'Parent', self.hShowMesh.ax);
				hold(self.hShowMesh.ax, 'on')
				for slice=1:size(self.hShowMesh.ROIs{1},1)
					self.hShowMesh.hROIs(slice,k) = plot3(self.hShowMesh.ax,self.hShowMesh.ROIs{1}{slice,k}(:,1),self.hShowMesh.ROIs{1}{slice,k}(:,2),self.hShowMesh.ROIs{1}{slice,k}(:,3),'DisplayName', [names{k},' Contours'],'LineWidth',1,'color', 'y');
				end
			end
			% Lighting
			set(self.hShowMesh.hmeshes, 'FaceLighting', 'none', 'EdgeLighting', 'none', 'BackFaceLighting', 'unlit');
% Test figure properites:
%{ 
figure
patch(self.hShowMesh.meshes{1}(1), 'FaceColor', 'w', 'FaceAlpha', 1, 'EdgeColor', 'r');
hold on; axis equal; xlabel('x'); ylabel('y'); zlabel('z'); grid on;
patch(self.hShowMesh.meshes{1}(3), 'FaceColor', 'w', 'FaceAlpha', 1, 'EdgeColor', 'b');
 %}
% Debug RV Mesh in the ref. frame:
%{ 
for k = 1:numel(self.hShowMesh.meshes{1})
	figure;
	axis equal; xlabel('x'); ylabel('y'); zlabel('z'); grid on;
	patch(self.hShowMesh.meshes{1}(k), 'FaceColor', 'w', 'FaceAlpha', 0.5);
end
 %}

			 
			%% playbar & listeners
			frames = self.viewerObj.Frames;
			self.hShowMesh.hplaybar = playbar(self.hShowMesh.fig);
			% place the playbar at the bottom/self.hShowMesh.center of the display panel
			self.hShowMesh.hplaybar.Units = 'normalized';
			self.hShowMesh.hplaybar.Position = [0.40 0.01 0.15 0.04];
			
			self.hShowMesh.hplaybar.Min = frames(1)-1;
			self.hShowMesh.hplaybar.Max = frames(2);
			self.hShowMesh.hplaybar.Enable = 'on';
			self.hShowMesh.hplaybar.Value = frames(1)-1;

			% Most Anonymous Functions have two input arguments by default:
			self.hShowMesh.hlisten_playbar = addlistener(self.hShowMesh.hplaybar,'NewValue',@(s,e)playbackFcn());%@(varargin)playbackFcn(varargin{1}.Value)
			self.hShowMesh.hlisten_playbar.Enabled = true;

			self.hShowMesh.delete = @deleteFcn;
			
			grid(self.hShowMesh.ax, 'on');
			% arrayfun(@(x)grid(x, 'on'),self.hShowMesh.ax)
			axis(self.hShowMesh.ax, 'equal', 'tight', 'manual');
			view(self.hShowMesh.ax, 3);
			% Set base zoom level
			zoom(self.hShowMesh.ax, 'reset');
			% Show or hide Plotted objects
			plotbrowser(self.hShowMesh.fig,'on');
			figure(self.hShowMesh.fig);
			% jFig = get(self.hShowMesh.fig,'JavaFrame');
			% jFig.getAxisComponent.requestFocus;
			plotedit(self.hShowMesh.fig, 'off');

			function deleteFcn()
				try delete(self.hShowMesh.fig(ishandle(self.hShowMesh.fig))); end
				set(self.Handles.hShowMesh,'Value',0);
				self.hShowMesh = [];
				
				try rmfield(self.viewerObj.cache,'regionalStrain'); end
				self.hPickSlice = [];
				set(self.Handles.hPickSlice,'Value',0);
				
				try delete(findall(0,'type','figure','-and','tag','WaitbarTimer')); end
			end

			function playbackFcn()
			%% Copyright (c) of the following section of the codes: 2016 DENSEanalysis Contributors
			% The codes of the following section are from the function "playbackFcn" in "DENSE3D.m".
				% Image Frames:
				self.hShowMesh.frame = self.hShowMesh.hplaybar.Value;
				% Referential self.hShowMesh.Frame + Image Frames:
				self.hShowMesh.Frame = self.hShowMesh.frame + 1;

				if numel(self.hShowMesh.meshes)<self.hShowMesh.Frame || isempty(self.hShowMesh.ROIs{self.hShowMesh.Frame})

					if isempty(self.dataObj.Interpolants)
						% Compute displacement splines as needed
						self.viewerObj.computeDisplacementSplines()
					end

					interp = self.dataObj.Interpolants(self.hShowMesh.frame);

					% Interpolate the epicardial surface
					M = [self.dataObj.EpicardialMesh, self.dataObj.EndocardialMesh];

					for k = 1:numel(M)
						M(k).vertices = M(k).vertices + interp.query(M(k).vertices);
					end

					self.hShowMesh.meshes{self.hShowMesh.Frame} = self.dataObj.rotateMesh(M, self.hShowMesh.rotmat, self.hShowMesh.center);
					
			%% Copyright (c) of the following section of the codes: Zhanqiu Liu (lafeir.lew@gmail.com)
			% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)
			% Last Modified: 19:32 June 27, 2018
			% Debug Meshes in the current frame:
					%% Drawn roi:
					tmp = [self.dataObj.Data.ROIInfo];
					self.hShowMesh.ROIs{self.hShowMesh.Frame} = {};
					for k = 1:numel(tmp)
						self.hShowMesh.ROIs{self.hShowMesh.Frame} = [self.hShowMesh.ROIs{self.hShowMesh.Frame}; tmp(k).Contour3D(self.hShowMesh.frame,:)];
					end
					self.hShowMesh.ROIs{self.hShowMesh.Frame} = cellfun(@(x)bsxfun(@plus, bsxfun(@minus, x, self.hShowMesh.center) * self.hShowMesh.rotmat, self.hShowMesh.center), self.hShowMesh.ROIs{self.hShowMesh.Frame},'UniformOutput', false);
					%% Calculated roi:
					%{ 
					for k = 1:numel(self.hShowMesh.meshes{self.hShowMesh.Frame})
						% self.hShowMesh.meshes{self.hShowMesh.Frame}:{epi,endoLV,endoRV}
						for slice=1:size(self.hShowMesh.rois,1)
							self.hShowMesh.ROIs{self.hShowMesh.Frame}{slice,k} = self.hShowMesh.rois{slice,k} + interp.query(self.hShowMesh.rois{slice,k});
							self.hShowMesh.ROIs{self.hShowMesh.Frame}{slice,k} = bsxfun(@plus, bsxfun(@minus, self.hShowMesh.ROIs{self.hShowMesh.Frame}{slice,k}, self.hShowMesh.center) * self.hShowMesh.rotmat, self.hShowMesh.center);
						end
					end
					 %}
				end
									
				for k = 1:numel(self.hShowMesh.meshes{self.hShowMesh.Frame})
				% self.hShowMesh.meshes{self.hShowMesh.Frame}:{epi,endoLV,endoRV}
					set(self.hShowMesh.hmeshes(k), self.hShowMesh.meshes{self.hShowMesh.Frame}(k))
					for slice=1:size(self.hShowMesh.ROIs{self.hShowMesh.Frame},1)
						set(self.hShowMesh.hROIs(slice,k), 'XData', self.hShowMesh.ROIs{self.hShowMesh.Frame}{slice,k}(:,1), 'YData', self.hShowMesh.ROIs{self.hShowMesh.Frame}{slice,k}(:,2), 'ZData', self.hShowMesh.ROIs{self.hShowMesh.Frame}{slice,k}(:,3));
					end
				end
				
				drawnow
				
				try
					if numel(self.hPickSlice.strainLocs)<self.hShowMesh.Frame || isempty(self.hPickSlice.strainLocs{self.hShowMesh.Frame})
						self.hPickSlice.strainLocs{self.hShowMesh.Frame} = self.hPickSlice.transStrainLocs();
					end
					self.hPickSlice.findNoises();
				end
			end
			
			function windowkeypress(~,evnt)
			% Last Modified: 8:37 PM Friday, November 20, 2015
			% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)
				if isempty(evnt.Modifier)
					key = evnt.Key;
				else
					key = [evnt.Modifier{1},evnt.Key];
					% modifiers = sort(evnt.Modifier);
					% key = strcat(sprintf('%s-', modifiers{:}), evnt.Key);
				end
				
				% avoid hplaybar being clear:
				playbar = self.hShowMesh.hplaybar;
				
				switch key
				% toolbar
					case 'equal'
						zoom(self.hShowMesh.ax, 2);
					case 'hyphen'
						zoom(self.hShowMesh.ax, 0.5);
					case 'altp'
					% Toggle for panning images
						% pan(self.hShowMesh.fig);
						% msgbox ok;
						pan(self.hShowMesh.fig);
						hmanager = uigetmodemanager(self.hShowMesh.fig);
						set(hmanager.WindowListenerHandles,'enable','off');
						set(self.hShowMesh.fig,'WindowKeyReleaseFcn',@(s,e)pan(self.hShowMesh.fig));
					case 'altd'
					% Toggle for data cursor
						datacursormode(self.hShowMesh.fig);
						hmanager = uigetmodemanager(self.hShowMesh.fig);
						set(hmanager.WindowListenerHandles,'enable','off');
						set(self.hShowMesh.fig,'WindowKeyReleaseFcn',@(s,e)pan(self.hShowMesh.fig));
					case 'altr'
					% Toggle for 3D rotation tool
						rotate3d(self.hShowMesh.fig);
						hmanager = uigetmodemanager(self.hShowMesh.fig);
						set(hmanager.WindowListenerHandles,'enable','off');
						set(self.hShowMesh.fig,'WindowKeyReleaseFcn',@(s,e)rotate3d(self.hShowMesh.fig));
					case 'altb'
					% Toggle for brush tool
						brush(self.hShowMesh.fig);
						hmanager = uigetmodemanager(self.hShowMesh.fig);
						set(hmanager.WindowListenerHandles,'enable','off');
						set(self.hShowMesh.fig,'WindowKeyReleaseFcn',@(s,e)brush(self.hShowMesh.fig));
					% 'command' for Macintosh computers
					case {'controls','command-s'}
						filemenufcn(self.hShowMesh.fig,'FileSave');
					case {'controlp','command-p'}
						printdlg(self.hShowMesh.fig);
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
		end

		function generateMeshes(self)
			% Generate the surface meshes for the endo- and epicardial contours
			if self.hShowMesh.isBiv
				self.genRVendoMesh();
				% self.dataObj.EndocardialMesh(2) = surfacemesh(cat(1, self.hShowMesh.rois{:,3}), 1); % adjust = true;
				if isempty(self.dataObj.EndocardialMesh); return; end
			end
			self.dataObj.EpicardialMesh = surfacemesh(cat(1, self.hShowMesh.rois{:,1}));
			tmp = surfacemesh(cat(1, self.hShowMesh.rois{:,2}));
			try
				self.dataObj.EndocardialMesh = [tmp,self.dataObj.EndocardialMesh];
			catch
				self.dataObj.EndocardialMesh = [struct('faces', tmp.faces, 'vertices', tmp.vertices), struct('faces', self.dataObj.EndocardialMesh.faces, 'vertices',self.dataObj.EndocardialMesh.vertices)];
			end

			self.dataObj.Apex = self.dataObj.autoApex();

			% Remove the mesh vertices that are above the basal slice
			self.dataObj.pruneMesh();
			%{ 
			figure;
			hold on; axis equal; xlabel('x'); ylabel('y'); zlabel('z'); grid on;
			patch(self.dataObj.EpicardialMesh, 'FaceColor', 'w', 'FaceAlpha', 0.5);
			patch(self.dataObj.EndocardialMesh(1), 'FaceColor', 'w', 'FaceAlpha', 0.5);
			patch(self.dataObj.EndocardialMesh(2), 'FaceColor', 'w', 'FaceAlpha', 0.5);
			 %}
			 
		end

		function genRVendoMesh(self)
		% Small animals (rodent, ect) usually have a much more crescent-shaped RV cavity
		% Handpick a good seed point where a small sphere will be placed
		% the location of the small sphere is not obtained by just taking the mean of the supplied contour points.
		%
		% Last Modified: 1:06 PM Thursday, October 15, 2015
		% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)
			
			import plugins.DENSE3D_Plugin_4CrescentOrgan.*
			% Use the up-to-dated Configuration file:
			% settingfile = Configuration(fullfile(fileparts(which(class(self.dataObj))), 'settings.json'));
			% Use the pre-loaded Configuration file:
			[center,initial_nor_R,final_nor_R,initial_tan_R,final_tan_R,~,~,~,~,~,fname] = SeedPoints(self.configObj);%settingfile
			
			fExit = false;
			%% pick a SeedPt
			if isempty(center)
				hFig = figure(...
				'Name', 'Handpick SeedPoint', ...
				'NumberTitle',   'off', ...
				'color', self.BackgroundColor, ...
				'units', 'normalized',...
				'outerposition',[0 0 1 1],...
				'renderer', 'opengl',...
				'DockControls', 'off',...
				'CloseRequestFcn',@(s,e)deleteFig(), ...
				'MenuBar',       'none', ...
				'Toolbar',       'figure');
				% 'WindowKeyPressFcn', @windowkeypress, ...
				hold on; axis equal; xlabel('x'); ylabel('y'); zlabel('z'); grid on;
				center = mean(cat(1, self.hShowMesh.rois{:,3}), 1);
				h = impoint(gca,center(1:2));
				hpt = scatter3(center(:,1),center(:,2),center(:,3),25,'b','x','LineWidth',2);
				
				for slc=1:size(self.hShowMesh.rois,1)
					plot3(self.hShowMesh.rois{slc,3}(:,1),self.hShowMesh.rois{slc,3}(:,2),self.hShowMesh.rois{slc,3}(:,3),'LineWidth',1.5,'color', 'b');
				end
				for slc=1:size(self.hShowMesh.rois,1)
					plot3(self.hShowMesh.rois{slc,2}(:,1),self.hShowMesh.rois{slc,2}(:,2),self.hShowMesh.rois{slc,2}(:,3),'LineWidth',1.5,'color', 'r');
				end
				
				% consists of two orthogonal vectors defining the imaging plane in real-world coordinates
				IOP = zeros(1,0); slc = 1;
				while numel(IOP) ~= 6
					IOP = self.dataObj.data(slc).SequenceInfo(1, 1).ImageOrientationPatient;
					slc = slc + 1;
				end
				nor = cross(IOP(1:3),IOP(4:6));
				tmp = 1/norm(nor); nor = nor * tmp;% sum(nor.^2)==1
				nor = reshape(nor, 1, []);
				view(nor);
				
				
				button = uicontrol('Style','Pushbutton','String','Confirm Placement!','units','Normalized','Position',[0.01 0.01 0.1 0.05],'Callback',@(s,e)adjustMesh()) ;
				
				annotation('textbox',[0.7 .85 0.2 0.15],...
                'Units', 'normalized',...
				'EdgeColor', 'none',...
				'FontSize', 8,...
				'String',{'Imported RestingContour3D:','* If RV Endocardial Surface is NOT longitudinally continued:','select different Seedpoints for different slice when running Phase Analysis','* DRAG the BLACK CROSS-CIRCLE POINT: place the Surface-Mesh-Fitting balloon','-Black cross point: projected point in the centroid','-Click Button: confirm placement'});
				% 'FontName','Arial',...
				% 'FitBoxToText', 'on',...
				% 'LineStyle','-',...
				% 'LineWidth',2,...
				% 'BackgroundColor',self.BackgroundColor,...
				% 'Color',[0.84 0.16 0]);
				% legend({'Imported RestingContour3D:','If RV Endocardial Surface is NOT longitudinally continued:','select different Seedpoints for different slice when running Phase Analysis','DRAG the BLACK CROSS-CIRCLE POINT: place the Surface-Mesh-Fitting balloon','Double-click Mouse Middle Button: confirm placement','Black cross point: projected point in the centroid'});%,'AutoUpdate','off'); %Starting in R2017a
				
				% Enforce boundary constraint
				fcn = makeConstrainToRectFcn('impoint',get(gca,'XLim'),get(gca,'YLim'));
				setPositionConstraintFcn(h,fcn);
				setColor(h,'k');
				% Update the projected Seedpoint with the new selection
				addNewPositionCallback(h,@(pos) proj(pos));

				set(self.Handles.hShowMesh,'Value',1);
				% Click Button to resume execution of the MATLAB command line
				waitfor(hFig, 'Name');
			end
			
			function CTR = proj(pos)
				SctPt = [pos 0];
				CTR = SctPt - dot(SctPt - center, nor) * nor;
				set(hpt,'XData',CTR(:,1),'YData',CTR(:,2),'ZData',CTR(:,3));
				title(sprintf('Balloon Location:(%1.0f,%1.0f,%1.0f)',CTR(:,1),CTR(:,2),CTR(:,3)));
				drawnow
			end
			
			function adjustMesh()
				position = getPosition(h);
				if isempty(position)
					keyboard
				else
					center = proj(position)
					set(hFig,'Name',sprintf('Seedpoint (%1.0f,%1.0f,%1.0f) Selected for Dataset %s',center(:,1),center(:,2),center(:,3),fname));
					delete(h); delete(button);
					drawnow
				end
			end

			function deleteFig()
				fExit = true;
				delete(hFig);
			end
			
			if fExit; return; end
			% Initialize waitbar timer
			hwait = waitbartimer();
			% Delete hwait ASA cleanupObj deleted: too costly
			% cleanupObj = onCleanup(@(x)delete(hwait(isvalid(hwait))));
			hwait.String = 'Rendering 3D Mesh...';
			hwait.WindowStyle = 'modal';
			hwait.AllowClose = false;
            hwait.Visible = 'on';
			hwait.start;
			drawnow
			opts.WindowStyle = 'normal';

			%% Initialize a RV endo mesh:
			points = cat(1, self.hShowMesh.rois{:,3});
			mesh = Generate_Mesh(points, 0, center);% adjust = true;

			%%% Adjust Mesh controlled by the Fitting parameters:
			% set the starting and ending parameter values that are lineraly varied (from the starting to the ending values) over the iterations. The parameters are "normal R", "tagent R", "smoothing weight", "points weight" and the number of interations.
			initial_smoothness_weight = .9;
			final_smoothness_weight = .9;
			initial_points_weight = .1;
			final_points_weight = .2;
			iterations = 1000;
			%% for human: 
			% initial_nor_R = 10; % [mm]
			% final_nor_R = 5; % [mm]
			% initial_tan_R = 5; % [mm]
			% final_tan_R = 5; % [mm]
			
			meshes.faces = [mesh.tri_n1;mesh.tri_n2;mesh.tri_n3]';

		if isempty(initial_nor_R)
			button = questdlg('Do you know the BEST VALUES for the fitting parameters: Initial Normal R & Final Normal R?','Iterative Calculation','No');
			if strcmpi(button,'no')
				while true
					% wait just as a waitbar timer doese
					tmp = inputdlg('Input a positive LOWER LIMIT for the distance from the mesh in the TANGENT direction:','Initial Tangent R',1,{'0.2'});
					if isempty(tmp); continue; end;
					tmp = sscanf(sprintf('%s*', tmp{:}), '%f*');
					if tmp > 0
						initial_tan_R = tmp
						break
					end		
				end
				while true
					tmp = inputdlg('Input a UPPER LIMIT for the distance from the mesh in the tangent direction:','Final TANGENT R',1,{'0.5'});
					if isempty(tmp); continue; end;
					tmp = sscanf(sprintf('%s*', tmp{:}), '%f*');
					if tmp >= initial_tan_R
						final_tan_R = tmp
						break
					end		
				end
				for final_nor_R = 0.1:0.1:1
					for initial_nor_R = round(3*final_nor_R)/10:0.1:round(7*final_nor_R)/10
						%% Debug Meshes in the current frame:
						meshes.vertices = Adjust_Mesh_Mex(mesh, points, iterations, [initial_smoothness_weight, final_smoothness_weight], [initial_points_weight, final_points_weight], [initial_nor_R, final_nor_R],  [initial_tan_R, final_tan_R]);

						figure;
						hold on; axis equal; xlabel('x'); ylabel('y'); zlabel('z'); grid on;
						title(['initial_nor_R=',num2str(initial_nor_R),';','final_nor_R=',num2str(final_nor_R),';','initial_tan_R=',num2str(initial_tan_R),';','final_tan_R=',num2str(final_tan_R)]);
						patch(meshes, 'FaceColor', 'w', 'FaceAlpha', 0.5);
						% scatter3(points(:,1),points(:,2),points(:,3),3,'k','*');
						for slice=1:size(self.hShowMesh.rois,1)
							plot3(self.hShowMesh.rois{slice,3}(:,1),self.hShowMesh.rois{slice,3}(:,2),self.hShowMesh.rois{slice,3}(:,3),'LineWidth',1,'color', 'b');
						end
					end
				end
				
				h = msgbox({'Compare a series of fittings of generated RV endo meshes to the actual RV endo contours','Write down the best fitting parameters','Click any button when you are ready!'},'Find the best fitting');
				waitfor(h);
			end
			
			while true
				tmp = inputdlg('Choose the best LOWER LIMIT for the distance from the mesh in the NORMAL direction based on a series of fittings of generated RV endo meshes to the actual RV endo contours:','Initial Normal R',1,{''},opts);
				if isempty(tmp); continue; end;
				tmp = sscanf(sprintf('%s*', tmp{:}), '%f*');
				if tmp > 0
					initial_nor_R = tmp
					break
				end		
			end
			while true
				tmp = inputdlg('Choose the best UPPER LIMIT for the distance from the mesh in the NORMAL direction based on a series of fittings of generated RV endo meshes to the actual RV endo contours:','Final Normal R',1,{''},opts);
				if isempty(tmp); continue; end;
				tmp = sscanf(sprintf('%s*', tmp{:}), '%f*');
				if tmp >= initial_nor_R
					final_nor_R = tmp
					break
				end		
			end
		end
				
		if isempty(initial_tan_R)
			button = questdlg('Do you know the BEST VALUES for the fitting parameters: Initial Tangent R & Final Tangent R?','Iterative Calculation','No');
			if strcmpi(button,'no')
				for final_tan_R = 0.1:0.1:1
					for initial_tan_R = 0.1:0.1:final_tan_R
						meshes.vertices = Adjust_Mesh_Mex(mesh, points, iterations, [initial_smoothness_weight, final_smoothness_weight], [initial_points_weight, final_points_weight], [initial_nor_R, final_nor_R],  [initial_tan_R, final_tan_R]);

						figure;
						hold on; axis equal; xlabel('x'); ylabel('y'); zlabel('z'); grid on;
						title(['initial_nor_R=',num2str(initial_nor_R),';','final_nor_R=',num2str(final_nor_R),';','initial_tan_R=',num2str(initial_tan_R),';','final_tan_R=',num2str(final_tan_R)]);
						patch(meshes, 'FaceColor', 'w', 'FaceAlpha', 0.5);
						% scatter3(points(:,1),points(:,2),points(:,3),3,'k','*');
						for slice=1:size(self.hShowMesh.rois,1)
							plot3(self.hShowMesh.rois{slice,3}(:,1),self.hShowMesh.rois{slice,3}(:,2),self.hShowMesh.rois{slice,3}(:,3),'LineWidth',1,'color', 'b');
						end
					end
				end

				h = msgbox({'Compare a series of fittings of generated RV endo meshes to the actual RV endo contours','Write down the best fitting parameters','Click any button when you are ready!'},'Find the best fitting');
				% wait until the object closes
				waitfor(h);
			end
			
			while true
				tmp = inputdlg('Choose the best LOWER LIMIT for the distance from the mesh in the TANGENT direction based on a series of fittings of generated RV endo meshes to the actual RV endo contours:','Initial Tangent R',1,{''},opts);
				if isempty(tmp); continue; end;
				tmp = sscanf(sprintf('%s*', tmp{:}), '%f*');
				if tmp > 0
					initial_tan_R = tmp
					break
				end		
			end
			while true
				tmp = inputdlg('Choose the best UPPER LIMIT for the distance from the mesh in the TANGENT direction based on a series of fittings of generated RV endo meshes to the actual RV endo contours:','Final Tangent R',1,{''},opts);
				if isempty(tmp); continue; end;
				tmp = sscanf(sprintf('%s*', tmp{:}), '%f*');
				if tmp >= initial_tan_R
					final_tan_R = tmp
					break
				end		
			end
		end
			
			meshes.vertices = Adjust_Mesh_Mex(mesh, points, iterations, [initial_smoothness_weight, final_smoothness_weight], [initial_points_weight, final_points_weight], [initial_nor_R, final_nor_R],  [initial_tan_R, final_tan_R]);
			
			import plugins.DENSE3D_Plugin_4CrescentOrgan.*
			self.dataObj.EndocardialMesh = struct('faces', [], 'vertices', []);
			% try
				self.dataObj.EndocardialMesh = wrapMesh(meshes);
			% catch
				% self.dataObj.EndocardialMesh(2) = struct('faces', meshes.faces, 'vertices', meshes.vertices);
			% end
			
			%% remove waitbar timer
			if exist('hwait','var'); hwait.stop; delete(hwait); end%delete(findall(0,'type','figure','-and','tag','WaitbarTimer'));
		end
		
		function pickSlice(self)
		% Last Modified: 23:03 June 29, 2018
		% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)
			
			if isempty(self.hShowMesh)
				self.showMesh();
			end
			
			import plugins.DENSE3D_Plugin_4CrescentOrgan.*
			% try deleteHandles(self.hPickSlice); end; try delete(self.hPickSlice.hPts); end

			if ~isempty(self.hPickSlice)
				set(self.Handles.hPickSlice,'Value',1);
				button = questdlg({'Previously Defined Slice of Interest Detected!','Choose what you want to do with it'},'Warning','Delete','Do nothing','Delete');
				if strcmpi(button,'Delete')
					self.deleteSOI();
				else
					figure(self.hShowMesh.fig);
					return
				end
			end
			
            % if isempty(self.dataObj.Data)
				% errordlg('Please add Files into the 3D Workspaces!','','modal');
				% Style='radiobutton': 'Selected'='off' means the object being selected
				% set(self.Handles.hPickSlice,'Value',0);
				% return			
			% end
			

			hwait = waitbartimer();
			hwait.WindowStyle = 'modal';
			hwait.AllowClose = false;
            hwait.Visible = 'on';
			hwait.start;
			drawnow
			
            if isempty(self.dataObj.Interpolants)
				hwait.String = 'Fitting Radial Basis Functions...';
                self.dataObj.interpolate();
            end

            hwait.String = 'Parameterizing the ventricle...';
			self.parameterize();

            % Assign normalized parameters corresponding to each vertice of EndocardialMesh
			% self.longitudinalParameterize();
			
			% Compute the strains
            hwait.String = 'Computing Strains...';
			points = self.dataObj.samplePoints(self.hPickSlice.length/30);
			computeStrains(self, points);
% Debug:
%{
self.dataObj.preview();
hmeshes = findall(allchild(gca),'-property','FaceColor');
set(hmeshes,'FaceColor', 'w');
set(0, 'CurrentFigure', 2);
hold on; axis equal; xlabel('x'); ylabel('y'); zlabel('z'); grid on;
scatter3(points(:,1),points(:,2),points(:,3),5,'k','o');
%}
						
		
            hwait.String = 'Setting up the Viewer...';

			% CALLBACK response to user interaction: push button clicks, slider movements, mouse button up, or Click another component/Enter
			[~,~,~,~,~,self.hShowMesh.frame,slice,ant,inf,tol,fname,parentDir] = SeedPoints(self.configObj);
			if isempty(slice); slice=[.5, .167]; ant=0.01; inf=0.01; tol=0; self.hShowMesh.frame=round(mean(self.viewerObj.Frames)); end
			% Avoid creating new object of playbar when callback: 
			tmp= self.hShowMesh.hplaybar;
			tmp.Value = self.hShowMesh.frame;
			
			self.hPickSlice.hText(1) = annotation(self.hShowMesh.fig,...
				'textbox', [0.65 0.97 0.16 0.02],...
                'Units', 'normalized', ...
				'EdgeColor', 'none',...
				'FontSize', 8,...
				'HorizontalAlignment', 'center',...
				'VerticalAlignment', 'middle',...
				'String',  'Location of the Slice');
			self.hPickSlice.hSlider(1) = uicontrol(...
				'parent', self.hShowMesh.fig,...
				'KeyPressFcn', @keypress,...
				'style', 'slider',...
                'Units', 'normalized', ...
				'Position', [0.65 0.95 0.12 0.02],...
                'Value', slice(1), ...
				'callback', @(s,e)sliderCB(1));
			self.hPickSlice.hEdit(1) = uicontrol(...
				'parent', self.hShowMesh.fig,...
				'KeyPressFcn', @keypress,...
				'style', 'edit',...
                'Units', 'normalized', ...
				'Position', [0.77 0.95 0.04 0.02],...
				'String', num2str(slice(1)), ...
				'callback', @(s,e)textCB(1));
			self.hPickSlice.hText(end+1) = annotation(self.hShowMesh.fig,...
				'textbox', [0.82 0.97 0.16 0.02],...
                'Units', 'normalized', ...
				'EdgeColor', 'none',...
				'FontSize', 8,...
				'HorizontalAlignment', 'center',...
				'VerticalAlignment', 'middle',...
				'String',  'Half Thickness of the Slice');
			self.hPickSlice.hSlider(2) = uicontrol(...
				'parent', self.hShowMesh.fig,...
				'KeyPressFcn', @keypress,...
				'style', 'slider',...
                'Units', 'normalized', ...
				'Position', [0.82 0.95 0.12 0.02],...
                'Value', slice(2), ...
				'callback', @(s,e)sliderCB(2));
                % 'Min', .01, ...
                % 'SliderStep', [0.01 0.10], ...%unit: fraction
			self.hPickSlice.hEdit(2) = uicontrol(...
				'parent', self.hShowMesh.fig,...
				'KeyPressFcn', @keypress,...
				'style', 'edit',...
                'Units', 'normalized', ...
				'Position', [0.94 0.95 0.04 0.02],...
				'String', num2str(slice(2)), ...
				'callback', @(s,e)textCB(2));
			
			self.hPickSlice.hButton = uicontrol(...
				'parent', self.hShowMesh.fig,...
				'KeyPressFcn', @keypress,...
				'Style', 'Pushbutton',...
				'String', 'Confirm New Slice of Interest!',...
				'units', 'Normalized',...
				'Position', [0.01 0.01 0.1 0.07],...
				'Callback', @(s,e)regionalStrains());
			self.hPickSlice.hText(end+1) = annotation(self.hShowMesh.fig,...
				'textbox', [0.01 .85 0.2 0.15],...
                'Units', 'normalized',...
				'EdgeColor', 'none',...
				'FontSize', 8,...
				'String', {'NOTICE:','1.Adjust the slice based on the shape of highlighted RV endocardial elements.','The lower limit of the slice should be large than 0 and the upper limit should be less than 1 !','2.Endo-surface of RV apex with crescent-shaped SA contour is hard to be reconsturcted:','Suggest to use LV endo-mesh to find cut-off longitudinal locations for dividing RV into equal thirds perpendicular to its long axis!','3.The highlighted epicardial surface is larger than the actual region used for strain analysis!','4.Deformations of Meshes calculated from DENSE phase data VS Drawn Contours','5.Frame #0 is the Referential Frame.','6.All shown numbers are normalized in a range of [0, 1].'});
			
			self.hPickSlice.hText(end+1) = annotation(self.hShowMesh.fig,...
				'textbox', [0.82 0.90 0.07 0.04],...
				'Units', 'normalized', ...
				'EdgeColor', 'none',...
				'FontSize', 8,...
				'VerticalAlignment', 'middle',...
				'String',  'Resect Insertion Tissue:');
			self.hPickSlice.hText(end+1) = annotation(self.hShowMesh.fig,...
				'textbox', [0.89 0.92 0.04 0.02],...
				'Units', 'normalized', ...
				'EdgeColor', 'none',...
				'FontSize', 8,...
				'VerticalAlignment', 'middle',...
				'String',  'Anterior');
			self.hPickSlice.hEdit(3) = uicontrol(...
				'parent', self.hShowMesh.fig,...
				'KeyPressFcn', @keypress,...
				'style', 'edit',...
				'Units', 'normalized', ...
				'Position', [0.89 0.90 0.04 0.02],...
				'String', num2str(ant), ...
				'callback', @(s,e)noiseAnt(true));
			self.hPickSlice.hText(end+1) = annotation(self.hShowMesh.fig,...
				'textbox', [0.94 0.92 0.04 0.02],...
				'Units', 'normalized', ...
				'EdgeColor', 'none',...
				'FontSize', 8,...
				'VerticalAlignment', 'middle',...
				'String',  'Inferior');
			self.hPickSlice.hEdit(4) = uicontrol(...
				'parent', self.hShowMesh.fig,...
				'KeyPressFcn', @keypress,...
				'style', 'edit',...
				'Units', 'normalized', ...
				'Position', [0.94 0.90 0.04 0.02],...
				'String', num2str(inf), ...
				'callback', @(s,e)noiseInf(true));
			
			self.hPickSlice.hText(end+1) = annotation(self.hShowMesh.fig,...
				'textbox', [0.82 0.87 0.12 0.02],...
				'Units', 'normalized', ...
				'EdgeColor', 'none',...
				'FontSize', 8,...
				'VerticalAlignment', 'middle',...
				'String',  '<Outliers> Show Ecc & Ell Strains Larger Than:');
			self.hPickSlice.hEdit(5) = uicontrol(...
				'parent', self.hShowMesh.fig,...
				'KeyPressFcn', @keypress,...
				'style', 'edit',...
				'Units', 'normalized', ...
				'Position', [0.94 0.87 0.04 0.02],...
				'String', num2str(tol), ...
				'callback', @(s,e)noiseTol());
			self.hPickSlice.hDisplay(1) = annotation(self.hShowMesh.fig,...
				'textbox', [0.84 0.84 0.14 0.02],...
				'Units', 'normalized', ...
				'EdgeColor', 'none',...
				'FontSize', 8,...
				'HorizontalAlignment', 'right',...
				'VerticalAlignment', 'middle');
			self.hPickSlice.hDisplay(2) = annotation(self.hShowMesh.fig,...
				'textbox', [0.84 0.81 0.14 0.02],...
				'Units', 'normalized', ...
				'EdgeColor', 'none',...
				'FontSize', 8,...
				'HorizontalAlignment', 'right',...
				'VerticalAlignment', 'middle');
			
			paramName = fieldnames(self.dataObj.Parameterization);
			POI = strncmpi(paramName,'Longitudinal',5);
			% funcParam = repmat({@returnTrue}, 1, numel(paramName));
			% funcParam{POI} = deal(@CparamWTinsertion);
			idxParam = find(POI);
			% sort "paramName" in a specific order:
			%{ 
			fields = {'Circumferential','Longitudinal','Radial'};
			tmp = cellfun(@(x)strncmpi(paramName,x,4), fields, 'UniformOutput', false);
			paramName = cellfun(@(x)cat(1, paramName{x}), tmp, 'UniformOutput', false);
			tmp = cellfun(@isempty, paramName); paramName(tmp) = [];
			idxParam = 2;
			 %}
			self.hPickSlice.hListbox(1) = uicontrol(...
				'parent', self.hShowMesh.fig,...
				'KeyPressFcn', @keypress,...
				'Style', 'listbox',...
				'String', paramName,...
				'Value', idxParam,...
				'units', 'Normalized',...
				'Position', [0.66 0.01 0.15 0.07],...
				'Callback', @(s,e)pickParam());
						
			nEndo = numel(self.dataObj.Strains);
			if nEndo > 2;set(self.hPickSlice.hButton,'Visible','off'); end% self.deleteSOI(); return;
			idxEndo = nEndo;
			% idxEndo = numel(self.dataObj.Parameterization);
			[idxP2M, cdata] = deal([]);
			self.hPickSlice.hListbox(2) = uicontrol(...
				'parent', self.hShowMesh.fig,...
				'KeyPressFcn', @keypress,...
				'Style', 'listbox',...
				'String', arrayfun(@(x)['Endocardial Mesh #', num2str(x)],1:idxEndo, 'UniformOutput', false),...% cellfun(@(x)['Endocardial Surface #', x],num2cell(1:idxEndo), 'UniformOutput', false)
				'Value', idxEndo,...
				'units', 'Normalized',...
				'Position', [0.82 0.01 0.15 0.07],...
				'Callback', @(s,e)pickEndo());

			% Strains([LVendo; RVendo] X [Septum=1, Freewall=2, Global=3])
			cirLoc = 2;
			self.hPickSlice.hListbox(3) = uicontrol(...
				'parent', self.hShowMesh.fig,...
				'Style', 'listbox',...
				'String', {'Septum', 'Freewall', 'Global'},...
				'Value', cirLoc,...
				'units', 'Normalized',...
				'Position', [0.82 0.1 0.15 0.07],...
				'Callback', @(s,e)pickCirLoc());
			flds={'CC','LL'}; colors={'g','c'}; markers={'x','o'}; str = ' E';%positive  funcCirLoc={@CparamWTinsertion,@CparamWTinsertion,@returnTrue};
			nFlds = numel(flds);
			[strains,idxStrainAnt,idxVertsAnt,idxStrainInf,idxVertsInf,idxTol,locs,param] = deal([]);%,locations
			set(self.hPickSlice.hEdit(3),'enable', 'on');
			for k = 1:nEndo
				self.hPickSlice.hPts(k) = scatter3(self.hShowMesh.ax,[],[],[],10,colors{k},markers{k},'DisplayName',['Noises@',flds{k}]);
			end

			pickParam();
			self.hPickSlice.findNoises = @findNoises;
			plotedit(self.hShowMesh.fig, 'off');
			
			if exist('hwait','var'); hwait.stop; delete(hwait); end
			
			function sliderCB(ind)
				val = get(self.hPickSlice.hSlider(ind),'value');
				set(self.hPickSlice.hEdit(ind),'String',num2str(val));
				slice(ind) = val;
				highlight();
			end
			
			function textCB(ind)
				val = str2double(get(self.hPickSlice.hEdit(ind),'string'));
				if isnan(val) || val<0 || val>1
					set(self.hPickSlice.hEdit(ind),'string',get(self.hPickSlice.hSlider(ind),'value'));
				else
					set(self.hPickSlice.hSlider(ind),'value',val);
					slice(ind) = val;
					highlight();
				end
			end

			function pickParam()
				idxParam = get(self.hPickSlice.hListbox(1),'value');
				if POI(idxParam)
					set(self.hPickSlice.hButton,'enable', 'on');
				else
					set(self.hPickSlice.hButton,'enable', 'off');
				end
				
				pickEndo();
			end

			function pickEndo()
				idxEndo = get(self.hPickSlice.hListbox(2),'value');
				tmp = size(self.dataObj.Parameterization(idxEndo).(paramName{idxParam}),1);
				idxP2M = cellfun(@(x)size(x,1)==tmp,{self.hShowMesh.meshes{1}.vertices});
				idxP2M = find(idxP2M);
				cdata = self.dataObj.Parameterization(idxEndo).(paramName{idxParam});
				% size(self.hShowMesh.meshes{1}(3).vertices,1)
				set(self.hShowMesh.hmeshes,'FaceColor', 'w', 'CDataMapping', 'scaled');%, 'FaceLighting', 'gouraud'
				
				pickCirLoc();
			end

			function pickCirLoc()
				cirLoc = get(self.hPickSlice.hListbox(3),'value');
				
				strains = self.hPickSlice.Strains(idxEndo, cirLoc);
				
				self.hPickSlice.strainLocs = {};
				self.hPickSlice.transStrainLocs = @(s,e,x,y)RefScanner2CurrentImage(strains.Locations, self.dataObj.Interpolants(self.hShowMesh.frame), self.hShowMesh.center, self.hShowMesh.rotmat);
				self.hPickSlice.strainLocs{self.hShowMesh.Frame} = self.hPickSlice.transStrainLocs();
				
				% reverse Circumferential Paramters for the same distance away: being taken care by 
				% set(self.hPickSlice.hEdit(3),'string',num2str(abs(ant-1)));
				% set(self.hPickSlice.hEdit(4),'string',num2str(abs(inf-1)));
				
				noiseAnt(); noiseInf(); noiseTol();
			end

			function noiseAnt(alter)
				val = str2double(get(self.hPickSlice.hEdit(3),'string'));
				if isnan(val) || val<0 || val>1
					set(self.hPickSlice.hEdit(3),'string',num2str(tol));
				else
					ant = val;
					idxStrainAnt = CparamWTinsertion(strains.Parameterization.Circumferential, ant);
					idxVertsAnt = CparamWTinsertion(self.dataObj.Parameterization(idxEndo).Circumferential, ant);
					% idxStrainAnt = funcCirLoc{cirLoc}(strains.Parameterization.Circumferential, ant);
					% idxVertsAnt = funcParam{idxParam}(self.dataObj.Parameterization(idxEndo).Circumferential, ant);
					
					if exist('alter','var'); findNoises(); end
				end
			end
			
			function noiseInf(alter)
				val = str2double(get(self.hPickSlice.hEdit(4),'string'));
				if isnan(val) || val<0 || val>1
					set(self.hPickSlice.hEdit(4),'string',num2str(tol));
				else
					inf = val;
					idxStrainInf = CparamWTinsertion(strains.Parameterization.reverseCircumferential, inf);
					idxVertsInf = CparamWTinsertion(self.dataObj.Parameterization(idxEndo).reverseCircumferential, inf);
					% idxStrainInf = funcCirLoc{cirLoc}(strains.Parameterization.reverseCircumferential, inf);
					% idxVertsInf = funcParam{idxParam}(self.dataObj.Parameterization(idxEndo).reverseCircumferential, inf);
					
					if exist('alter','var'); findNoises(); end
				end
			end
			
			function noiseTol()
				val = str2double(get(self.hPickSlice.hEdit(5),'string'));
				if isnan(val) || val<-1 || val>1
					set(self.hPickSlice.hEdit(5),'string',num2str(tol));
				else
					tol = val;
					idxTol = cellfun(@(x)strains.(x)>tol, flds,'UniformOutput', false);
					findNoises();
				end
			end
			
			function findNoises
				for k = 1:nFlds
					idxLocs = idxTol{k}(:,self.hShowMesh.Frame)&idxStrainAnt&idxStrainInf;
					locs{k} = self.hPickSlice.strainLocs{self.hShowMesh.Frame}(idxLocs,:);
					param{k} = strains.Parameterization.(paramName{idxParam})(idxLocs);
				end
%{
strains.Parameterization.(paramName{idxParam})(idxLocs);
abs(ans-slice(1)) < slice(2);
k=find(ans);
tmp=find(idxLocs);
strains.Parameterization.Circumferential(tmp(k(1)))
find(strains.Parameterization.Circumferential==1);
%}
				
				highlight();
			end

			function highlight()
				ind = abs(cdata-slice(1)) < slice(2)&idxVertsAnt&idxVertsInf;
				set(self.hShowMesh.hmeshes(idxP2M), 'FaceColor', 'flat', 'FaceVertexCData', int8(ind), 'CDataMapping', 'scaled');
				
				if POI(idxParam)
					points = get(self.hShowMesh.hmeshes(idxP2M),'Vertices');
					points = points(ind,3);
					pts = get(self.hShowMesh.hmeshes(1),'Vertices');
					pts = pts(:,3);
					try ind = (pts>min(points)) & (pts<max(points)); end
					set(self.hShowMesh.hmeshes(1), 'FaceColor', 'flat', 'FaceVertexCData', int8(ind));
				end
				
				for k = 1:nEndo
					ind = abs(param{k}-slice(1)) < slice(2);
					points = locs{k}(ind,:);
					set(self.hPickSlice.hPts(k),'XData',points(:,1),'YData',points(:,2),'ZData',points(:,3));
					tmp = size(points,1);
					set(self.hPickSlice.hDisplay(k),'String',[num2str(tmp),str,flds{k},' (marker: ',markers{k},')']);
				end
				
				drawnow;
			end
						
			function regionalStrains()
				%{ 
				limit = [slice(1)-slice(2),slice(1)+slice(2)];
				if limit(1)<0 || limit(2)>1
					errordlg('The lower limit of the slice should be large than 0 and the upper limit should be less than 1 !','Incorrect Setting of the Slice of Interest','modal');
					return
				end
				 %}
				hwait = waitbartimer();
				cleanupObj = onCleanup(@(x)delete(hwait(isvalid(hwait))));
				hwait.String = 'Setting Up the Slice of Interest...';
				hwait.WindowStyle = 'modal';
				hwait.AllowClose = false;
				hwait.start;	
								
				%% find mapping: queryParams -> StrainElement.Parameterization.Longitudinal
				%{ 
				queryParams = [0.167861786754660 0.213685102311070 0.262646080618122 0.392233662306029 0.752770068547287 0.981737868949697 0.421596620482415 0.784201050416771 0.565422944685641 0.256833049127323 0.273272803232508]';
				for k = 1:nEndo
					% test if tol is small enough for a unique finding: 1e-5 is barely enough
					% indices = arrayfun(@(x)abs(self.dataObj.Parameterization(k).Longitudinal-x) <= 1e-5,queryParams,'UniformOutput', false);
					indices = arrayfun(@(x)abs(self.dataObj.Strains(k).Parameterization.Longitudinal-x) <= 1e-6,queryParams,'UniformOutput', false);
					inds{k} = cellfun(@find,indices,'UniformOutput', false);
				end
				 %}
				
				%% considering insertion pt shift when longitudinal location varys:
				Strains = self.hPickSlice.Strains(idxEndo, 3);
				fields = fieldnames(Strains);
				ind = CparamWTinsertion(Strains.Parameterization.Circumferential, ant) & CparamWTinsertion(Strains.Parameterization.reverseCircumferential, inf);
				for k = 1:numel(fields)
					value = Strains.(fields{k});
					if isstruct(value)
						tmp = fieldnames(value);
						for ii = 1:numel(tmp)
							Strains.(fields{k}).(tmp{ii}) = value.(tmp{ii})(ind,:);
						end
					else
						Strains.(fields{k}) = value(ind,:);
					end
				end
				
				lsegments = [0-eps,slice(1)-slice(2),slice(1)+slice(2)];
				% Break it into 3 segments longitudinally using bsxfun
				L = bsxfun(@gt, Strains.Parameterization.Longitudinal, lsegments);
				Lseg = max(cumsum(L, 2), [], 2);

				% A rectangle with a largest width approximate the actual wedge
				%{
				tmp = Strains.Parameterization.(paramName{idxParam})(ind);
				tmp = abs(tmp-slice(1)) < slice(2);
				csegments = Strains.Parameterization.Circumferential(ind);
				csegments = csegments(tmp);
				csegments = [min(csegments), max(csegments)];
				 %}
				clear tmp;
				tmp(1,:) = interp1(self.hPickSlice.insertion{idxEndo}(:,1), self.hPickSlice.insertion{idxEndo}(:,2), [slice(1)-slice(2),slice(1)+slice(2)]);
				tmp(2,:) = interp1(self.hPickSlice.insertion{idxEndo}(:,1), self.hPickSlice.insertion{idxEndo}(:,3), [slice(1)-slice(2),slice(1)+slice(2)]);
				clear csegments;
				% A rectangle with a smallest width approximate the actual wedge
				if idxEndo < 2
					csegments(1,:) = [linspace(ant, min(tmp(1,:))-inf, 3), linspace(max(tmp(1,:))+inf, 1-ant, 5)];
					csegments(2,:) = [linspace(inf, min(tmp(2,:))-ant, 5), linspace(max(tmp(2,:))+ant, 1-inf, 3)];
				else
					csegments(1,:) = [linspace(ant, min(tmp(1,:))-inf, 5), linspace(max(tmp(1,:))+inf, 1-ant, 3)];
					csegments(2,:) = [linspace(inf, min(tmp(2,:))-ant, 3), linspace(max(tmp(2,:))+ant, 1-inf, 5)];
					% csegments = [0-eps, csegments{1}, csegments{2}];
				end
				% compare ONE Parameterization value of each element to: 24 Parameterization values of all circumferential segments @ same longitudinal location -> assign a logical value for each element
				C = bsxfun(@gt, Strains.Parameterization.Circumferential, csegments(1,:));%Greater than
				Cseg = max(cumsum(C, 2), [], 2);
				indices{1} = accumarray([Cseg(:), Lseg(:)], 1:size(L,1), [], @(x){x(:).'});
				C = bsxfun(@gt, Strains.Parameterization.reverseCircumferential, csegments(2,:));
				Cseg = max(cumsum(C, 2), [], 2);
				indices{2} = accumarray([Cseg(:), Lseg(:)], 1:size(L,1), [], @(x){x(:).'});
				
%{ 
points = Strains.Locations(cat(2,segments{1}),:);
points = Strains.Locations(indices{7,  2},:);
tmp = Strains.Parameterization.Circumferential(indices{7,  2},:);
points = bsxfun(@plus, bsxfun(@minus, points, self.hShowMesh.center) * self.hShowMesh.rotmat, self.hShowMesh.center);
scatter3(self.hShowMesh.ax,points(:,1),points(:,2),points(:,3),10,'g','x');
 %}
			
				% nSegments is: 17 for LV, 12 for RV (default)
				nSegs = [17,12];
				nFrames = size(Strains.(fields{1}),2);
				
				% Add the folder plugins>dense3D_plugin into search path:
				import plugins.dense3D_plugin.*
				% Use the mean LV contraction as a reference for DelayTimes
				tmp = self.dataObj.Strains(1).p2;
				reference = mean(tmp, 1);
				RD = RegionalDyssynchrony(permute(Strains.p2, [1 3 2]));
				delays = RD.computeRegionalDelays(reference);
	
				if idxEndo < 2
					% Group Strains based upon their parameterization
					% if cirLoc < 2
					% segments = {
						% [unique([indices{1}{1, 2},indices{2}{6, 2}])]
						% [unique([indices{1}{2, 2},indices{2}{7, 2}])]
					% };
					% ind = [8,9];					
					% elseif cirLoc < 3
					% segments = {
						% [unique([indices{1}{4, 2},indices{2}{1, 2}])]
						% [unique([indices{1}{5, 2},indices{2}{2, 2}])]
						% [unique([indices{1}{6, 2},indices{2}{3, 2}])]
						% [unique([indices{1}{7, 2},indices{2}{4, 2}])]
					% };
					% ind = [10:12,7];					
					% else
					segments = {
						[unique([indices{1}{1, 2},indices{2}{6, 2}])]
						[unique([indices{1}{2, 2},indices{2}{7, 2}])]
						[unique([indices{1}{4, 2},indices{2}{1, 2}])]
						[unique([indices{1}{5, 2},indices{2}{2, 2}])]
						[unique([indices{1}{6, 2},indices{2}{3, 2}])]
						[unique([indices{1}{7, 2},indices{2}{4, 2}])]
					};
					ind = [8:12,7];
					% end
					
					for k = 1:numel(fields)
						value = Strains.(fields{k});
						if isstruct(value); continue; end
						output(1).(fields{k}) = NaN(nSegs(1),size(value,2));
						func = @(x)mean(Strains.(fields{k})(x,:), 1);
						results.(fields{k}) = cellfun(func, segments, 'UniformOutput', 0);
						output(1).(fields{k})(ind,:) = cat(1, results.(fields{k}){:});
						if nEndo > 1
							output(2).(fields{k}) = NaN(nSegs(2),size(value,2));
						end
					end
					[output.CURE,output.RURE,output.LURE]=deal(NaN(3,nFrames));
					output(1).CURE(2,:) = CURE(output(1).CC(ind,:).');
					output(1).RURE(2,:) = CURE(output(1).RR(ind,:).');
					output(1).LURE(2,:) = CURE(output(1).LL(ind,:).');
					output(1).CLShearAngle = rad2deg(self.dataObj.torsion(output(1)));
					output(1).DelayTimes = NaN(nSegs(1),1);
					output(1).DelayTimes(ind,:) = cellfun(@(x)mean(delays(x)), segments);
					if nEndo > 1
						output(2).CLShearAngle = NaN(nSegs(2),nFrames);
						output(2).DelayTimes = NaN(nSegs(2),1);
					end
				else
					segments = {
						[unique([indices{1}{1, 2},indices{2}{4, 2}])]
						[unique([indices{1}{2, 2},indices{2}{5, 2}])]
						[unique([indices{1}{3, 2},indices{2}{6, 2}])]
						[unique([indices{1}{4, 2},indices{2}{7, 2}])]
						[unique([indices{1}{6, 2},indices{2}{1, 2}])]
						[unique([indices{1}{7, 2},indices{2}{2, 2}])]
					};
					
					for k = 1:numel(fields)
						value = Strains.(fields{k});
						if isstruct(value); continue; end
						output(1).(fields{k}) = NaN(nSegs(1),size(value,2));
						output(2).(fields{k}) = NaN(nSegs(2),size(value,2));
						func = @(x)mean(Strains.(fields{k})(x,:), 1);
						results.(fields{k}) = cellfun(func, segments, 'UniformOutput', 0);
						output(1).(fields{k})([9,8],:) = cat(1, results.(fields{k}){[5,6]});
						output(2).(fields{k})(5:8,:) = cat(1, results.(fields{k}){1:4});
					end
					[output.CURE,output.RURE,output.LURE]=deal(NaN(3,nFrames));
					%% RV Freewall w.t. Septum:
					%{ 
					output(2).CURE(2,:) = CURE([output(2).CC(5:8,:);output(1).CC([9,8],:)].');
					% Septal nSegments = 2 < 3
					output(1).CURE(2,:) = CURE(output(1).CC([9,8],:).');
					output(1).RURE(2,:) = CURE(output(1).RR([9,8],:).');
					output(1).LURE(2,:) = CURE(output(1).LL([9,8],:).');
					 %}
					output(2).CURE(2,:) = CURE(output(2).CC(5:8,:).');
					output(2).RURE(2,:) = CURE(output(2).RR(5:8,:).');
					output(2).LURE(2,:) = CURE(output(2).LL(5:8,:).');
					output(1).CLShearAngle = rad2deg(self.dataObj.torsion(output(1)));
					output(2).CLShearAngle = rad2deg(self.dataObj.torsion(output(2)));
					output(1).DelayTimes = NaN(nSegs(1),1);
					output(1).DelayTimes([9,8],:) = cellfun(@(x)mean(delays(x)), segments([5,6]));
					output(2).DelayTimes = NaN(nSegs(2),1);
					output(2).DelayTimes(5:8,:) = cellfun(@(x)mean(delays(x)), segments(1:4));
				end

				% Create output struct to save everything in
				output(idxEndo).Segmentation = segments;

				self.viewerObj.cache.regionalStrain = output;
				
				button = questdlg('Do you wanna save the result?','Yes');
				if strcmpi(button,'Yes')				
					fields = {'RR','CC','LL','CL'};
					% fields = {'CL'};
					workbookNm =fullfile(parentDir,'RV_polar_strains_Baseline.xlsx');
					for k = 1:numel(fields)
						sheetNm = ['E',lower(fields{k})];
						tmp = cat(1, results.(fields{k}){:});
						tmp = num2cell(tmp(:,self.hShowMesh.Frame));
						try
							[~,~, existingData] = xlsread(workbookNm,sheetNm);
						catch
							existingData = [{'Dataset';'self.hShowMesh.Frame No';'Slice'};num2cell([1:6]')];
						end
						try
							xlswrite(workbookNm,[existingData,[{fname};{self.hShowMesh.frame};sprintf('+%.3g',slice);tmp]],sheetNm);
						catch
							keyboard;
						end
					end
					savePath = get(self.configObj, 'Location.LoadWorkspace', pwd);
					save(fullfile(savePath,['RV_polar_strains_',sprintf('s%.3g',slice),'_f',num2str(self.hShowMesh.frame),'.mat']),'results');
					hgexport(self.hShowMesh.fig,fullfile(savePath,['RV_polar_strains_',sprintf('s%.3g',slice),'_f',num2str(self.hShowMesh.frame),'.jpg']),hgexport('factorystyle'), 'Format', 'jpeg');
					msg = {'Results saved in the file:',workbookNm};
				else
					msg = {};
				end
				
				msgbox({'Definition Completed for Slice of Interest:',['SeedPoint Variable: Slice=',num2str(slice)],'Please seclect any option for Bulleseye Plot to see results!',msg{:}},'Slice of Interest for Strain Analysis');%,'modal'
			end
		end

        function parameterize(self)
            import plugins.DENSE3D_Plugin_4CrescentOrgan.*
			
			rois = [self.dataObj.Data.ROIInfo];
			SAslices = strwcmpi({rois.ROIType}, '*sa*');
			idxBase = find(SAslices, 1, 'first');
			idxApex = find(SAslices, 1, 'last');
			
			if ~self.hShowMesh.isBiv
				sqInfo  = cat(3, self.dataObj.Data.SequenceInfo);
				sqInfo  = squeeze(sqInfo(1,1,:));
				self.hPickSlice.imagePos = catField(sqInfo,'ImagePositionPatient');
				% self.hPickSlice.imagePos = cat(2,sqInfo.ImagePositionPatient)';
				% self.hPickSlice.sliceLocs = catField(sqInfo,'SliceLocation');
				tmp = [sqInfo.tform]; tmp = [tmp.tdata];
				self.hPickSlice.tdata.T = cat(3,tmp.T);
				self.hPickSlice.tdata.Tinv = cat(3,tmp.Tinv);

				anInfo = [self.dataObj.Data.AnalysisInfo];
%{
files = DENSE2DobjectFileName()
anInfo = cellfun(@(x)load(x, 'AnalysisInfo', '-mat'), files);
anInfo = [anInfo.AnalysisInfo];
%}
				try self.hPickSlice.PositionC = catField(anInfo,'SegDistribution'); end
				try self.hPickSlice.PositionA = catField(anInfo,'PositionA'); end
				
				str = 'Nseg'; multiplier = 1;
				if ~isfield(anInfo, str); str = 'Nmodel'; multiplier = 10; end
				tmp = catField(anInfo,str);
				self.hPickSlice.Nseg = tmp(1) * multiplier;
				
				str = 'PositionB';
				if ~isfield(anInfo, str)
					anInfo = [self.dataObj.Data.DENSEInfo];
					str = 'SegPos';
				end
				self.hPickSlice.PositionB = catField(anInfo,str);
				% self.hPickSlice.PositionB = cat(1,anInfo.(str));
				

			end
			
			basalSlice = self.dataObj.Data(idxBase).SequenceInfo(1);
			apicalSlice = self.dataObj.Data(idxApex).SequenceInfo(1);
			self.hPickSlice.normal = cross(apicalSlice.ImageOrientationPatient(1:3), apicalSlice.ImageOrientationPatient(4:6));
			try
				self.hPickSlice.length = point2planeDistance(basalSlice.ImagePositionPatient.', apicalSlice.ImagePositionPatient, self.hPickSlice.normal);
				% make sure that the normal points towards the base
				if self.hPickSlice.length < 0; self.hPickSlice.normal = -self.hPickSlice.normal; end
			catch
				locs = cellfun(@(x)x(1).SliceLocation, {self.dataObj.Data.SequenceInfo});
				self.hPickSlice.length = max(locs)-min(locs);
			end
			% pt2planeDist = bsxfun(@minus, points, pt(:)');
			% pt2planeDist = pt2planeDist * (self.hPickSlice.normal(:) ./ norm(self.hPickSlice.normal));
			self.hPickSlice.gdApex = 1/(numel(rois)-1);
			
            % Compute the parameterization for each endocardial mesh
			nEndo = numel(self.dataObj.EndocardialMesh);
            arrayfun(@(x)self.longitudinalParameterize(x, apicalSlice.ImagePositionPatient), 1:nEndo);
            arrayfun(@(x)circumferentialParameterize(self, x), 1:nEndo);
        end
        
		
		function longitudinalParameterize(self, index, apicalSlice)
            import plugins.DENSE3D_Plugin_4CrescentOrgan.*
			verts = self.dataObj.EndocardialMesh(index).vertices;
			faces = self.dataObj.EndocardialMesh(index).faces;
            apex = self.dataObj.Apex(index,:);
			
            %% Longitudinal parameterization
			apex2apicalSlice = point2planeDistance(apex, apicalSlice, self.hPickSlice.normal);
			
			%% Acceptable auto-generated apex: extend more than 10% of the cavity length below the Apical Slice
			% button = questdlg({'Previously Defined Slice of Interest Detected!','Choose what you want to do with it'},'Warning','Delete','Do nothing','Delete');
			if abs(apex2apicalSlice/self.hPickSlice.length) > self.hPickSlice.gdApex
				% Compute the distance from all verts to the apex
				dist2apex = sqrt(sum(bsxfun(@minus, apex, verts).^2,2));
			else
				% Find seed points closed to the apical slice
				apicalPoints = intersectVerts(verts, faces, apicalSlice, self.hPickSlice.normal);
% Debug:
%{
[~,~,apicalPoints] = half_space_intersect(verts, faces, apicalSlice, self.hPickSlice.normal, 'Cap', false);
hold on; scatter3(verts(apicalPoints,1),verts(apicalPoints,2),verts(apicalPoints,3),10,'b','*');
 %}
			
				% Another quick way: use image coordinates after rotateMesh
				%{ 
				apicalPoints = abs(self.hShowMesh.meshes{1}(index+1).vertices(:,end)-self.hShowMesh.ROIs{1}{idxApex,end}(1,end)) < 1e-1;
				%}
				% Now find the minimum distance to a apical point
				lparams = sqrt(sum(bsxfun(@minus, apicalPoints.', permute(verts, [2 3 1])).^2, 1));
				dist2apex = min(lparams, [], 2);
								
				dists = point2planeDistance(verts, apicalSlice, self.hPickSlice.normal);
				apexPoints = dists < 0;
				% numel(find(apexPoints))%normal direction
				dist2apex(apexPoints) = deal(0);
			end

			% Find the open edge of the mesh to use as the base seed points
			basalPoints = unique(outline(faces));
% Debug:
%{
self.dataObj.preview();
hold on; scatter3(verts(basalPoints,1),verts(basalPoints,2),verts(basalPoints,3),10,'b','*');

ind = zeros(size(faces,1),1); ind(basalPoints) = 1;
hmeshes = findall(allchild(gca),'-property','FaceColor');
set(hmeshes,'FaceColor', 'w');
set(hmeshes(2), 'FaceColor', 'flat', 'FaceVertexCData', ind, 'CDataMapping', 'scaled');	
 %}
			
			lparams = sqrt(sum(bsxfun(@minus, verts(basalPoints,:).', permute(verts, [2 3 1])).^2, 1));
			dist2base = min(lparams, [], 2);

			% Normalize the distance between the base and apex
			self.dataObj.Parameterization(index).Longitudinal = dist2apex(:) ./ (dist2apex(:) + dist2base(:));
        end
		
        function insertion = anteriorInsertion(self, location, longParam)% index, 
			if self.hShowMesh.isBiv
				import plugins.DENSE3D_Plugin_4CrescentOrgan.*
				epiOutline = intersectVerts(self.dataObj.EpicardialMesh.vertices, self.dataObj.EpicardialMesh.faces, location, self.hPickSlice.normal, true);
				nEndo = numel(self.dataObj.EndocardialMesh);
				
				%% produce the INTERSECTED LINE for finding  insertion points
				% Alternative method: bug @ endoOutlines{1} has to be LVendo while endoOutlines{1} for RVendo, otherwise the anteriorInsertion will be flipped
				%{ 
				if size(location, 1) > 1
				% "isoline" is delivered here:
					endoOutlines{1} = location;
					location = mean(endoOutlines{1});
					if nEndo > 1
						ind = [1:nEndo]; ind(index) = [];
						tmp = arrayfun(@(x)intersectVerts(self.dataObj.EndocardialMesh(x).vertices, self.dataObj.EndocardialMesh(x).faces, location, self.hPickSlice.normal), ind, 'UniformOutput', false);
						if any(cellfun(@isempty, tmp)); insertion = []; return; end
						endoOutlines = cat(2, endoOutlines, tmp{:});
					end
					alternative = false;
				else
					alternative = true;
				end
				 %}
%{ 
plot3(isoline(:,1),isoline(:,2),isoline(:,3),'LineWidth',1.5,'color', 'b')
scatter3(endoOutlines{2}(:,1),endoOutlines{2}(:,2),endoOutlines{2}(:,3),10,'b','o');
scatter3(epiOutline(:,1),epiOutline(:,2),epiOutline(:,3),10,'r','x');
%}

				%% Copyright (c) of the following section of the codes: 2016 DENSEanalysis Contributors
				% The following section of the codes are from the function "playbackFcn" in "DENSE3D.m".
				% Find the EPI point that is equidistance from each contour
				[epidist, epidistind] = deal(nan(size(epiOutline, 1), nEndo));
				for k = 1:nEndo
					% if alternative 
						endoOutlines{k} = intersectVerts(self.dataObj.EndocardialMesh(k).vertices, self.dataObj.EndocardialMesh(k).faces, location, self.hPickSlice.normal);
						if isempty(endoOutlines{k}); insertion = []; return; end
					% end
					% Compute the point-wise distance between each Endo and the EPI
					distances = bsxfun(@minus, epiOutline.', permute(endoOutlines{k}, [2 3 1]));
					distances = squeeze(sum(distances.^2, 1));

					[epidist(:,k), epidistind(:,k)] = min(distances, [], 2);
				end

				% Repeat the last point
				epidist     = epidist([1:end 1],:);
				epidistind  = epidistind([1:end 1],:);

				% Find where the endos cross which will give us the
				% indices of the anterior and inferior insertion points
				sgn = sign(diff(epidist, [], 2));
				insertions = logical(diff(sgn));

				% Make sure that we only picked up two of them
				tmp = sum(insertions);
				if  tmp ~= 2
					errordlg({['Found ',num2str(tmp),' insertion points for the contour line with the longitudinal parameter: ',num2str(longParam)],'Try to use alternative method in the codes above!'},'Debug Mode','modal');
					keyboard;
				end

				% Now compute the XYZ coordinates of each of these two
				% points by averaging the closest points on the endos and
				% epicardial boundary
				ind = num2cell(epidistind(insertions,:));
				epipoints = epiOutline(insertions,:);
				tmp = cellfun(@(x,y)x(y,:), endoOutlines, ind(1,:), 'uni', 0);
				one = mean(cat(1, epipoints(1,:), tmp{:}));
				tmp = cellfun(@(x,y)x(y,:), endoOutlines, ind(2,:), 'uni', 0);
				two = mean(cat(1, epipoints(2,:), tmp{:}));

				% Figure out which one is actually the anterior one by
				% using the cross product with a vector to the apex
				centers = cellfun(@(x)mean(x, 1), endoOutlines, 'uni', 0);
				N = cross(centers{2} - centers{1}, ...
						  self.dataObj.Apex(1,:) - centers{1});
				% Positive distance is going to be anterior
				side = sign(point2planeDistance([one; two], centers{1}, N));
				
				%% Copyright (c) of the following section of the codes: Zhanqiu Liu (lafeir.lew@gmail.com)
				% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)
				% Last Modified: 19:32 June 27, 2018				
				tmp = {'INFERIOR', 'ANTERIOR'};
				if  prod(side) > 0
					errordlg({['Two insertion points are at the same ',tmp{side(1)<0+1},' side for the contour line with the longitudinal parameter: ',num2str(longParam)],'Try to use alternative method in the codes above!'},'Debug Mode','modal');
					keyboard;
				end
				
				if side(1)<0; insertion=[two; one]; else; insertion=[one; two]; end;
			
            else
            %% LV
				distances = bsxfun(@minus,self.hPickSlice.imagePos,location);
				nROIs = size(distances,1);
				distances = dot(distances, repmat(reshape(self.hPickSlice.normal,1,[]), nROIs, 1), 2);
				sgn = sign(distances); idxLocs = find(diff(sgn));
				tmp = numel(idxLocs);
				switch tmp
					case 1
						idxLocs = [idxLocs, idxLocs+1];
					case 0
						idxLocs = [nROIs, nROIs];
					otherwise
						if tmp~=2 || ~sum(distances(idxLocs)==0)
							errordlg({'There is duplicate input present in the 3D Workspaces!','Or the inputs are not ordered descendingly!'},'Debug Mode','modal');
							keyboard;
						end
				end
				
				%% weights: weigh more when closer
				% tmp = abs(distances(idxLocs(end:-1:1)));% slow
				tmp = abs(distances(idxLocs)); tmp = flip(tmp);
				weight = reshape(tmp/(tmp(1)+tmp(end)),[],1);

				%% ANTERIOR insertion point:
				one = dot(self.hPickSlice.PositionB(idxLocs,:), repmat(weight,1,2), 1);
				
				%% INFERIOR insertion point:
				try
					tmp = cat(1,self.hPickSlice.PositionC{:,2});
					idx = isnan(tmp); ind = idx(idxLocs);
					switch sum(ind)
						case 0
							tmp = tmp(idxLocs);
						case 1
							tmp = tmp(idxLocs); tmp(ind) = tmp(~ind);
						case 2
							tmp = repmat(mean(tmp(~idx)),2,1);
						otherwise
							errordlg({'There is duplicate input present in the 3D Workspaces!','Or the inputs are not ordered descendingly!'},'Debug Mode','modal');
							keyboard;
					end
					
					two = dot(tmp, weight, 1);
					theta = (two/self.hPickSlice.Nseg) * 2*pi;
				catch
					% Rotate the anterior insertion point by pi/3 to get the inferior insertion point:
					theta = 20*pi / self.hPickSlice.Nseg;
				end
                % Rotate the other way if we need it to go clockwise
                if self.dataObj.Data(1).AnalysisInfo.Clockwise
                    theta = -theta;
                end
				
                % Create a rotation matrix to shift the anterior insertion
                % to the correct location
                R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
                % Approximated center of intersected contour line:
				center = dot(self.hPickSlice.PositionA(idxLocs,:), repmat(weight,1,2), 1);
                two = ((one - center) * R) + center;
				
				%% forward spatial transformation
				tform = self.dataObj.Data(1).SequenceInfo(1).tform;
				tmp = arrayfun(@(x)repmat(x, 4, 4),weight,'UniformOutput',false);
				tmp = cat(3,tmp{:});
				tform.tdata.T = dot(self.hPickSlice.tdata.T(:,:,idxLocs), tmp, 3);
				tform.tdata.Tinv = dot(self.hPickSlice.tdata.Tinv(:,:,idxLocs), tmp, 3);
				
                one = tformfwd(tform, [one, 0]);
% scatter3(ans(:,1),ans(:,2),ans(:,3),30,'k','o');
				% one = ([one, 0] - self.hShowMesh.center) * self.hShowMesh.rotmat + self.hShowMesh.center;
                two = tformfwd(tform, [two, 0]);
				insertion=[one; two];				
% scatter3(insertion(:,1),insertion(:,2),insertion(:,3),30,'b','d');
            end
        end		
		
		function saveGIF(self)
		% Last Modified: 23:03 June 29, 2018
		% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)
			
			if isempty(self.hShowMesh); return; end
			
			import plugins.DENSE3D_Plugin_4CrescentOrgan.*

			hwait = waitbartimer();
			cleanupObj = onCleanup(@(x)delete(hwait(isvalid(hwait))));
			hwait.WindowStyle = 'modal';
			hwait.AllowClose = false;
            hwait.Visible = 'on';
            hwait.String = 'Saving the movie...';
			hwait.start;
			drawnow
			
			[uifile, uipath, uipopup] = uiputfile({'*.gif','GIF movie(*.gif)'},'Select a path to save',pwd);
			if ~uipopup; return; end

			button = questdlg('Choose the type of face rendering:','Save GIF','Smooth','Silly Putty','Wire Frame','Smooth');%{,'(close this window to reset the viwer)'}
			
			% LnVis = get(self.hShowMesh.hROIs, 'Visible');
			set(self.hShowMesh.hROIs, 'Visible', 'off');
			% nSurfaces = numel(self.hShowMesh.meshes{1});
			set(self.hShowMesh.fig, 'color', 'k');
            switch lower(button)
                case 'smooth'
					set(self.hShowMesh.hmeshes, 'FaceLighting', 'gouraud', 'EdgeLighting', 'gouraud', 'BackFaceLighting', 'reverselit', 'AmbientStrength', .6);
					set(self.hShowMesh.hmeshes(1), 'FaceColor', 'w', 'FaceAlpha', 0.5, 'EdgeColor', 'g', 'EdgeAlpha', 0.25);
					set(self.hShowMesh.hmeshes(2), 'FaceColor', 'b', 'FaceAlpha', 0.75, 'EdgeColor', 'k', 'EdgeAlpha', 0.5);
					try set(self.hShowMesh.hmeshes(3), 'FaceColor', 'r', 'FaceAlpha', 0.75, 'EdgeColor', 'k', 'EdgeAlpha', 0.5); end
					light
                case 'silly putty'
					set(self.hShowMesh.hmeshes, 'FaceLighting', 'gouraud', 'EdgeLighting', 'gouraud', 'BackFaceLighting', 'reverselit', 'AmbientStrength', .6);
					set(self.hShowMesh.hmeshes, 'FaceColor', [1 2/3 2/3], 'EdgeColor', 'none');
					light
                otherwise
                % case 'wire frame'
					set(self.hShowMesh.hmeshes, 'FaceLighting', 'flat', 'EdgeLighting', 'flat', 'BackFaceLighting', 'none');
					set(self.hShowMesh.hmeshes, 'FaceColor', 'none', 'EdgeColor', 'w', 'EdgeAlpha', 0.1);
            end

			gif = animatedgif(fullfile(uipath,uifile), 'LoopCount', Inf, 'DelayTime', 0.1);
			for k = 0:self.hShowMesh.hplaybar.Max
				self.hShowMesh.hplaybar.Value = k;
                gif.addFrame(getframe(self.hShowMesh.ax));
            end
			
			%% reset the viewer
			% tmp = arrayfun(@mat2cell,self.hShowMesh.hROIs);
			% arrayfun(@(s,e)set(s, 'Visible', e), tmp{:}, LnVis{:})
			delete(findobj(self.hShowMesh.ax, 'type', 'light'));
			set(self.hShowMesh.fig, 'color', self.BackgroundColor);
			set(self.hShowMesh.hmeshes, 'FaceLighting', 'none', 'EdgeLighting', 'none', 'BackFaceLighting', 'unlit', 'AmbientStrength', .3);
			set(self.hShowMesh.hmeshes, 'FaceColor', 'w', 'FaceAlpha', 1, 'EdgeColor', 'k', 'EdgeAlpha', 1);
		end
		
    end
	
    methods (Access = 'protected')

        function resize(self)
            resize@DataViewer(self);
            if ~isempty(self.api.resizeFcn)
                try
                    self.api.resizeFcn();
                catch
                    self.api.resizeFcn = [];
                end
            end
        end

        function playback(self)
            if ~isempty(self.api.playbackFcn)
                self.api.playbackFcn();
            end
        end

        function reset(self)
 		% Copyright (c) 2016, Jonathan Suever
           self.hlisten_playbar.Enabled = false;
            stop(self.hplaybar);
            self.hplaybar.Min = 1;
            self.hplaybar.Max = 0;
            self.hplaybar.Visible = 'off';

            if ~isempty(self.api.deleteFcn)
                self.api.deleteFcn();
            end

            self.api = self.template;

            self.exportaxes = false;
            self.exportrect = [];
            self.isAllowExportImage = false;
            self.isAllowExportVideo = false;
        end
    end
	
end


function deleteHandles(struc)
% Last Modified: 17:28 August 7, 2018
% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)
	fields = fieldnames(struc);
	ind = cellfun(@(x)ishandle(struc.(x)),fields,'UniformOutput',false);
	ind = logical(cellfun(@(x)prod(x(:)),ind));
	cellfun(@(x)delete(struc.(x)),fields(ind));
end

function array = catField(structure, fieldName)
% Last Modified: 17:28 August 7, 2018
% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)
	if size(structure(1).(fieldName),1) == 1
		array = cat(1,structure.(fieldName));				
	else
		array = cat(2,structure.(fieldName))';
	end
end

function ind = CparamWTinsertion(array, tol)
	% ALWAYS NOT included: Cparam = 0 & 1->kill tons of strains @ insertion
	ind = array>tol & (1-array)>tol;
end

function locations = RefScanner2CurrentImage(locations, interpolants, center, rotmat)
	locations = locations + interpolants.query(locations);
	locations = bsxfun(@plus, bsxfun(@minus, locations, center) * rotmat, center);
end

% function out = returnTrue(varargin)
	% out = true;
% end

% Empty KeyPressFcn: no keystrokes go through to the parent figure.
function keypress(~,~)
	return
end