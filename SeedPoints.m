function [center,initial_nor_R,final_nor_R,initial_tan_R,final_tan_R,frame,slice,ant,inf,tol,fname,parentDir] = SeedPoints(settingfile)
	%
	% INPUTS:
	%   self:   Object,  Instance of the TestPlugin object
	%   handles:    struct, guidata structure from the DENSEanalysis
	%               GUI. This provides access to all data and controls 
	%               in the GUI.
	% Last Modified: 12:34 July 12, 2018
	% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)

	center = [];
	initial_nor_R = []; % [mm]
	final_nor_R = []; % [mm]
	initial_tan_R = []; % [mm]
	final_tan_R = []; % [mm]
	frame = [];
	slice = [];
	ant = [];
	inf = [];
	tol = [];

	%% Read pre-set SeedPt & fitting parameters:
	fname = get(settingfile, 'Location.LoadWorkspace', pwd);
	[parentDir,fname] = fileparts(fname);
	if isempty(fname); [parentDir,fname]=fileparts(parentDir); end
	fname = regexp(fname,'\d*','Match');
	if iscell(fname); fname = fname{:}; end;
	varargout
	if exist('fname','var')	
		switch lower(fname)
		case '20150628'%% dataset 20150628
			center = [2.51346206948696	7.08576724017031	-3.12295877567370];
			%% dataset 20150628discontinued
			initial_nor_R = .3 % [mm]
			final_nor_R = 1.1 % [mm]
			initial_tan_R = .1 % [mm]
			final_tan_R = .5 % [mm]
			frame = 10;
			%% dataset 20150628
			% initial_nor_R = .4 % [mm]
			% final_nor_R = .8 % [mm]
			% initial_tan_R = .3 % [mm]
			% final_tan_R = .5 % [mm]
		case '20150823'%% dataset 20150823
			center = [-1.88966019628451	9.06825070793825	-8.74126683666250];
			initial_nor_R = .4 % [mm]
			final_nor_R = .6 % [mm]
			initial_tan_R = .3 % [mm]
			final_tan_R = .3 % [mm]				
			frame = 11;
		case '20150830'%% dataset 20150830
			center = [9.7667 5.7960 -3.0269];
			initial_nor_R = .3 % [mm]
			final_nor_R = .8 % [mm]
			initial_tan_R = .2 % [mm]
			final_tan_R = .4 % [mm]				
			frame = 11;
		case '20150913'%% dataset 20150913
			center = [-2.698042507113268,11.490590391752633,-4.684133112167247];
			initial_nor_R = .6 % [mm]
			final_nor_R = .8 % [mm]
			initial_tan_R = .2 % [mm]
			final_tan_R = .5 % [mm]
			frame = 10;
		case '20150920'%% dataset 20150920
			center = [ 3.4982   10.8276   -4.4754];
			initial_nor_R = .5 % [mm]
			final_nor_R = .7 % [mm]
			initial_tan_R = .2 % [mm]
			final_tan_R = .3 % [mm]
			frame = 10;
		case '20151008'
			center = [5.3413 7.4723 1.7974];
			initial_nor_R = .3 % [mm]
			final_nor_R = .8 % [mm]
			initial_tan_R = .3 % [mm]
			final_tan_R = .3 % [mm]
			frame = 11;
		case '20151018'%%  Late Systolic
			center = [11.2793 5.6288 8.7150];
			initial_nor_R = .2 % [mm]
			final_nor_R = .6 % [mm]
			initial_tan_R = .3 % [mm]
			final_tan_R = .3 % [mm]
			frame = 11;
		% case '20151018'%% Early Systolic
			% center = [11.2793 5.6288 8.7150];
			% initial_nor_R = .2 % [mm]
			% final_nor_R = .6 % [mm]
			% initial_tan_R = .2 % [mm]
			% final_tan_R = .3 % [mm]
		case '20151101'% seeds:4+3+3+3+5+4
			center = [-1.0529 5.9606 -6.4824];
			initial_nor_R = .6; % [mm]
			final_nor_R = .9; % [mm]
			initial_tan_R = .4; % [mm]
			final_tan_R = .4; % [mm]
			frame = 9;
		% case '20151101'%% seeds:1+2+5+3+5+4
			% center = [-1.1265 6.0892 -6.4609];
			% initial_nor_R = .5; % [mm]
			% final_nor_R = 1; % [mm]
			% initial_tan_R = .1; % [mm]
			% final_tan_R = .3; % [mm]
		% case '20151101'%% seeds:3+3+5+3+5+2
			% center = [-1.3411 6.3238 -6.3429];
			% initial_nor_R = .5; % [mm]
			% final_nor_R = .9; % [mm]
			% initial_tan_R = .3; % [mm]
			% final_tan_R = .4; % [mm]
		% case '20151103'%% seeds:4+2+3+3+3
			% center = [-3.7081 4.6878 -3.2771]
		case '20151103'%% seeds:4+1+3+3+3
			center = [-3.7764 4.7990 -3.2689];
			initial_nor_R = .5; % [mm]
			final_nor_R = .7; % [mm]
			initial_tan_R = .2; % [mm]
			final_tan_R = .2; % [mm]
			% frame = ;
		case '20151104'%% seeds:5+3+3+3+3+3
			center = [2.9737 5.0854 -4.8901];
			initial_nor_R = .4; % [mm]
			final_nor_R = .6; % [mm]
			initial_tan_R = .2; % [mm]
			final_tan_R = .5; % [mm]
			frame = 11;
		end
	end
	
	% if nargout;varargout = {, }; end	
end
