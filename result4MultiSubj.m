function result4MultiSubj(varargin)
% Last Modified: 1:06 PM Thursday, October 15, 2015
% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)

% set currunt Dir.
import plugins.DENSE3D_Plugin_4CrescentOrgan.*


%% Define some Program Constants
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% fields = {'XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ', 'RR', 'RC', 'RL', 'CR', 'CC', 'CL', 'LR', 'LC', 'LL', 'p1', 'p2', 'p3', 'CLShearAngle', 'J'};
Marker = {'x','+','*','o','^','s','v','d','>','p','<','+','*','o','s','d','x'};
Color = {'b','r','k','g','y','m','c','b','r','k','g','y','m','c'};
% LineStyle = {'-',':','-.','--'};

% commandwindow;
%%% Analysis Type:
tmp = {'1.Average Strain Curves (all-in-one)','2.Strain curves plotted curve-wise','3.Global Strain curves plotted curve-wise Base, Mid-ventricular, Apex and the Whole LV','4.Global Torsion','5.Transmural trend','6.Figures in Thesis Format','7.Segmentational strain curves with errorbar(Matlab 2013b or below)','7.1.Regional strain curves with errorbar(Matlab 2013b or below)','7.2.Transmural strain curves with errorbar(Matlab 2013b or below)','7.3.Global strain curves at Base and Apex plotted with errorbar(Matlab 2013b or below)','7.4.Mid-ventricular Global strain curves with errorbar','8.Reproducibility assessment(Report the numbers of pairwise datasets)','8.1.Reproducibility assessment for global strains: Base VS Mid VS Apex','9.Two-Way ANOVA','10.Output raw data for statistical test + Std. Err. & Summary of Peak Strains','11.Global Peak strains of 3 SA slices: Base vs Mid vs Apex','','----------Control VS Isoproterenol----------','12.Control VS Isoproterenol: Strain curves plotted curve-wise','13.Control VS Isoproterenol: Average Strain Curves (all-in-one)','14. Table and Bar Graph of Global Peak strains of 3 SA slices: Base vs Mid vs Apex(Matlab 2013b or below)','15.Segmentational strain curves with errorbar(Matlab 2013b or below)','15.1.Regional strain curves with errorbar(Matlab 2013b or below)','15.2.Transmural strain curves with errorbar(Matlab 2013b or below)','15.3.Global strain curves with errorbar : Base vs Mid vs Apex (Matlab 2013b or below)','16.Two-sample t-test with unpooled variance & wilcoxon rank sum test','17.Output raw data for statistical test + Std. Err. & Summary of Peak Strains','Select the number of output option:'}; 
option = input(strcat([sprintf('%s\n',tmp{:})]));
% tmp = {'1.JPEG','2.FIG(compatible for Matlab)','Select the format of figures to be exported:'}; 
% figtype = input(strcat([sprintf('%s\n',tmp{:})]));

%%% Files selection:
% button = questdlg('Do you wish to use a User Interface to select Datasets of interest?','Options','No');
% if strcmp(button,'No')
button = input('Press any key for using a User Interface to select Datasets of interest.\nOtherwise leave it blank!\n','s');
if isempty(button)
	% Upper level of path:
	% path = strsplit(pwd, '\'); path = {path{1:end-1}}; path = strjoin(path, '\');
	tmp = dir(pwd); folders = {tmp.name}; folders = folders([tmp.isdir] == 1); folders = {folders{3:end}};
	tmp = arrayfun(@num2str, 1:numel(folders), 'UniformOutput', false);
	tmp = strcat(tmp, '.', folders);
	ind = input(strcat([sprintf('%s\n',tmp{:}),'Input a MATRIX of the numbers of Datasets of interest[use Space/Colon/Semicolon(;) to seperate]:']));
	ind = reshape(ind,[],1); ind = round(abs(ind)); ind(ind==0) = [];
	% Default input datasets:
	if isempty(ind)
		% if option == 8
		if floor(option) == 8
			% 1103 1101 0920 1008 1018
			ind = [9 19 8 18 5 15 6 16 7 17]
			% ind = [1:3, 5:6,8:9;11:13, 15:16,18:19]
			% tmp = [1:10;11:20]; ind = tmp(:);
		else
			ind = [11:20]
		end
	end
	if option >= 12
		nBaseline = numel(ind);
		tmp = input(strcat([sprintf('%s\n',tmp{:}),'Input a MATRIX of the numbers of Isoproterenol Datasets[use Space/Space/Semicolon(;) to seperate]:']));
		if isempty(tmp)
			tmp = [21:29]
		end
		ind = [ind, tmp];
	end
	folders = folders(ind);
	files = fullfile(pwd,folders,strcat('DENSE3Dobject','.mat'));
else
	clear files folders;
	ii = 0;
	while true
		[uifile, uipath, uipopup] = uigetfile({'*.mat', 'Select DENSE3D workspaces (*.mat)'},'Open',pwd);
		if ~uipopup %isequal(uifile,0) || isequal(uipath,0)
			break
		else
			ii = ii + 1;
			files{ii} = fullfile(uipath,uifile);
			tmp = fileparts(uipath);
			[~,folders{ii}] = fileparts(tmp);
		end
	end
end


%%% Read files:
nFiles = numel(files);
clear data;
% data = cellfun(@(x)load(x, 'DENSE3Dobject', '-mat'), files);
for ii = 1:nFiles
	if ~exist(files{ii},'file')
		% Variable names: files-DENSE3D objects V.S. Files-SA slices
		Files = DENSE2DobjectFileName(fullfile(pwd,folders{ii}),'SA')
		% Files = fullfile(pwd,folders{ii},Files);
		% Files = DENSE2DobjectFileName()
		dbstop if error
		filename = fullfile(pwd,folders{ii},['DENSE3Dobject']);
		DENSE3Dobject = DENSE3D(Files,'Layers',Layers,'Interp',Interp,'Regions',Regions,'Insertion','SegDistro','Gif',filename)%,'Output','Type','Linear','Constant',0,'Smooth',0)
			%   The Following parameter/value pairs are allowed:
			%
			%   Output:     String, Path to a MAT file to store the results
			%   Layers:     Integer, Indicates how many transmural layers
			%               (Default = 5)
			%   Interp:     Integer, Amount of interpolation to use in z
			%               direction (Default = 1)
			%   Regions:    Integer, indicates how many elements we want
			%               circumferentially (Default = 64)
			%   Type:       Function Handle or String, Type of RBF
			%               interpolation to use (Default = 'Linear')
			%   Constant:   Scalar, RBF Constant (Default = 0)
			%   Smooth:     Scalar, RBF Smoothing Factor (Default = 0)
		% DENSE3Dobject = DENSE3D(Files)
		% DENSE3Dobject = DENSE3D(Files,'Layers',3,'Interp',0,'Regions',60,'Insertion','average','SegDistro','average')
		% save(strcat('DENSE3Dobject'),'DENSE3Dobject');
		save(filename,'DENSE3Dobject');
	end
	tmp = load(files{ii}); %, 'DENSE3Dobject', '-mat');
	try
		nFrames(ii) = tmp.DENSE3Dobject.Data(1).SequenceInfo(1, 1).CardiacNumberOfImages;
	catch
		nFrames(ii) = size(tmp.DENSE3Dobject.Mesh.regionalstrains.CC{1},2);
	end
	data{ii} = tmp.DENSE3Dobject.Mesh.regionalstrains;
	
	% switch option
		% case 4
	% torsion_global{ii} = mean(tmp.DENSE3Dobject.Mesh.regionalstrains.CLShearAngle{1},1);
	torsion_global{ii} = mean(tmp.DENSE3Dobject.Mesh.strains.CLShearAngle,1);
	
		% case 11
	for jj = 1:numel(fields)
		strain_SliceGlobal{ii}.(fields{jj}) = mean(cat(1,data{ii}.(fields{jj}){:}),1);
		strain_BaseGlobal{ii}.(fields{jj}) = mean(tmp.DENSE3Dobject.Mesh.strains.(fields{jj})(tmp.DENSE3Dobject.Mesh.longInd>2,:),1);
		strain_ApexGlobal{ii}.(fields{jj}) = mean(tmp.DENSE3Dobject.Mesh.strains.(fields{jj})(tmp.DENSE3Dobject.Mesh.longInd<2,:),1);
		strain_Global{ii}.(fields{jj}) = mean(tmp.DENSE3Dobject.Mesh.strains.(fields{jj}),1);		
		%{ 
		% 4th layer for global results:
		strain_SliceGlobal{ii}.(fields{jj}) = mean(data{ii}.(fields{jj}){1,4},1);
		strain_BaseGlobal{ii}.(fields{jj}) = mean(tmp.DENSE3Dobject.Mesh.strains.(fields{jj})(tmp.DENSE3Dobject.Mesh.absind==1,:),1);
		temp = max(tmp.DENSE3Dobject.Mesh.absind);
		strain_ApexGlobal{ii}.(fields{jj}) = mean(tmp.DENSE3Dobject.Mesh.strains.(fields{jj})(tmp.DENSE3Dobject.Mesh.absind==temp,:),1);
		 %}
	end
end


% Validate inputed datasets:
[ind, ~] = size(data{1}.(fields{1}));
[ind2, ~] = size(data{1}.(fields{1}){1,1});
for ii = 1:nFiles
	[nSlices, nLayers] = size(data{ii}.(fields{1}));
	if nSlices ~= ind
		error(strcat('The number of Slices for each Inputed Datasets is different:',sprintf(' %d',nSlices),' NOT EQUAL TO',sprintf(' %d',ind)));
		return	
	end
	[nRegions, ~] = size(data{ii}.(fields{1}){1,1});
	if nRegions > 6 || nRegions ~= ind2
		error(strcat('The number of regions for Input #',sprintf('%d',ii),' is invalid:',sprintf('%d',nRegions),'----Only support that myocardium is divided into 4 or 6 regions.'));
		return
	end
	if nLayers > 4
		if mod((nLayers-1),2)~=0
			layer = [1,ceil((nLayers-1)/2),nLayers-1,nLayers];
			for jj = 1:numel(fields)					
				for slice = 1:nSlices		
					data{ii}.(fields{jj}) = data{ii}.(fields{jj})(slice,layer);
				end
			end
		else
			clear tmp2;
			for jj = 1:numel(fields)					
				for slice = 1:nSlices		
					tmp = vertcat(data{ii}.(fields{jj}){slice,(nLayers-1)/2},data{ii}.(fields{jj}){slice,(nLayers-1)/2+1});
					tmp2{slice,1} = data{ii}.(fields{jj}){slice, 1};
					tmp2{slice,2}(1,:) = mean(tmp(1:4:nLayers,:));
					tmp2{slice,2}(2,:) = mean(tmp(2:4:nLayers,:));
					tmp2{slice,2}(3,:) = mean(tmp(3:4:nLayers,:));
					tmp2{slice,2}(4,:) = mean(tmp(4:4:nLayers,:));
					tmp2{slice,3} = data{ii}.(fields{jj}){slice, nLayers-1};
					tmp2{slice,4} = data{ii}.(fields{jj}){slice, nLayers};
				end
				data{ii}.(fields{jj}) = tmp2;
			end
		end
	end
end

%% Define some Program Constants
if nRegions == 4
	regions = {'Septal','Inferior','Lateral','Anterior'};
elseif nRegions == 6
	regions = {'Anteroseptal','Inferoseptal','Inferior','Inferolateral','Anterolateral','Anterior'};
else
	error(strcat('The number of regions for Input #',sprintf('%d',ii),' is invalid:',sprintf('%d',nRegions),'----Only support that myocardium is divided into 4 or 6 regions.'));
	return
end
if nLayers == 1
	transmural = {' Myocardium'};
elseif nLayers == 4
	transmural = {' Sub-endocardium',' Mid-wall',' Sub-epicardium',' Myocardium'};	
else
	error(strcat('The number of layers for Input #',sprintf('%d',ii),' is invalid:',sprintf('%d',nRegions),'----Only support that myocardium is divided into 1 or 4 layers.'));
	return
end

% Normalized:
minFrame = min(nFrames);
queryPoints = [1/minFrame:1/minFrame:1]*100;
clear normalData;
for ii = 1:nFiles
	if nFrames(ii) == minFrame
		% normalDataL{ii} = data{ii};
		% normalDataS{ii} = data{ii};
		normalData{ii} = data{ii};
		torsionH(ii,:) = torsion_global{ii};
		%{ 
		%% Time curves for basal & apical global strains:
		for jj = 1:numel(fields)					
			normalSliceGlobalData.(fields{jj})(ii,:) = strain_SliceGlobal{ii}.(fields{jj});
			normalBaseGlobalData.(fields{jj})(ii,:) = strain_BaseGlobal{ii}.(fields{jj});
			normalApexGlobalData.(fields{jj})(ii,:) = strain_ApexGlobal{ii}.(fields{jj});
			normalGlobalData.(fields{jj})(ii,:) = strain_Global{ii}.(fields{jj});
		end
		 %}
	else
		samplePoints = [1/nFrames(ii):1/nFrames(ii):1]*100;
		for jj = 1:numel(fields)					
			for slice = 1:nSlices
				for layer = 1:4
					for region = 1:nRegions	
						sampleValues = data{ii}.(fields{jj}){slice,layer}(region,:);
						% normalData{ii}.(fields{jj}){slice,layer}(region,:) = interp1(samplePoints, sampleValues, queryPoints);
						normalData{ii}.(fields{jj}){slice,layer}(region,:) = spline(samplePoints, sampleValues, queryPoints);
						% normalData{ii}.(fields{jj}){slice,layer}(region,:) = pchip(samplePoints, sampleValues, queryPoints);
					end
				end
			end
			
			%% Time curves for basal & apical global strains:
			normalSliceGlobalData.(fields{jj})(ii,:) = spline(samplePoints, strain_SliceGlobal{ii}.(fields{jj}), queryPoints);
			normalBaseGlobalData.(fields{jj})(ii,:) = spline(samplePoints, strain_BaseGlobal{ii}.(fields{jj}), queryPoints);
			normalApexGlobalData.(fields{jj})(ii,:) = spline(samplePoints, strain_ApexGlobal{ii}.(fields{jj}), queryPoints);
			normalGlobalData.(fields{jj})(ii,:) = spline(samplePoints, strain_Global{ii}.(fields{jj}), queryPoints);
		end
		torsionH(ii,:) = spline(samplePoints, torsion_global{ii}, queryPoints);
	end
end

% Interpolation Validation:
%{
%% Interpolate points exactly in a line connecting  nearby 2 points: conservative
figure
plot([1/nFrames(ii):1/nFrames(ii):1],data{ii}.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','b','Marker','x','MarkerSize',5); 
hold on;
plot([1/minFrame:1/minFrame:1],normalDataL{ii}.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color', 'r','Marker', '+', 'MarkerSize',5); 
%% Preserve Max value:
figure %bold
plot([1/nFrames(ii):1/nFrames(ii):1],data{ii}.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','b','Marker','x','MarkerSize',5); 
hold on;
plot([1/minFrame:1/minFrame:1],normalDataS{ii}.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color', 'r','Marker', '+', 'MarkerSize',5); 
%% Interpolate points closely in a line connecting  nearby 2 points:
 %}
if nFrames(ii) ~= minFrame
	figure
	plot([1/nFrames(ii):1/nFrames(ii):1],data{ii}.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','b','Marker','x','MarkerSize',5); 
	hold on;
	plot([1/minFrame:1/minFrame:1],normalData{ii}.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color', 'r','Marker', '+', 'MarkerSize',5); 
	title('Validate Interpolation!');
end

% Leading name of figures:
if numel(folders) < 4
	fname = strjoin(folders,'VS');
else
	fname = strcat(num2str(numel(folders)),'subjects');
end

%%% 3D Regional and Transmural Strain
switch option
% switch floor(option)
case 1
% Average Strain Curves (all-in-one):
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 1]);
		for region = 1:nRegions	
			for layer = 1:nLayers
				tmp = [];
				for ii = 1:nFiles
					tmp = [tmp; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
				end
				% E.(fields{jj}){slice,layer}(region,:)
				tmp2 = mean(tmp);
				% stderr.(fields{jj}){slice,layer}(region,:)
				tmp3 = sem(tmp);
				ave.(fields{jj}){slice,layer}(region,:) = tmp2;
				%% Tpeak is for sign preserve!!!
				[~,Tpeak] = max(abs(tmp2));
				Epeak.(fields{jj})(region,layer) = tmp2(Tpeak);
				stderr.(fields{jj})(region,layer) = tmp3(Tpeak);
				CoV = tmp3(Tpeak)/abs(tmp2(Tpeak))*100;
				
				tEpeak{nLayers*(jj-1)+layer,region} = sprintf(['%0.2f',177,'%0.2f'], [tmp2(Tpeak);tmp3(Tpeak)]);
				
				subplot(nRegions,nLayers,nLayers*(region-1)+layer);
				title(strcat(regions{region},transmural{layer},'-CoV@peak:',num2str(CoV),'%'));
				% xlabel('Cardiac Cycle(%)');
				xlabel(strcat('Cardiac Cycle(%):Peak=',num2str(round(queryPoints(Tpeak))),'%'));
				% ylabel(strcat(regions{region},'-',transmural{layer}));
				% ylabel(strcat('E',lower(fields{jj})));
				ylabel(strcat('E',lower(fields{jj}),':Peak@',tEpeak{nLayers*(jj-1)+layer,region}));
				grid on;
				hold on;
				plot(queryPoints,tmp2,'LineStyle',':','LineWidth',1,'color','b','Marker','x','MarkerSize', 5);
				plot(queryPoints,tmp2-tmp3,'LineStyle',':','LineWidth',0.5,'color','^','Marker','+','MarkerSize', 3);
				plot(queryPoints,tmp2+tmp3,'LineStyle',':','LineWidth',0.5,'color','v','Marker','+','MarkerSize', 3);
				set(gca,'XMinorTick','on','YMinorTick','on');
			end
		end
		axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
		text(0.5, 1,{strcat('slice#',num2str(slice),'-Peak Strain E',lower(fields{jj})),strcat('Distrib. of Circumf. Elemt:',', #Segments/slice:',', Transmural Interp:',', Z Interp:'),'Dark:Mean VS Light:±Standard Deviation'},'HorizontalAlignment','center','VerticalAlignment', 'top');
		hgexport(gcf,strcat('StrainCurves_AllInOne_',fname,'_E',lower(fields{jj}),num2str(slice)),hgexport('factorystyle'), 'Format', 'jpeg');
		close gcf;
	end
end
% numel(fields)*nLayers
% table = [{'Components','PeakStrain','Septal','Inferior','Lateral','Anterior'};[reshape([fields;repmat({''},1,numel(fields));repmat(repmat({''},1,numel(fields)),2,1)],[],1),repmat({'Endo','Mid','Epi','Myocardium'}',numel(fields),1),tEpeak]];
% xlswrite(strcat('Summary_PeakStrain_',fname),table);
save(strcat('Var_PeakStrain_',fname),'Epeak','stderr','CoV');
		
case 2
% Strain curves plotted curve-wise:
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 1]);
		% title(strcat('Peak Strain E',lower(fields{jj})));
		for region = 1:nRegions	
			for layer = 1:nLayers
				% subplot(nRegions,nLayers,nLayers*(region-1)+layer);
				subplot(nRegions,nLayers,nLayers*(region-1)+layer);
				title(strcat(regions{region},transmural{layer}));
				hold on;
				for ii = 1:nFiles
					plot(normalData{ii}.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color', Color{ii},'Marker', Marker{ii},'MarkerSize',5); %,'MarkerEdgeColor','b');
				end
				set(gca,'XMinorTick','on','YMinorTick','on');
				% axis equal;
				xlabel('Cardiac Cycle(frames)');
				% ylabel(strcat('Peak Strain E',lower(fields{jj}),' at',regions{region},transmural{layer}));
				% ylabel(strcat(regions{region},transmural{layer}));
				ylabel(strcat('E',lower(fields{jj})));
				grid on;
			end
		end
		% hold off;
		axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
		text(0.5, 1,[strcat(' slice#',num2str(slice),'-Peak Strain E',lower(fields{jj})),strcat('Distrib. of Circumf. Elemt:',', #Segments/slice:',', Transmural Interp:',', Z Interp:'),strcat(Color(1:nFiles), Marker(1:nFiles),'.', folders)],'HorizontalAlignment','center','VerticalAlignment', 'top');
		%saveas(gcf,strcat(fname,'_',num2str(slice),'_E',lower(fields{jj})),'jpg');
		tmp = strcat(fname,'_','_E',lower(fields{jj}),num2str(slice));
		hgexport(gcf,tmp,hgexport('factorystyle'), 'Format', 'jpeg');
		hgsave(gcf,tmp);
		close gcf;
	end
end

case 3
% Global Strain curves plotted curve-wise at Base, Mid-ventricular, Apex and the Whole LV:
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 1]);

		subplot(2,2,1);
		title('Apex');
		hold on;
		for ii = 1:nFiles
			plot(normalApexGlobalData.(fields{jj})(ii,:),'LineWidth',1,'color', Color{ii},'Marker', Marker{ii},'MarkerSize',5); %,'MarkerEdgeColor','b');
		end
		set(gca,'XMinorTick','on','YMinorTick','on');
		% axis equal;
		xlabel('Cardiac Cycle(frames)');
		% ylabel(strcat('Peak Strain E',lower(fields{jj}),' at',regions{region},transmural{layer}));
		% ylabel(strcat(regions{region},transmural{layer}));
		ylabel(strcat('E',lower(fields{jj})));
		grid on;
		
		subplot(2,2,2);
		title('Base');
		hold on;
		for ii = 1:nFiles
			plot(normalBaseGlobalData.(fields{jj})(ii,:),'LineWidth',1,'color', Color{ii},'Marker', Marker{ii},'MarkerSize',5); %,'MarkerEdgeColor','b');
		end
		set(gca,'XMinorTick','on','YMinorTick','on');
		% axis equal;
		xlabel('Cardiac Cycle(frames)');
		% ylabel(strcat('Peak Strain E',lower(fields{jj}),' at',regions{region},transmural{layer}));
		% ylabel(strcat(regions{region},transmural{layer}));
		ylabel(strcat('E',lower(fields{jj})));
		grid on;

		subplot(2,2,3);
		title('Mid-ventricular');
		hold on;
		for ii = 1:nFiles
			plot(normalSliceGlobalData.(fields{jj})(ii,:),'LineWidth',1,'color', Color{ii},'Marker', Marker{ii},'MarkerSize',5); %,'MarkerEdgeColor','b');
		end
		set(gca,'XMinorTick','on','YMinorTick','on');
		% axis equal;
		xlabel('Cardiac Cycle(frames)');
		% ylabel(strcat('Peak Strain E',lower(fields{jj}),' at',regions{region},transmural{layer}));
		% ylabel(strcat(regions{region},transmural{layer}));
		ylabel(strcat('E',lower(fields{jj})));
		grid on;

		subplot(2,2,4);
		title('Whole LV');
		hold on;
		for ii = 1:nFiles
			plot(normalGlobalData.(fields{jj})(ii,:),'LineWidth',1,'color', Color{ii},'Marker', Marker{ii},'MarkerSize',5); %,'MarkerEdgeColor','b');
		end
		set(gca,'XMinorTick','on','YMinorTick','on');
		% axis equal;
		xlabel('Cardiac Cycle(frames)');
		% ylabel(strcat('Peak Strain E',lower(fields{jj}),' at',regions{region},transmural{layer}));
		% ylabel(strcat(regions{region},transmural{layer}));
		ylabel(strcat('E',lower(fields{jj})));
		grid on;
		
		axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
		text(0.5, 1,[strcat(' slice#',num2str(slice),'-Peak Strain E',lower(fields{jj})),strcat('Distrib. of Circumf. Elemt:',', #Segments/slice:',', Transmural Interp:',', Z Interp:'),strcat(Color(1:nFiles), Marker(1:nFiles),'.', folders)],'HorizontalAlignment','center','VerticalAlignment', 'top');
		%saveas(gcf,strcat(fname,'_',num2str(slice),'_E',lower(fields{jj})),'jpg');
		tmp = strcat(fname,'Global_','_E',lower(fields{jj}),num2str(slice));
		hgexport(gcf,tmp,hgexport('factorystyle'), 'Format', 'jpeg');
		hgsave(gcf,tmp);
		close gcf;
	end
end

case 4
% Global Torsion:
clear tmp;
figure('units','normalized','outerposition',[0 0 1 1]);
% title(strcat('Peak Strain E',lower(fields{jj})));
% axis equal;
hold on;
for ii = 1:nFiles
	tmp{ii} = plot(queryPoints,torsionH(ii,:),'LineWidth',1,'color', Color{ii},'Marker', Marker{ii},'MarkerSize',5); %,'MarkerEdgeColor','b');
end
xlabel('Cardiac Cycle(%)');
ylabel('Circumferential-longitudinal shear angle(°/cm^2)');
grid on;
set(gca,'XMinorTick','on','YMinorTick','on');
legend([tmp{:}], folders);
hgexport(gcf,strcat('Torsion_'),hgexport('factorystyle'), 'Format', 'jpeg');
close gcf;
% Mean+Std. Dev.:
torsionAve = mean(torsionH);
torsionStd = sem(torsionH);
figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
plot(queryPoints,torsionAve,'LineWidth',1,'color','b','Marker','x','MarkerSize', 5);
plot(queryPoints,torsionAve-torsionStd,'LineWidth',0.5,'color','c','Marker','+','MarkerSize', 3);
plot(queryPoints,torsionAve+torsionStd,'LineWidth',0.5,'color','c','Marker','+','MarkerSize', 3);
xlabel('Cardiac Cycle(%)');
ylabel('Circumferential-longitudinal shear angle(°/cm^2)');
grid on;
set(gca,'XMinorTick','on','YMinorTick','on');
hgexport(gcf,strcat('TorsionStd_'),hgexport('factorystyle'), 'Format', 'jpeg');
close gcf;

case 5 
% Transmural trend:
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
frame = 9;
for jj = 1 : 3
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 0.5]);
		for region = 1:nRegions	
			tmp = []; tmp2 = [];
			for layer = 1:3
				tmp = [tmp,ave.(fields{jj}){slice,layer}(region,frame)];
				tmp2 = [tmp2,stderr.(fields{jj}){slice,layer}(region,frame)];
			end
			subplot(1,nRegions,region);
			title(strcat(regions{region},' segments'));
			xlabel('Sub-endocardium        Mid-wall        Sub-epicardium');
			set(gca,'XTickLabel', {'Sub-endocardium','Mid-wall','Sub-epicardium'});
			ylabel(strcat('E',lower(fields{jj})));
			grid on;
			hold on;
			plot(tmp,'LineWidth',1,'color','b','Marker','x','MarkerSize', 5);
			plot(tmp-tmp2,'LineWidth',0.5,'color','c','Marker','+','MarkerSize', 3);
			plot(tmp+tmp2,'LineWidth',0.5,'color','c','Marker','+','MarkerSize', 3);
			% set(gca,'XMinorTick','off','YMinorTick','on');
			set(gca,'XTick',[1 2 3],'YMinorTick','on');
		end
		hgexport(gcf,strcat('TransmuralTrend_','_E',lower(fields{jj}),num2str(slice)),hgexport('factorystyle'), 'Format', 'jpeg');	
		close gcf;
	end
end

case 6
% Figures in Thesis Format: 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 1]);
		for region = 1:nRegions	
			for layer = 1:nLayers
				subplot(nRegions,nLayers,nLayers*(region-1)+layer);
				title(strcat(regions{region},transmural{layer}));
				xlabel('Cardiac Cycle(%)');
				ylabel(strcat('E',lower(fields{jj})));
				grid on;
				hold on;
				tmp = [];
				for ii = 1:nFiles
					tmp = [tmp; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
				end
				tmp2 = mean(tmp);
				tmp3 = sem(tmp);
				ave.(fields{jj}){slice,layer}(region,:) = tmp2;
				stderr.(fields{jj}){slice,layer}(region,:) = tmp3;
				plot(queryPoints,tmp2,'LineWidth',1,'color','b','Marker','x','MarkerSize', 5);
				plot(queryPoints,tmp2-tmp3,'LineWidth',0.5,'color','c','Marker','+','MarkerSize', 3);
				plot(queryPoints,tmp2+tmp3,'LineWidth',0.5,'color','c','Marker','+','MarkerSize', 3);
				set(gca,'XMinorTick','on','YMinorTick','on');
				hgexport(gcf,strcat('Strain','_E',lower(fields{jj}),num2str(slice)),hgexport('factorystyle'), 'Format', 'jpeg');
				close gcf;
			end
		end
	end
end

case 7
% Segmentational strain curves with errorbar: 
fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		for region = 1:nRegions	
			figure('units','normalized','outerposition',[0 0 1 1]);
			for layer = 1:nLayers-1
				tmp = [];
				for ii = 1:nFiles
					tmp = [tmp; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
				end
				tmp2 = mean(tmp);
				tmp3 = sem(tmp);
				ave.(fields{jj}){slice,layer}(region,:) = tmp2;
				stderr.(fields{jj}){slice,layer}(region,:) = tmp3;
				hb{layer} = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{layer+3},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
				% Adjust error bar width: only work in 2013b
				errorbar_tick(hb{layer}, 150);%140-(-1)^layer*20
				hold on;
			end
			% title(strcat(regions{region},transmural{layer}));
			xlabel('Cardiac Cycle(%)');
			ylabel(strcat('E',lower(fields{jj})));
			set(gca,'XMinorTick','on','YMinorTick','on');
			grid on;
			legend(transmural{1:nLayers-1});
			tmp = strcat('StrainCurves_Segmentational_',fname,'_E',lower(fields{jj}),num2str(slice),'_',regions(region));
			hgexport(gcf,tmp{:},hgexport('factorystyle'), 'Format', 'jpeg');
			hgsave(gcf,tmp{:});
			close gcf;
		end
	end
end
save(strcat('Var_errorbar_',fname),'aveBase','errBase','aveIso','errIso');

case 7.1
% Regional strain curves with errorbar: 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
layer = nLayers;
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 1]);
		for region = 1:nRegions	
			tmp = [];
			for ii = 1:nFiles
				tmp = [tmp; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
			end
			tmp2 = mean(tmp);
			tmp3 = sem(tmp);
			ave.(fields{jj}){slice,layer}(region,:) = tmp2;
			stderr.(fields{jj}){slice,layer}(region,:) = tmp3;
			hb = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{region+3},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
			% Adjust error bar width: only work in 2013b
			errorbar_tick(hb, 150);%140-(-1)^layer*20
			hold on;
		end
		% title(strcat(regions{region},transmural{layer}));
		xlabel('Cardiac Cycle(%)');
		ylabel(strcat('E',lower(fields{jj})));
		set(gca,'XMinorTick','on','YMinorTick','on');
		grid on;
		legend(regions);
		tmp = strcat('StrainCurves_Regional_',fname,'_E',lower(fields{jj}),num2str(slice));
		hgexport(gcf,tmp,hgexport('factorystyle'), 'Format', 'jpeg');
		hgsave(gcf,tmp);
		close gcf;
	end
end
save(strcat('Var_errorbar_Regional_',fname),'ave','stderr');

case 7.2
% Transmural strain curves with errorbar: 
fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 1]);
		for layer = 1:nLayers-1
			tmp = [];
			for ii = 1:nFiles
				tmp = [tmp; mean(normalData{ii}.(fields{jj}){slice,layer})];
			end
			tmp2 = mean(tmp);
			tmp3 = sem(tmp);
			ave.(fields{jj}){slice,layer} = tmp2;
			stderr.(fields{jj}){slice,layer} = tmp3;
			hb{layer} = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{layer+3},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
			% Adjust error bar width: only work in 2013b
			errorbar_tick(hb{layer}, 150);%140-(-1)^layer*20
			hold on;
		end
		% title(strcat(regions{region},transmural{layer}));
		xlabel('Cardiac Cycle(%)');
		ylabel(strcat('E',lower(fields{jj})));
		set(gca,'XMinorTick','on','YMinorTick','on');
		grid on;
		legend(transmural{1:nLayers-1});
		tmp = strcat('StrainCurves_Transmural_',fname,'_E',lower(fields{jj}),num2str(slice));
		hgexport(gcf,tmp{:},hgexport('factorystyle'), 'Format', 'jpeg');
		hgsave(gcf,tmp{:});
		close gcf;
	end
end
save(strcat('Var_errorbar_Transmural_',fname),'aveBase','errBase','aveIso','errIso');

case 7.3
% Global strain curves at Base and Apex in Paper Format: 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	figure('units','normalized','outerposition',[0 0 1 1]);
	tmp = [];

	tmp2 = mean(normalBaseGlobalData.(fields{jj}));
	tmp3 = sem(normalBaseGlobalData.(fields{jj}));
	hb{1} = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{4},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
	% Adjust error bar width: only work in 2013b
	errorbar_tick(hb{1}, 150);
	hold on;

	tmp2 = mean(normalApexGlobalData.(fields{jj}));
	tmp3 = sem(normalApexGlobalData.(fields{jj}));
	hb{2} = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{5},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
	% Adjust error bar width: only work in 2013b
	errorbar_tick(hb{2}, 150);
	
	% title(strcat(regions{region},transmural{layer}));
	xlabel('Cardiac Cycle(%)');
	ylabel(strcat('E',lower(fields{jj})));
	set(gca,'XMinorTick','on','YMinorTick','on');
	grid on;
	legend({' Base',' Apex'});
	tmp = strcat('StrainCurves_Global_BASEvsAPEX_',fname,'_E',lower(fields{jj}));
	hgexport(gcf,tmp,hgexport('factorystyle'), 'Format', 'jpeg');
	hgsave(gcf,tmp);
	close gcf;
end


case 7.4
% Mid-ventricular Global strain curves with errorbar: 
fields = {'RR', 'RC', 'RL'};
for jj = 1:numel(fields)					
	figure('units','normalized','outerposition',[0 0 1 1]);
		tmp2 = mean(normalSliceGlobalData.(fields{jj}));
		tmp3 = sem(normalSliceGlobalData.(fields{jj}));
		hb = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{1},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
		% Adjust error bar width: only work in 2013b
		errorbar_tick(hb, 150);
		hold on;
	% title(strcat(regions{region},transmural{layer}));
	xlabel('Cardiac Cycle(%)');
	ylabel(strcat('E',lower(fields{jj})));
	set(gca,'XMinorTick','on','YMinorTick','on');
	grid on;
	tmp = strcat('StrainCurves_Global_Mid_','_','_E',lower(fields{jj}));
	hgexport(gcf,tmp{:},hgexport('factorystyle'), 'Format', 'jpeg');
	hgsave(gcf,tmp{:});
	close gcf;
end

case 8
% Reproducibility assessment: [1 12 2 13 3 14 4 15 5 16 6 17 7 18 8 19 9 20 10 21] [8 19 7 18 4 15 5 16 6 17]
% 1101 1018 0913 0920 1008
if mod(nFiles,2)
	error('The number of datasets to be compared should be even');
	return
end
%{ 
if option = 8.1
strain_SliceGlobal = strain_BaseGlobal;
elseif option = 8.2
strain_SliceGlobal = strain_ApexGlobal;
end
 %}
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% pick my own bin locations
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		for layer = 1:nLayers
			for region = 1:nRegions	
				for ii = 1:nFiles
					[~,Tpeak] = max(abs(data{ii}.(fields{jj}){slice,layer}(region,:)));
					Epeak.(fields{jj}){slice,layer}(region,ii) = data{ii}.(fields{jj}){slice,layer}(region,Tpeak);
				end
				tmp = [Epeak.(fields{jj}){slice,layer}(region,1:2:nFiles);Epeak.(fields{jj}){slice,layer}(region,2:2:nFiles)];
				CoV.(fields{jj}){slice}(region,layer) = sum(std(tmp))/abs(sum(mean(tmp)))*100;
				tCoV(nLayers*(jj-1)+layer,region) = CoV.(fields{jj}){slice}(region,layer);
			end
			
			%%%% The fifth section of CoV:Endo -> Mid -> Epi -> Global
			for ii = 1:nFiles
				tmp3 = mean(data{ii}.(fields{jj}){slice,layer}(:,:));
				[~,Tpeak] = max(abs(tmp3));
				Epeak.(fields{jj}){slice,layer}(nRegions+1,ii) = tmp3(Tpeak);
			end
			tmp = [Epeak.(fields{jj}){slice,layer}(nRegions+1,1:2:nFiles);Epeak.(fields{jj}){slice,layer}(nRegions+1,2:2:nFiles)];
			CoV.(fields{jj}){slice}(nRegions+1,layer) = sum(std(tmp))/abs(sum(mean(tmp)))*100;
			tCoV(nLayers*(jj-1)+layer,nRegions+1) = CoV.(fields{jj}){slice}(nRegions+1,layer);
		end
		
		%{ 
		figure('units','normalized','outerposition',[0 0 1 1]);		
		bar(hCoV);
		title(strcat(regions{region},transmural{layer},'Inter-Oberver LV Strains'));
		xlabel('LV Segments');
		set(gca,'XTickLabel', regions);
		ylabel('Coefficient of Variation(%)');
		set(gca,'XMinorTick','off','YMinorTick','on');
		grid on;
		legend(transmural);
		applyhatch_pluscolor(gcf, '\/x-');
		% applyhatch_pluscolor(gcf, '\/x-', 0, [1 0 1 0]);, jet(4)
		hgexport(gcf,strcat('CoV_',fname,'_E',lower(fields{jj}),num2str(slice)),hgexport('factorystyle'), 'Format', 'jpeg');
		hgsave(gcf,strcat('CoV_',fname,'_E',lower(fields{jj}),num2str(slice)));
		close gcf;
		close gcf;
		table = [{'CoV','Endo','Mid','Epi','Myocardium'};[{'Septal','Inferior','Lateral','Anterior'}',num2cell(hCoV)]];
		xlswrite(strcat('CoV_',fname),table, strcat('E_',lower(fields{jj})));
		%}
		figure('units','normalized','outerposition',[0 0 1 1]);		
		bar(CoV.(fields{jj}){slice});
		title(strcat(regions{region},transmural{layer},' Inter-Oberver LV Strains'));
		xlabel('LV Segments');
		set(gca,'XTickLabel', {regions{:},'Global'});
		ylabel('Coefficient of Variation(%)');
		set(gca,'XMinorTick','off','YMinorTick','on');
		grid on;
		legend(transmural);
		applyhatch_pluscolor(gcf, '\/x-');
		% set(gcf,'Visible','on');
		hgexport(gcf,strcat('CoV_',fname,'_E',lower(fields{jj}),num2str(slice)),hgexport('factorystyle'), 'Format', 'jpeg');
		hgsave(gcf,strcat('CoV_',fname,'_E',lower(fields{jj}),num2str(slice)));
		close gcf;
		close gcf;
		table = [[{'CoV'},transmural];[[regions,{'Global'}]',num2cell(CoV.(fields{jj}){slice})]];
		xlswrite(strcat('Summary_CoV_',fname),table, strcat('E_',lower(fields{jj})));
	end
end
table = [[{'Components','CoV'},regions,{'Global'}];[reshape([fields;repmat(repmat({''},1,numel(fields)),nLayers-1,1)],[],1),repmat(transmural',numel(fields),1),num2cell(tCoV)]];
xlswrite(strcat('Summary_CoV_AllInOne',fname),table,'Langragian Strains');
save(strcat('Var_CoV_',fname),'CoV','tCoV');

case 8.1
%% CoV for global strains: Base VS Mid VS Apex 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
clear hCoV;
for jj = 1:numel(fields)			
	hCoV.(fields{jj}){1,1} = 'Base';
	hCoV.(fields{jj}){1,2} = 'Mid';
	hCoV.(fields{jj}){1,3} = 'Apex';
	tmpBase3 = []; tmpMid3 = []; tmpApex3 = [];
	for ii = 1:nFiles
		[~,Tpeak] = max(abs(strain_BaseGlobal{ii}.(fields{jj})));
		tmpBase3{ii} = strain_BaseGlobal{ii}.(fields{jj})(Tpeak);
		[~,Tpeak] = max(abs(strain_SliceGlobal{ii}.(fields{jj})));
		tmpMid3{ii} = strain_SliceGlobal{ii}.(fields{jj})(Tpeak);
		[~,Tpeak] = max(abs(strain_ApexGlobal{ii}.(fields{jj})));
		tmpApex3{ii} = strain_ApexGlobal{ii}.(fields{jj})(Tpeak);
	end
	tmpBase = []; tmpBase1 = 0; tmpBase2 = 0; 
	tmpMid = []; tmpMid1 = 0; tmpMid2 = 0; 
	tmpApex = []; tmpApex1 = 0; tmpApex2 = 0; 
	for ii = 1:2:nFiles
		tmpBase = [tmpBase3{ii},tmpBase3{ii+1}];
		tmpBase1 = tmpBase1 + mean(tmpBase);
		tmpBase2 = tmpBase2 + std(tmpBase);
		tmpMid = [tmpMid3{ii},tmpMid3{ii+1}];
		tmpMid1 = tmpMid1 + mean(tmpMid);
		tmpMid2 = tmpMid2 + std(tmpMid);
		tmpApex = [tmpApex3{ii},tmpApex3{ii+1}];
		tmpApex1 = tmpApex1 + mean(tmpApex);
		tmpApex2 = tmpApex2 + std(tmpApex);
	end
	hCoV.(fields{jj}){2,1} = tmpBase2/abs(tmpBase1)*100;
	hCoV.(fields{jj}){2,2} = tmpMid2/abs(tmpMid1)*100;
	hCoV.(fields{jj}){2,3} = tmpApex2/abs(tmpApex1)*100;		
	xlswrite(strcat('CoV_baseVSmidVSapex_',fname),hCoV.(fields{jj}), strcat('E_',lower(fields{jj})));
end

case 9
%% Two-Way ANOVA: Columns-Transmural Layers; Rows-Circumferential Segments
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		ANOVAin = [];
		for region = 1:nRegions	
			clear Epeak;
			for ii = 1:nFiles
				for layer = 1:nLayers-1
					[~,Tpeak] = max(abs(data{ii}.(fields{jj}){slice,layer}(region,:)));
					Epeak(ii,layer) = data{ii}.(fields{jj}){slice,layer}(region,Tpeak);
				end
			end
		ANOVAin = [ANOVAin;Epeak];
		end
	end
	[p.(fields{jj}),tbl.(fields{jj}),stats.(fields{jj})] = anova2(ANOVAin,nFiles,'off');
	xlswrite(strcat('ttest_2ANOVA.xls'),tbl.(fields{jj}), strcat('E', lower(fields{jj})));
end

case 10
% Output raw data to Excel for statistical test + Std. Dev. & Summary of Peak Strains:
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
clear Epeak;
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		for layer = 1:nLayers
			for region = 1:nRegions	
				tmp = [];
				for ii = 1:nFiles
					[~,Tpeak] = max(abs(data{ii}.(fields{jj}){slice,layer}(region,:)));
					tmp = [tmp, data{ii}.(fields{jj}){slice,layer}(region,Tpeak)];
				end
				Epeak.(fields{jj}){region,layer} = tmp;
				xlswrite(strcat('Statistics_Raw_',fields{jj},'.xls'),tmp, strcat(regions{region}, transmural{layer}));

				tEpeak{nLayers*(jj-1)+layer,region} = sprintf(['%0.2f',177,'%0.2f'], [mean(tmp);sem(tmp)]);
				% xlswrite(strcat(fname,'_peak_strains_mean.xls'),mean(tmp), strcat('E',lower(fields{jj})));
				% xlswrite(strcat(fname,'_peak_strains_stdev.xls'),sem(tmp), strcat('E',lower(fields{jj})));
			end
			
			for ii = 1:nFiles
				tmp3 = mean(data{ii}.(fields{jj}){slice,layer}(:,:));
				[~,Tpeak] = max(abs(tmp3));
				Epeak.(fields{jj}){slice,layer}(nRegions+1,ii) = tmp3(Tpeak);					
			end			
			xlswrite(strcat('Statistics_Raw_',fields{jj},'.xls'),Epeak.(fields{jj}){slice,layer}(nRegions+1,:), strcat('Global', transmural{layer}));
			tEpeak{nLayers*(jj-1)+layer,nRegions+1} = sprintf(['%0.2f',177,'%0.2f'], [mean(Epeak.(fields{jj}){slice,layer}(nRegions+1,:));sem(Epeak.(fields{jj}){slice,layer}(nRegions+1,:))]);			
		end
	end
end
table = [{'Components','PeakStrain','Septal','Inferior','Lateral','Anterior','Global'};[reshape([fields;repmat({''},1,numel(fields));repmat(repmat({''},1,numel(fields)),2,1)],[],1),repmat({'Endo','Mid','Epi','Myocardium'}',numel(fields),1),tEpeak]];
xlswrite(strcat('Summary_PeakStrain_',fname),table);
save(strcat('Var_PeakStrain_',fname),'Epeak','tEpeak');


case 11
% Global Peak strains of 3 SA slices: Base vs Mid vs Apex
clear strain_GlobalPeak;
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)				
	tmpBase = []; tmpMid = []; tmpApex = [];
	for ii = 1:nFiles
		strain_GlobalPeak.(fields{jj}){1,4+ii} = num2str(ii);
		[~,Tpeak] = max(abs(strain_BaseGlobal{ii}.(fields{jj})));
		tmpBase(ii) = strain_BaseGlobal{ii}.(fields{jj})(Tpeak);
		strain_GlobalPeak.(fields{jj}){2,4+ii} = tmpBase(ii);
		[~,Tpeak] = max(abs(strain_SliceGlobal{ii}.(fields{jj})));
		tmpMid(ii) = strain_SliceGlobal{ii}.(fields{jj})(Tpeak);
		strain_GlobalPeak.(fields{jj}){3,4+ii} = tmpMid(ii);
		[~,Tpeak] = max(abs(strain_ApexGlobal{ii}.(fields{jj})));
		tmpApex(ii) = strain_ApexGlobal{ii}.(fields{jj})(Tpeak);
		strain_GlobalPeak.(fields{jj}){4,4+ii} = tmpApex(ii);
	end
	strain_GlobalPeak.(fields{jj}){1,2} = 'Mean';
	strain_GlobalPeak.(fields{jj}){1,3} = 'StdErr';
	strain_GlobalPeak.(fields{jj}){1,4} = 'StdDev';
	strain_GlobalPeak.(fields{jj}){2,1} = 'Base';
	strain_GlobalPeak.(fields{jj}){3,1} = 'Mid';
	strain_GlobalPeak.(fields{jj}){4,1} = 'Apex';
	strain_GlobalPeak.(fields{jj}){2,2} = mean(tmpBase);
	strain_GlobalPeak.(fields{jj}){2,3} = sem(tmpBase)
	strain_GlobalPeak.(fields{jj}){2,4} = std(tmpBase);
	strain_GlobalPeak.(fields{jj}){3,2} = mean(tmpMid);
	strain_GlobalPeak.(fields{jj}){3,3} = sem(tmpMid);
	strain_GlobalPeak.(fields{jj}){3,4} = std(tmpMid);
	strain_GlobalPeak.(fields{jj}){4,2} = mean(tmpApex);
	strain_GlobalPeak.(fields{jj}){4,3} = sem(tmpApex);
	strain_GlobalPeak.(fields{jj}){4,4} = std(tmpApex);
	xlswrite(strcat('GlobalPeakStrain_',fname),strain_GlobalPeak.(fields{jj}), strcat('E_',lower(fields{jj})));	
end
		
case 12
% Control VS Isoproterenol: Strain curves plotted curve-wise
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 1]);
		for region = 1:nRegions	
			for layer = 1:nLayers
				tmp = [];
				for ii = 1:nFiles
					tmp = [tmp; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
				end
				tmp2 = mean(tmp);
				tmp3 = sem(tmp);
				ave.(fields{jj}){slice,layer}(region,:) = tmp2;
				[~,Tpeak] = max(abs(tmp2));
				stderr.(fields{jj}){slice,layer}(region,:) = tmp3;
				CoV = tmp3(Tpeak)/tmp2(Tpeak)*100;
				subplot(nRegions,nLayers,nLayers*(region-1)+layer);
				title(strcat(regions{region},transmural{layer},':Peak@',num2str(CoV),'%'));
				% xlabel('Cardiac Cycle(%)');
				xlabel(strcat('Cardiac Cycle(%):Peak@',num2str(round(queryPoints(Tpeak))),'%'));
				% ylabel(strcat(regions{region},'-',transmural{layer}));
				% ylabel(strcat('E',lower(fields{jj})));
				ylabel(strcat('E',lower(fields{jj}),':Peak@',num2str(tmp2(Tpeak))));
				grid on;
				hold on;
				plot(queryPoints,tmp2,'LineWidth',1,'color','b','Marker','x','MarkerSize', 5);
				plot(queryPoints,tmp2-tmp3,'LineWidth',0.5,'color','c','Marker','+','MarkerSize', 3);
				plot(queryPoints,tmp2+tmp3,'LineWidth',0.5,'color','c','Marker','+','MarkerSize', 3);
				plot(queryPoints,normalData{1}.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color', 'k','Marker', '+','MarkerSize',5); %,'MarkerEdgeColor','b');
				set(gca,'XMinorTick','on','YMinorTick','on');
			end
		end
		axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
		text(0.5, 1,{strcat('slice#',num2str(slice),'-Peak Strain E',lower(fields{jj})),strcat('Distrib. of Circumf. Elemt:',', #Segments/slice:',', Transmural Interp:',', Z Interp:'),'Dark:Mean VS Light:±Standard Deviation'},'HorizontalAlignment','center','VerticalAlignment', 'top');
		hgexport(gcf,strcat('CONTROLvsISO','_E',lower(fields{jj}),num2str(slice)),hgexport('factorystyle'), 'Format', 'jpeg');
		close gcf;
	end
end

case 13
% Control VS Isoproterenol: Average Strain Curves (all-in-one)
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 1]);
		for region = 1:nRegions	
			for layer = 1:nLayers
				tmp = []; tmp2 = [];
				for ii = 1:nBaseline
					tmp = [tmp; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
				end
				for ii = nBaseline+1:nFiles
					tmp2 = [tmp2; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
				end
				aveBase.(fields{jj}){slice,layer}(region,:) = mean(tmp);
				errBase.(fields{jj}){slice,layer}(region,:) = sem(tmp);
				aveIso.(fields{jj}){slice,layer}(region,:) = mean(tmp2);
				errIso.(fields{jj}){slice,layer}(region,:) = sem(tmp2);
				
				%% Tpeak is for sign preserve!!!
				[~,Tpeak] = max(abs(aveIso.(fields{jj}){slice,layer}(region,:)));
				% Epeak.(fields{jj})(region,layer) = aveIso.(fields{jj}){slice,layer}(region,Tpeak);
				% stderr.(fields{jj})(region,layer) = errIso.(fields{jj}){slice,layer}(region,Tpeak);
				CoV = errIso.(fields{jj}){slice,layer}(region,Tpeak)/abs(aveIso.(fields{jj}){slice,layer}(region,Tpeak))*100;
				tEpeak{nLayers*(jj-1)+layer,region} = sprintf(['%0.2f',177,'%0.2f'], [aveIso.(fields{jj}){slice,layer}(region,Tpeak);errIso.(fields{jj}){slice,layer}(region,Tpeak)]);
							
				subplot(nRegions,nLayers,nLayers*(region-1)+layer);
				title(strcat(regions{region},transmural{layer},'-CoV@peak:',num2str(CoV),'%'));
				% xlabel('Cardiac Cycle(%)');
				xlabel(strcat('Cardiac Cycle(%):Peak=',num2str(round(queryPoints(Tpeak))),'%'));
				% ylabel(strcat(regions{region},'-',transmural{layer}));
				% ylabel(strcat('E',lower(fields{jj})));
				ylabel(strcat('E',lower(fields{jj}),':Peak@',tEpeak{nLayers*(jj-1)+layer,region}));
				grid on;
				hold on;
				plot(queryPoints,aveBase.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','b','Marker','x','MarkerSize', 5);
				plot(queryPoints,aveBase.(fields{jj}){slice,layer}(region,:)-errBase.(fields{jj}){slice,layer}(region,:),'LineStyle',':','LineWidth',0.5,'color','c','Marker','^','MarkerSize', 3);
				plot(queryPoints,aveBase.(fields{jj}){slice,layer}(region,:)+errBase.(fields{jj}){slice,layer}(region,:),'LineStyle',':','LineWidth',0.5,'color','c','Marker','v','MarkerSize', 3);
				plot(queryPoints,aveIso.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','k','Marker','*','MarkerSize', 5);
				plot(queryPoints,aveIso.(fields{jj}){slice,layer}(region,:)-errIso.(fields{jj}){slice,layer}(region,:),'LineStyle',':','LineWidth',0.5,'color','g','Marker','^','MarkerSize', 3);
				plot(queryPoints,aveIso.(fields{jj}){slice,layer}(region,:)+errIso.(fields{jj}){slice,layer}(region,:),'LineStyle',':','LineWidth',0.5,'color','g','Marker','v','MarkerSize', 3);
				set(gca,'XMinorTick','on','YMinorTick','on');
			end
		end
		axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
		text(0.5, 1,{strcat('slice#',num2str(slice),'-Peak Strain E:Control VS Isoproterenol',lower(fields{jj})),strcat('Distrib. of Circumf. Elemt:',', #Segments/slice:',', Transmural Interp:',', Z Interp:'),'Dark:Mean VS Light:±Standard Deviation'},'HorizontalAlignment','center','VerticalAlignment', 'top');
		hgexport(gcf,strcat('CONTROLvsISO_StrainCurves_DottedLines',fname,'_','_E',lower(fields{jj}),num2str(slice)),hgexport('factorystyle'), 'Format', 'jpeg');
		close gcf;
	end
end
save(strcat('Var_CONTROLvsISO_',fname),'aveBase','errBase','aveIso','errIso');

case 14
% Control VS Isoproterenol: Table and Bar Graph of Global Peak strains of 3 SA slices (Base vs Mid vs Apex)
clear strain_GlobalPeak_Base strain_GlobalPeak_Iso hStrain err;
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)				
	tmpBase = []; tmpMid = []; tmpApex = [];
	for ii = 1:nBaseline
		strain_GlobalPeak_Base.(fields{jj}){1,4+ii} = num2str(ii);
		[~,Tpeak] = max(abs(strain_BaseGlobal{ii}.(fields{jj})));
		tmpBase(ii) = strain_BaseGlobal{ii}.(fields{jj})(Tpeak);
		strain_GlobalPeak_Base.(fields{jj}){2,4+ii} = tmpBase(ii);
		[~,Tpeak] = max(abs(strain_SliceGlobal{ii}.(fields{jj})));
		tmpMid(ii) = strain_SliceGlobal{ii}.(fields{jj})(Tpeak);
		strain_GlobalPeak_Base.(fields{jj}){3,4+ii} = tmpMid(ii);
		[~,Tpeak] = max(abs(strain_ApexGlobal{ii}.(fields{jj})));
		tmpApex(ii) = strain_ApexGlobal{ii}.(fields{jj})(Tpeak);
		strain_GlobalPeak_Base.(fields{jj}){4,4+ii} = tmpApex(ii);
	end
	strain_GlobalPeak_Base.(fields{jj}){1,2} = 'Mean';
	strain_GlobalPeak_Base.(fields{jj}){1,3} = 'StdErr';
	strain_GlobalPeak_Base.(fields{jj}){1,4} = 'StdDev';
	strain_GlobalPeak_Base.(fields{jj}){2,1} = 'Base';
	strain_GlobalPeak_Base.(fields{jj}){3,1} = 'Mid';
	strain_GlobalPeak_Base.(fields{jj}){4,1} = 'Apex';
	hStrain.(fields{jj})(1,1) = mean(tmpBase);
	strain_GlobalPeak_Base.(fields{jj}){2,2} = hStrain.(fields{jj})(1,1);
	err.(fields{jj})(1,1) = sem(tmpBase);
	strain_GlobalPeak_Base.(fields{jj}){2,3} = err.(fields{jj})(1,1);
	strain_GlobalPeak_Base.(fields{jj}){2,4} = std(tmpBase);
	hStrain.(fields{jj})(2,1) = mean(tmpMid);
	err.(fields{jj})(2,1) = sem(tmpMid);
	strain_GlobalPeak_Base.(fields{jj}){3,2} = hStrain.(fields{jj})(2,1);
	strain_GlobalPeak_Base.(fields{jj}){3,3} = err.(fields{jj})(2,1);
	strain_GlobalPeak_Base.(fields{jj}){3,4} = std(tmpMid);
	hStrain.(fields{jj})(3,1) = mean(tmpApex);
	err.(fields{jj})(3,1) = sem(tmpApex);
	strain_GlobalPeak_Base.(fields{jj}){4,2} = hStrain.(fields{jj})(3,1);
	strain_GlobalPeak_Base.(fields{jj}){4,3} = err.(fields{jj})(3,1);
	strain_GlobalPeak_Base.(fields{jj}){4,4} = std(tmpApex);
	xlswrite(strcat('GlobalPeakStrain_Base_',fname),strain_GlobalPeak_Base.(fields{jj}), strcat('E_',lower(fields{jj})));
	
	tmpBase = []; tmpMid = []; tmpApex = [];
	for ii = 1:nFiles-nBaseline
		strain_GlobalPeak_Iso.(fields{jj}){1,4+ii} = num2str(ii);
		[~,Tpeak] = max(abs(strain_BaseGlobal{ii+nBaseline}.(fields{jj})));
		tmpBase(ii) = strain_BaseGlobal{ii+nBaseline}.(fields{jj})(Tpeak);
		strain_GlobalPeak_Iso.(fields{jj}){2,4+ii} = tmpBase(ii);
		[~,Tpeak] = max(abs(strain_SliceGlobal{ii+nBaseline}.(fields{jj})));
		tmpMid(ii) = strain_SliceGlobal{ii+nBaseline}.(fields{jj})(Tpeak);
		strain_GlobalPeak_Iso.(fields{jj}){3,4+ii} = tmpMid(ii);
		[~,Tpeak] = max(abs(strain_ApexGlobal{ii+nBaseline}.(fields{jj})));
		tmpApex(ii) = strain_ApexGlobal{ii+nBaseline}.(fields{jj})(Tpeak);
		strain_GlobalPeak_Iso.(fields{jj}){4,4+ii} = tmpApex(ii);
	end
	strain_GlobalPeak_Iso.(fields{jj}){1,2} = 'Mean';
	strain_GlobalPeak_Iso.(fields{jj}){1,3} = 'StdErr';
	strain_GlobalPeak_Iso.(fields{jj}){1,4} = 'StdDev';
	strain_GlobalPeak_Iso.(fields{jj}){2,1} = 'Base';
	strain_GlobalPeak_Iso.(fields{jj}){3,1} = 'Mid';
	strain_GlobalPeak_Iso.(fields{jj}){4,1} = 'Apex';
	hStrain.(fields{jj})(1,2) = mean(tmpBase);
	err.(fields{jj})(1,2) = sem(tmpBase);
	strain_GlobalPeak_Iso.(fields{jj}){2,2} = hStrain.(fields{jj})(1,2);
	strain_GlobalPeak_Iso.(fields{jj}){2,3} = err.(fields{jj})(1,2);
	strain_GlobalPeak_Iso.(fields{jj}){2,4} = std(tmpBase);
	hStrain.(fields{jj})(2,2) = mean(tmpMid);
	err.(fields{jj})(2,2) = sem(tmpMid);
	strain_GlobalPeak_Iso.(fields{jj}){3,2} = hStrain.(fields{jj})(2,2);
	strain_GlobalPeak_Iso.(fields{jj}){3,3} = err.(fields{jj})(2,2);
	strain_GlobalPeak_Iso.(fields{jj}){3,4} = std(tmpMid);
	hStrain.(fields{jj})(3,2) = mean(tmpApex);
	err.(fields{jj})(3,2) = sem(tmpApex);
	strain_GlobalPeak_Iso.(fields{jj}){4,2} = hStrain.(fields{jj})(3,2);
	strain_GlobalPeak_Iso.(fields{jj}){4,3} = err.(fields{jj})(3,2);
	strain_GlobalPeak_Iso.(fields{jj}){4,4} = std(tmpApex);
	xlswrite(strcat('GlobalPeakStrain_Iso_',fname),strain_GlobalPeak_Iso.(fields{jj}), strcat('E_',lower(fields{jj})));
	
	figure('units','normalized','outerposition',[0 0 0.5 0.5]);		
	h = bar(hStrain.(fields{jj}));
	% title(strcat(regions{region},transmural{layer},' Inter-Oberver LV Strains'));
	% xlabel('LV Segments');
	set(gca,'XTickLabel', {'Base','Mid-LV','Apex'});
	ylabel(['Peak Systolic',' E',lower(fields{jj}),' Strains']);
	set(gca,'XMinorTick','off','YMinorTick','on');
	grid on;
	legend({'CTL','ISO'},'Location','northeastoutside');
	%% In MATLAB R2014b and below:	
	for i=1:2
	  tmp = get(get(h(i),'children'),'xdata');
	  barsx(:,i)=mean(tmp,1);
	end
	% h2 = applyhatch_pluscolor(h1, '\-');
	% h2 = applyhatch_pluscolor(h1, '\-', err.(fields{jj}));
	% set(h2,'Visible','on');
	hold all;
	[~,Tpeak] = max(abs(tmpMid));
	errorbar(barsx,hStrain.(fields{jj}),err.(fields{jj}),'.');
	%{ 
	if  tmpMid(Tpeak)> 0
		errorbar(barsx,hStrain.(fields{jj}),zeros(3,2),err.(fields{jj}),'.');
	else
		errorbar(barsx,hStrain.(fields{jj}),err.(fields{jj}),zeros(3,2),'.');	
	end
	 %}
	%% In MATLAB R2016b:
	% h.CapSize = 12; 
	% errorbar(barsx,hStrain.(fields{jj}),zeros(3,2),err.(fields{jj}),'.','CapSize',10);
	%% Fail to get xdata:
	% errorbar_tick(h, 150);	
	hgexport(gcf,strcat('GlobalPeakStrain_CTLvsISO_',fname,'_E',lower(fields{jj})),hgexport('factorystyle'), 'Format', 'jpeg');
	hgsave(gcf,strcat('GlobalPeakStrain_CTLvsISO_',fname,'_E',lower(fields{jj})));
	close gcf;
end


case 15
% Control VS Isoproterenol: Segmentational strain curves with errorbar: 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 1]);
		for region = 1:nRegions	
			for layer = 1:nLayers-1
				subplot(nRegions,nLayers-1,(nLayers-1)*(region-1)+layer);
				tmp = []; tmp2 = [];
				for ii = 1:nBaseline
					tmp = [tmp; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
				end
				for ii = nBaseline+1:nFiles
					tmp2 = [tmp2; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
				end
				aveBase.(fields{jj}){slice,layer}(region,:) = mean(tmp);
				errBase.(fields{jj}){slice,layer}(region,:) = sem(tmp);
				aveIso.(fields{jj}){slice,layer}(region,:) = mean(tmp2);
				errIso.(fields{jj}){slice,layer}(region,:) = sem(tmp2);
				
				hbBase{layer} = errorbar(queryPoints,aveBase.(fields{jj}){slice,layer}(region,:),errBase.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','k','Marker','o','MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
				hold on;				
				hbIso{layer} = errorbar(queryPoints,aveIso.(fields{jj}){slice,layer}(region,:),errIso.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','b','Marker','^','MarkerSize', 8);
				% Adjust error bar width: only work in 2013b
				errorbar_tick(hbBase{layer}, 150);%140-(-1)^layer*20
				errorbar_tick(hbIso{layer}, 150);
				
				title(strcat(regions{region},transmural{layer},':p-value=',num2str(pval*100),'%'));
				xlabel('Cardiac Cycle(%)');
				ylabel(strcat('E',lower(fields{jj})));
				set(gca,'XMinorTick','on','YMinorTick','on');
				grid on;
				legend({' Baseline',' ISO'});
			end
		end
		tmp = strcat('CONTROLvsISO_StrainCurves_Segmentational_',fname,'_E',lower(fields{jj}),num2str(slice));
		hgexport(gcf,tmp,hgexport('factorystyle'), 'Format', 'jpeg');
		hgsave(gcf,tmp);
		close gcf;
	end
end
save(strcat('Var_CONTROLvsISO_errorbar_',fname),'aveBase','errBase','aveIso','errIso');

case 15.1
% Control VS Isoproterenol: Regional strain curves with errorbar: 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
layer = nLayers;
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 0.5]);
		for region = 1:nRegions	
			subplot(1,nRegions,region);
			tmp = []; tmp2 = [];
			for ii = 1:nBaseline
				tmp = [tmp; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
			end
			for ii = nBaseline+1:nFiles
				tmp2 = [tmp2; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
			end

			aveBase.(fields{jj}){slice,layer}(region,:) = mean(tmp);
			errBase.(fields{jj}){slice,layer}(region,:) = sem(tmp);
			aveIso.(fields{jj}){slice,layer}(region,:) = mean(tmp2);
			errIso.(fields{jj}){slice,layer}(region,:) = sem(tmp2);
			
			hbBase{layer} = errorbar(queryPoints,aveBase.(fields{jj}){slice,layer}(region,:),errBase.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','k','Marker','o','MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
			hold on;
			hbIso{layer} = errorbar(queryPoints,aveIso.(fields{jj}){slice,layer}(region,:),errIso.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','b','Marker','^','MarkerSize', 8);
			
			% Adjust error bar width: only work in 2013b
			errorbar_tick(hbBase{layer}, 150);%140-(-1)^layer*20
			errorbar_tick(hbIso{layer}, 150);
			
			title(strcat(regions{region},' segments'));
			xlabel('Cardiac Cycle(%)');
			ylabel(strcat('E',lower(fields{jj})));
			set(gca,'XMinorTick','on','YMinorTick','on');
			grid on;
			legend({' Baseline',' ISO'});

		end
		tmp = strcat('CONTROLvsISO_StrainCurves_Regional_',fname,'_E',lower(fields{jj}),num2str(slice));
		hgexport(gcf,tmp,hgexport('factorystyle'), 'Format', 'jpeg');
		hgsave(gcf,tmp);
		close gcf;
	end
end
save(strcat('Var_CONTROLvsISO_errorbar_Regional_',fname),'aveBase','errBase','aveIso','errIso');

case 15.2
% Control VS Isoproterenol: Transmural strain curves with errorbar: 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 0.5]);
		for layer = 1:nLayers-1
			subplot(1,nLayers-1,layer);
			tmp = []; tmp2 = [];
			for ii = 1:nBaseline
				tmp = [tmp; mean(normalData{ii}.(fields{jj}){slice,layer})];
			end
			for ii = nBaseline+1:nFiles
				tmp2 = [tmp2; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
			end

			aveBase.(fields{jj}){slice,layer} = mean(tmp);
			errBase.(fields{jj}){slice,layer} = sem(tmp);
			aveIso.(fields{jj}){slice,layer} = mean(tmp2);
			errIso.(fields{jj}){slice,layer} = sem(tmp2);

			hbBase{layer} = errorbar(queryPoints,aveBase.(fields{jj}){slice,layer},errBase.(fields{jj}){slice,layer},'LineWidth',1,'color','k','Marker','o','MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
			hold on;
			hbIso{layer} = errorbar(queryPoints,aveIso.(fields{jj}){slice,layer},errIso.(fields{jj}){slice,layer},'LineWidth',1,'color','b','Marker','^','MarkerSize', 8);
			
			% Adjust error bar width: only work in 2013b
			errorbar_tick(hbBase{layer}, 150);%140-(-1)^layer*20
			errorbar_tick(hbIso{layer}, 150);
			
			title(strcat(transmural{layer},' segments'));
			xlabel('Cardiac Cycle(%)');
			ylabel(strcat('E',lower(fields{jj})));
			set(gca,'XMinorTick','on','YMinorTick','on');
			grid on;
			legend({' Baseline',' ISO'});

		end
		tmp = strcat('CONTROLvsISO_StrainCurves_Transmural_',fname,'_E',lower(fields{jj}),num2str(slice));
		hgexport(gcf,tmp,hgexport('factorystyle'), 'Format', 'jpeg');
		hgsave(gcf,tmp);
		close gcf;
	end
end
save(strcat('Var_CONTROLvsISO_errorbar_Transmural_',fname),'aveBase','errBase','aveIso','errIso');

case 15.3
return
% Control VS Isoproterenol: Global strain curves with errorbar : Base vs Mid vs Apex: 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	figure('units','normalized','outerposition',[0 0 1 0.5]);

	subplot(1,3,1);
	tmp2 = mean(normalBaseGlobalData.(fields{jj})(1:nBaseline,:));
	tmp3 = sem(normalBaseGlobalData.(fields{jj})(1:nBaseline,:));
	hb = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{1},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
	% Adjust error bar width: only work in 2013b
	errorbar_tick(hb, 150);
	hold on;
	tmp2 = mean(normalBaseGlobalData.(fields{jj})(nBaseline+1:nFiles,:));
	tmp3 = sem(normalBaseGlobalData.(fields{jj})(nBaseline+1:nFiles,:));
	hb = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{4},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
	% Adjust error bar width: only work in 2013b
	errorbar_tick(hb, 150);
	title('Base');
	xlabel('Cardiac Cycle(%)');
	ylabel(strcat('E',lower(fields{jj})));
	set(gca,'XMinorTick','on','YMinorTick','on');
	grid on;
	legend({' CTL',' ISO'});
	
	subplot(1,3,2);
	tmp2 = mean(normalSliceGlobalData.(fields{jj})(1:nBaseline,:));
	tmp3 = sem(normalSliceGlobalData.(fields{jj})(1:nBaseline,:));
	hb = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{1},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
	% Adjust error bar width: only work in 2013b
	errorbar_tick(hb, 150);
	hold on;
	tmp2 = mean(normalSliceGlobalData.(fields{jj})(nBaseline+1:nFiles,:));
	tmp3 = sem(normalSliceGlobalData.(fields{jj})(nBaseline+1:nFiles,:));
	hb = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{4},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
	% Adjust error bar width: only work in 2013b
	errorbar_tick(hb, 150);
	title('Mid');
	xlabel('Cardiac Cycle(%)');
	ylabel(strcat('E',lower(fields{jj})));
	set(gca,'XMinorTick','on','YMinorTick','on');
	grid on;
	legend({' CTL',' ISO'});
	
	subplot(1,3,3);
	tmp2 = mean(normalApexGlobalData.(fields{jj})(1:nBaseline,:));
	tmp3 = sem(normalApexGlobalData.(fields{jj})(1:nBaseline,:));
	hb = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{1},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
	% Adjust error bar width: only work in 2013b
	errorbar_tick(hb, 150);
	hold on;
	tmp2 = mean(normalApexGlobalData.(fields{jj})(nBaseline+1:nFiles,:));
	tmp3 = sem(normalApexGlobalData.(fields{jj})(nBaseline+1:nFiles,:));
	hb = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{4},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
	% Adjust error bar width: only work in 2013b
	errorbar_tick(hb, 150);
	title('Apex');
	xlabel('Cardiac Cycle(%)');
	ylabel(strcat('E',lower(fields{jj})));
	set(gca,'XMinorTick','on','YMinorTick','on');
	grid on;
	legend({' CTL',' ISO'});
	
	tmp = strcat('CONTROLvsISO_StrainCurves_Global_',fname,'_E',lower(fields{jj}));
	hgexport(gcf,tmp,hgexport('factorystyle'), 'Format', 'jpeg');
	hgsave(gcf,tmp);
	close gcf;
end

case 15.4
% ALL-IN-ONE!!!
% Control VS Isoproterenol: Segmentational strain curves with errorbar: 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		for region = 1:nRegions	
			figure('units','normalized','outerposition',[0 0 1 1]);
			for layer = 1:nLayers-1
				tmp = []; tmp2 = [];
				for ii = 1:nBaseline
					tmp = [tmp; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
				end
				for ii = nBaseline+1:nFiles
					tmp2 = [tmp2; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
				end
				aveBase.(fields{jj}){slice,layer}(region,:) = mean(tmp);
				errBase.(fields{jj}){slice,layer}(region,:) = sem(tmp);
				aveIso.(fields{jj}){slice,layer}(region,:) = mean(tmp2);
				errIso.(fields{jj}){slice,layer}(region,:) = sem(tmp2);

				hbBase{layer} = errorbar(queryPoints,aveBase.(fields{jj}){slice,layer}(region,:),errBase.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','k','Marker',Marker{layer+3},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
				hold on;				
				hbIso{layer} = errorbar(queryPoints,aveIso.(fields{jj}){slice,layer}(region,:),errIso.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','b','Marker',Marker{layer+3},'MarkerSize', 8);
				% Adjust error bar width: only work in 2013b
				errorbar_tick(hbBase{layer}, 150);%140-(-1)^layer*20
				errorbar_tick(hbIso{layer}, 150);
			end
			% title(strcat(regions{region},transmural{layer}));
			xlabel('Cardiac Cycle(%)');
			ylabel(strcat('E',lower(fields{jj})));
			set(gca,'XMinorTick','on','YMinorTick','on');
			grid on;
			legend(transmural{1:nLayers-1});
			tmp = strcat('CONTROLvsISO_StrainCurves_Segmentational_',fname,'_E',lower(fields{jj}),num2str(slice),'_',regions(region));
			hgexport(gcf,tmp{:},hgexport('factorystyle'), 'Format', 'jpeg');
			hgsave(gcf,tmp{:});
			close gcf;
		end
	end
end
save(strcat('Var_CONTROLvsISO_errorbar_',fname),'aveBase','errBase','aveIso','errIso');

case 15.5
% ALL-IN-ONE!!!
% Control VS Isoproterenol: Regional strain curves with errorbar: 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
layer = nLayers;
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 1]);
		for region = 1:nRegions	
			tmp = []; tmp2 = [];
			for ii = 1:nBaseline
				tmp = [tmp; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
			end
			for ii = nBaseline+1:nFiles
				tmp2 = [tmp2; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
			end

			aveBase.(fields{jj}){slice,layer}(region,:) = mean(tmp);
			errBase.(fields{jj}){slice,layer}(region,:) = sem(tmp);
			aveIso.(fields{jj}){slice,layer}(region,:) = mean(tmp2);
			errIso.(fields{jj}){slice,layer}(region,:) = sem(tmp2);
			
			hbBase{layer} = errorbar(queryPoints,aveBase.(fields{jj}){slice,layer}(region,:),errBase.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','k','Marker',Marker{layer+3},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
			hold on;
			hbIso{layer} = errorbar(queryPoints,aveIso.(fields{jj}){slice,layer}(region,:),errIso.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','b','Marker',Marker{layer+3},'MarkerSize', 8);
			
			% Adjust error bar width: only work in 2013b
			errorbar_tick(hbBase{layer}, 150);%140-(-1)^layer*20
			errorbar_tick(hbIso{layer}, 150);
		end
		% title(strcat(regions{region},transmural{layer}));
		xlabel('Cardiac Cycle(%)');
		ylabel(strcat('E',lower(fields{jj})));
		set(gca,'XMinorTick','on','YMinorTick','on');
		grid on;
		legend(regions);
		tmp = strcat('CONTROLvsISO_StrainCurves_Regional_',fname,'_E',lower(fields{jj}),num2str(slice));
		hgexport(gcf,tmp,hgexport('factorystyle'), 'Format', 'jpeg');
		hgsave(gcf,tmp);
		close gcf;
	end
end
save(strcat('Var_CONTROLvsISO_errorbar_Regional_',fname),'aveBase','errBase','aveIso','errIso');

case 15.6
% ALL-IN-ONE!!!
% Control VS Isoproterenol: Transmural strain curves with errorbar: 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		figure('units','normalized','outerposition',[0 0 1 1]);
		for layer = 1:nLayers-1
			tmp = []; tmp2 = [];
			for ii = 1:nBaseline
				tmp = [tmp; mean(normalData{ii}.(fields{jj}){slice,layer})];
			end
			for ii = nBaseline+1:nFiles
				tmp2 = [tmp2; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
			end

			aveBase.(fields{jj}){slice,layer} = mean(tmp);
			errBase.(fields{jj}){slice,layer} = sem(tmp);
			aveIso.(fields{jj}){slice,layer} = mean(tmp2);
			errIso.(fields{jj}){slice,layer} = sem(tmp2);

			hbBase{layer} = errorbar(queryPoints,aveBase.(fields{jj}){slice,layer},errBase.(fields{jj}){slice,layer},'LineWidth',1,'color','k','Marker',Marker{layer+3},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
			hold on;
			hbIso{layer} = errorbar(queryPoints,aveIso.(fields{jj}){slice,layer},errIso.(fields{jj}){slice,layer},'LineWidth',1,'color','b','Marker',Marker{layer+3},'MarkerSize', 8);
			
			% Adjust error bar width: only work in 2013b
			errorbar_tick(hbBase{layer}, 150);%140-(-1)^layer*20
			errorbar_tick(hbIso{layer}, 150);	
		end
		% title(strcat(regions{region},transmural{layer}));
		xlabel('Cardiac Cycle(%)');
		ylabel(strcat('E',lower(fields{jj})));
		set(gca,'XMinorTick','on','YMinorTick','on');
		grid on;
		legend(transmural{1:nLayers-1});
		tmp = strcat('CONTROLvsISO_StrainCurves_Transmural_',fname,'_E',lower(fields{jj}),num2str(slice));
		hgexport(gcf,tmp{:},hgexport('factorystyle'), 'Format', 'jpeg');
		hgsave(gcf,tmp{:});
		close gcf;
	end
end
save(strcat('Var_CONTROLvsISO_errorbar_Transmural_',fname),'aveBase','errBase','aveIso','errIso');

case 15.7
% Plotted Separately!!!
% Control VS Isoproterenol: Regional strain curves with errorbar: 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
layer = nLayers;
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		for region = 1:nRegions	
			figure('units','normalized','outerposition',[0 0 1 1]);
			tmp = []; tmp2 = [];
			for ii = 1:nBaseline
				tmp = [tmp; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
			end
			for ii = nBaseline+1:nFiles
				tmp2 = [tmp2; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
			end

			aveBase.(fields{jj}){slice,layer}(region,:) = mean(tmp);
			errBase.(fields{jj}){slice,layer}(region,:) = sem(tmp);
			aveIso.(fields{jj}){slice,layer}(region,:) = mean(tmp2);
			errIso.(fields{jj}){slice,layer}(region,:) = sem(tmp2);
			
			hbBase{layer} = errorbar(queryPoints,aveBase.(fields{jj}){slice,layer}(region,:),errBase.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','k','Marker','o','MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
			hold on;
			hbIso{layer} = errorbar(queryPoints,aveIso.(fields{jj}){slice,layer}(region,:),errIso.(fields{jj}){slice,layer}(region,:),'LineWidth',1,'color','b','Marker','^','MarkerSize', 8);
			
			% Adjust error bar width: only work in 2013b
			errorbar_tick(hbBase{layer}, 150);%140-(-1)^layer*20
			errorbar_tick(hbIso{layer}, 150);
			
			% title(strcat(regions{region},transmural{layer}));
			xlabel('Cardiac Cycle(%)');
			ylabel(strcat('E',lower(fields{jj})));
			set(gca,'XMinorTick','on','YMinorTick','on');
			grid on;
			legend({' Baseline',' ISO'});
			tmp = strcat('CONTROLvsISO_StrainCurves_Regional_',fname,'_E',lower(fields{jj}),num2str(slice),'_',regions(region));
			hgexport(gcf,tmp,hgexport('factorystyle'), 'Format', 'jpeg');
			hgsave(gcf,tmp);
			close gcf;
		end
	end
end
save(strcat('Var_CONTROLvsISO_errorbar_Regional_',fname),'aveBase','errBase','aveIso','errIso');

case 15.8
% Plotted Separately!!!
% Control VS Isoproterenol: Transmural strain curves with errorbar: 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		for layer = 1:nLayers-1
			figure('units','normalized','outerposition',[0 0 1 1]);
			tmp = []; tmp2 = [];
			for ii = 1:nBaseline
				tmp = [tmp; mean(normalData{ii}.(fields{jj}){slice,layer})];
			end
			for ii = nBaseline+1:nFiles
				tmp2 = [tmp2; normalData{ii}.(fields{jj}){slice,layer}(region,:)];
			end

			aveBase.(fields{jj}){slice,layer} = mean(tmp);
			errBase.(fields{jj}){slice,layer} = sem(tmp);
			aveIso.(fields{jj}){slice,layer} = mean(tmp2);
			errIso.(fields{jj}){slice,layer} = sem(tmp2);

			hbBase{layer} = errorbar(queryPoints,aveBase.(fields{jj}){slice,layer},errBase.(fields{jj}){slice,layer},'LineWidth',1,'color','k','Marker','o','MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
			hold on;
			hbIso{layer} = errorbar(queryPoints,aveIso.(fields{jj}){slice,layer},errIso.(fields{jj}){slice,layer},'LineWidth',1,'color','b','Marker','^','MarkerSize', 8);
			
			% Adjust error bar width: only work in 2013b
			errorbar_tick(hbBase{layer}, 150);%140-(-1)^layer*20
			errorbar_tick(hbIso{layer}, 150);
			
			% title(strcat(regions{region},transmural{layer}));
			xlabel('Cardiac Cycle(%)');
			ylabel(strcat('E',lower(fields{jj})));
			set(gca,'XMinorTick','on','YMinorTick','on');
			grid on;
			legend({' Baseline',' ISO'});
			tmp = strcat('CONTROLvsISO_StrainCurves_Transmural_',fname,'_E',lower(fields{jj}),num2str(slice),'_',transmural(layer));
			hgexport(gcf,tmp{:},hgexport('factorystyle'), 'Format', 'jpeg');
			hgsave(gcf,tmp{:});
			close gcf;
		end
	end
end
save(strcat('Var_CONTROLvsISO_errorbar_Transmural_',fname),'aveBase','errBase','aveIso','errIso');

case 15.9
% Plotted Separately!!!
% Control VS Isoproterenol: Global strain curves with errorbar : Base vs Mid vs Apex: 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	figure('units','normalized','outerposition',[0 0 1 1]);

	tmp2 = mean(normalBaseGlobalData.(fields{jj}));
	tmp3 = sem(normalBaseGlobalData.(fields{jj}));
	hb = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{4},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
	% Adjust error bar width: only work in 2013b
	errorbar_tick(hb, 150);
	hold on;

	tmp2 = mean(normalSliceGlobalData.(fields{jj}));
	tmp3 = sem(normalSliceGlobalData.(fields{jj}));
	hb = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{1},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
	% Adjust error bar width: only work in 2013b
	errorbar_tick(hb, 150);
	
	tmp2 = mean(normalApexGlobalData.(fields{jj}));
	tmp3 = sem(normalApexGlobalData.(fields{jj}));
	hb = errorbar(queryPoints,tmp2,tmp3,'LineWidth',1,'color','k','Marker',Marker{5},'MarkerSize', 8);% 'MarkerFaceColor', [.3 1 .3],
	% Adjust error bar width: only work in 2013b
	errorbar_tick(hb, 150);
	
	% title(strcat(regions{region},transmural{layer}));
	xlabel('Cardiac Cycle(%)');
	ylabel(strcat('E',lower(fields{jj})));
	set(gca,'XMinorTick','on','YMinorTick','on');
	grid on;
	legend({' Base',' Mid',' Apex'});
	tmp = strcat('CONTROLvsISO_StrainCurves_Global_',fname,'_E',lower(fields{jj}));
	hgexport(gcf,tmp,hgexport('factorystyle'), 'Format', 'jpeg');
	hgsave(gcf,tmp);
	close gcf;
end

case 16
%% Control VS Isoproterenol: two-sample t-test with unpooled variance & wilcoxon rank sum test
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		% Control VS Isoproterenol: Two-Sample t-Test for Segmentational strains: 
		for layer = 1:nLayers
			% layer = 4 for Transmural strains: 		
			for region = 1:nRegions	
				for ii = 1:nFiles
					% Tpeak is for sign preserve!!!
					[~,Tpeak] = max(abs(data{ii}.(fields{jj}){slice,layer}(region,:)));
					Epeak.(fields{jj}){slice,layer}(region,ii) = data{ii}.(fields{jj}){slice,layer}(region,Tpeak);
				end
				
				tmp = Epeak.(fields{jj}){slice,layer}(region,1:nBaseline);
				tmp2 = Epeak.(fields{jj}){slice,layer}(region,nBaseline+1:nFiles);

				% two-sample t-test with unpooled variance
				[~,Pval.(fields{jj}){slice}(region,layer),ci.(fields{jj}){slice,layer}(region,:),stats.(fields{jj}){slice,layer}{region}] = ttest2(tmp,tmp2,'Vartype','unequal');
				tPval(nLayers*(jj-1)+layer,region) = Pval.(fields{jj}){slice}(region,layer);
				% wilcoxon rank sum test
				[PVAL.(fields{jj}){slice}(region,layer),~,STATS.(fields{jj}){slice,layer}{region}] = ranksum(tmp,tmp2);
				wPval(nLayers*(jj-1)+layer,region) = PVAL.(fields{jj}){slice}(region,layer);
			end
			
			%% Control VS Isoproterenol: Two-Sample t-Test for Regional strains:				
			for ii = 1:nFiles
				tmp3 = mean(data{ii}.(fields{jj}){slice,layer}(:,:));
				[~,Tpeak] = max(abs(tmp3));
				Epeak.(fields{jj}){slice,layer}(nRegions+1,ii) = tmp3(Tpeak);
			end
			% The fifth section of Pval:Endo -> Mid -> Epi -> Global
			tmp = Epeak.(fields{jj}){slice,layer}(nRegions+1,1:nBaseline);
			tmp2 = Epeak.(fields{jj}){slice,layer}(nRegions+1,nBaseline+1:nFiles);
			[~,Pval.(fields{jj}){slice}(nRegions+1,layer),ci.(fields{jj}){slice,layer}(nRegions+1,:),stats.(fields{jj}){slice,layer}{nRegions+1}] = ttest2(tmp,tmp2,'Vartype','unequal');
			tPval(nLayers*(jj-1)+layer,nRegions+1) = Pval.(fields{jj}){slice}(nRegions+1,layer);
			[PVAL.(fields{jj}){slice}(nRegions+1,layer),~,STATS.(fields{jj}){slice,layer}{nRegions+1}] = ranksum(tmp,tmp2);
			wPval(nLayers*(jj-1)+layer,nRegions+1) = PVAL.(fields{jj}){slice}(nRegions+1,layer);			
		end
		
		%{ 
		figure('units','normalized','outerposition',[0 0 1 1]);		
		bar(Pval.(fields{jj}){slice});
		title(strcat(regions{region},transmural{layer},' Inter-Oberver LV Strains'));
		xlabel('LV Segments');
		set(gca,'XTickLabel', {regions{:},'Global'});
		ylabel('Coefficient of Variation(%)');
		set(gca,'XMinorTick','off','YMinorTick','on');
		grid on;
		legend(transmural);
		applyhatch_pluscolor(gcf, '\/x-');
		% set(gcf,'Visible','on');
		hgexport(gcf,strcat('Pval_',fname,'_E',lower(fields{jj}),num2str(slice)),hgexport('factorystyle'), 'Format', 'jpeg');
		hgsave(gcf,strcat('Pval_',fname,'_E',lower(fields{jj}),num2str(slice)));
		close gcf;
		close gcf;
		table = [{'Pval','Endo','Mid','Epi','Myocardium'};[{'Septal','Inferior','Lateral','Anterior','Global'}',num2cell(Pval.(fields{jj}){slice})]];
		xlswrite(strcat('Summary_Pval_',fname),table, strcat('E_',lower(fields{jj})));
		%}
	end
end
table = [{'two-sample t-test with unpooled variance','Pval','Septal','Inferior','Lateral','Anterior','Global'};[reshape([fields;repmat({''},1,numel(fields));repmat(repmat({''},1,numel(fields)),2,1)],[],1),repmat({'Endo','Mid','Epi','Myocardium'}',numel(fields),1),num2cell(tPval)]];
xlswrite(strcat('Summary_Pval_t-test_',fname),table, 'Langragian Strains');
table = [{'wilcoxon rank sum test','Pval','Septal','Inferior','Lateral','Anterior','Global'};[reshape([fields;repmat({''},1,numel(fields));repmat(repmat({''},1,numel(fields)),2,1)],[],1),repmat({'Endo','Mid','Epi','Myocardium'}',numel(fields),1),num2cell(wPval)]];
xlswrite(strcat('Summary_Pval_wilcoxon_',fname),table,'Langragian Strains');
save(strcat('Var_Pval_',fname),'Pval','tPval','PVAL','wPval');

case 16.1
%% CoV for global strains: Base VS Mid VS Apex 
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
clear hCoV;
for jj = 1:numel(fields)			
	hCoV.(fields{jj}){1,1} = 'Base';
	hCoV.(fields{jj}){1,2} = 'Mid';
	hCoV.(fields{jj}){1,3} = 'Apex';
	tmpBase3 = []; tmpMid3 = []; tmpApex3 = [];
	for ii = 1:nFiles
		[~,Tpeak] = max(abs(strain_BaseGlobal{ii}.(fields{jj})));
		tmpBase3{ii} = strain_BaseGlobal{ii}.(fields{jj})(Tpeak);
		[~,Tpeak] = max(abs(strain_SliceGlobal{ii}.(fields{jj})));
		tmpMid3{ii} = strain_SliceGlobal{ii}.(fields{jj})(Tpeak);
		[~,Tpeak] = max(abs(strain_ApexGlobal{ii}.(fields{jj})));
		tmpApex3{ii} = strain_ApexGlobal{ii}.(fields{jj})(Tpeak);
	end
	tmpBase = []; tmpBase1 = 0; tmpBase2 = 0; 
	tmpMid = []; tmpMid1 = 0; tmpMid2 = 0; 
	tmpApex = []; tmpApex1 = 0; tmpApex2 = 0; 
	for ii = 1:2:nFiles
		tmpBase = [tmpBase3{ii},tmpBase3{ii+1}];
		tmpBase1 = tmpBase1 + mean(tmpBase);
		tmpBase2 = tmpBase2 + std(tmpBase);
		tmpMid = [tmpMid3{ii},tmpMid3{ii+1}];
		tmpMid1 = tmpMid1 + mean(tmpMid);
		tmpMid2 = tmpMid2 + std(tmpMid);
		tmpApex = [tmpApex3{ii},tmpApex3{ii+1}];
		tmpApex1 = tmpApex1 + mean(tmpApex);
		tmpApex2 = tmpApex2 + std(tmpApex);
	end
	hCoV.(fields{jj}){2,1} = tmpBase2/abs(tmpBase1)*100;
	hCoV.(fields{jj}){2,2} = tmpMid2/abs(tmpMid1)*100;
	hCoV.(fields{jj}){2,3} = tmpApex2/abs(tmpApex1)*100;		
	xlswrite(strcat('CoV_baseVSmidVSapex_',fname),hCoV.(fields{jj}), strcat('E_',lower(fields{jj})));
end

case 17
% Output raw data to Excel for statistical test + Std. Dev. & Summary of Peak Strains:
fields = {'RR', 'CC', 'LL', 'RC', 'RL', 'CL', 'CLShearAngle'};
% fields = {'CC', 'LL', 'CL', 'CLShearAngle'};
clear Epeak;
for jj = 1:numel(fields)					
	for slice = 1:nSlices
		for layer = 1:nLayers
			for region = 1:nRegions	
				for ii = 1:nBaseline
					[~,Tpeak] = max(abs(normalData{ii}.(fields{jj}){slice,layer}(region,:)));
					Epeak.(fields{jj}){slice,layer}(region,ii) = data{ii}.(fields{jj}){slice,layer}(region,Tpeak);
				end
				xlswrite(strcat('Statistics_Raw_Control_',fields{jj},'.xls'),Epeak.(fields{jj}){slice,layer}(region,1:nBaseline), strcat(regions{region}, transmural{layer}));
				
				tMean(nLayers*(jj-1)+layer,region) = mean(Epeak.(fields{jj}){slice,layer}(region,1:nBaseline));
				tSEM(nLayers*(jj-1)+layer,region) = sem(Epeak.(fields{jj}){slice,layer}(region,1:nBaseline));
				
				tEpeak{nLayers*(jj-1)+layer,region} = sprintf(['%0.2f',177,'%0.2f'], [tMean(nLayers*(jj-1)+layer,region);tSEM(nLayers*(jj-1)+layer,region)]);
			end
			
			for ii = 1:nBaseline
				tmp3 = mean(data{ii}.(fields{jj}){slice,layer}(:,:));
				[~,Tpeak] = max(abs(tmp3));
				Epeak.(fields{jj}){slice,layer}(nRegions+1,ii) = tmp3(Tpeak);					
			end			
			xlswrite(strcat('Statistics_Raw_Control_',fields{jj},'.xls'),Epeak.(fields{jj}){slice,layer}(nRegions+1,1:nBaseline), strcat('Global', transmural{layer}));
			
			tMean(nLayers*(jj-1)+layer,nRegions+1) = mean(Epeak.(fields{jj}){slice,layer}(nRegions+1,1:nBaseline));
			tSEM(nLayers*(jj-1)+layer,nRegions+1) = sem(Epeak.(fields{jj}){slice,layer}(nRegions+1,1:nBaseline));

			tEpeak{nLayers*(jj-1)+layer,nRegions+1} = sprintf(['%0.2f',177,'%0.2f'], [tMean(nLayers*(jj-1)+layer,nRegions+1);tSEM(nLayers*(jj-1)+layer,nRegions+1)]);
		end
	end
end
table = [{'Components','PeakStrain','Septal','Inferior','Lateral','Anterior','Global'};[reshape([fields;repmat({''},1,numel(fields));repmat(repmat({''},1,numel(fields)),2,1)],[],1),repmat({'Endo','Mid','Epi','Myocardium'}',numel(fields),1),tEpeak]];
xlswrite(strcat('Summary_PeakStrain_Control_',fname),table);

table = [{'Components','PeakStrain','Septal','Inferior','Lateral','Anterior','Global'};[reshape([fields;repmat({''},1,numel(fields));repmat(repmat({''},1,numel(fields)),2,1)],[],1),repmat({'Endo','Mid','Epi','Myocardium'}',numel(fields),1),num2cell(tMean)]];
xlswrite(strcat('Summary_PeakStrain_Control__Mean',fname),table);

table = [{'Components','PeakStrain','Septal','Inferior','Lateral','Anterior','Global'};[reshape([fields;repmat({''},1,numel(fields));repmat(repmat({''},1,numel(fields)),2,1)],[],1),repmat({'Endo','Mid','Epi','Myocardium'}',numel(fields),1),num2cell(tSEM)]];
xlswrite(strcat('Summary_PeakStrain_Control__SEM',fname),table);

save(strcat('Var_PeakStrain_Control_',fname),'tEpeak','tMean','tSEM');

for jj = 1:numel(fields)					
	for slice = 1:nSlices
		for layer = 1:nLayers
			for region = 1:nRegions	
				tmp = [];
				for ii = nBaseline+1:nFiles
					[~,Tpeak] = max(abs(data{ii}.(fields{jj}){slice,layer}(region,:)));
					Epeak.(fields{jj}){slice,layer}(region,ii) = data{ii}.(fields{jj}){slice,layer}(region,Tpeak);
				end
				xlswrite(strcat('Statistics_Raw_Iso_',fields{jj},'.xls'),Epeak.(fields{jj}){slice,layer}(region,nBaseline+1:nFiles), strcat(regions{region}, transmural{layer}));
				
				tMean(nLayers*(jj-1)+layer,region) = mean(Epeak.(fields{jj}){slice,layer}(region,nBaseline+1:nFiles));
				tSEM(nLayers*(jj-1)+layer,region) = sem(Epeak.(fields{jj}){slice,layer}(region,nBaseline+1:nFiles));
				
				tEpeak{nLayers*(jj-1)+layer,region} = sprintf(['%0.2f',177,'%0.2f'], [tMean(nLayers*(jj-1)+layer,region);tSEM(nLayers*(jj-1)+layer,region)]);
			end
			
			for ii = nBaseline+1:nFiles
				tmp3 = mean(data{ii}.(fields{jj}){slice,layer}(:,:));
				[~,Tpeak] = max(abs(tmp3));
				Epeak.(fields{jj}){slice,layer}(nRegions+1,ii) = tmp3(Tpeak);					
			end			
			xlswrite(strcat('Statistics_Raw_Iso_',fields{jj},'.xls'),Epeak.(fields{jj}){slice,layer}(nRegions+1,nBaseline+1:nFiles), strcat('Global', transmural{layer}));
			
			tMean(nLayers*(jj-1)+layer,nRegions+1) = mean(Epeak.(fields{jj}){slice,layer}(nRegions+1,nBaseline+1:nFiles));
			tSEM(nLayers*(jj-1)+layer,nRegions+1) = sem(Epeak.(fields{jj}){slice,layer}(nRegions+1,nBaseline+1:nFiles));

			tEpeak{nLayers*(jj-1)+layer,nRegions+1} = sprintf(['%0.2f',177,'%0.2f'], [tMean(nLayers*(jj-1)+layer,nRegions+1);tSEM(nLayers*(jj-1)+layer,nRegions+1)]);
		end
	end
end
table = [{'Components','PeakStrain','Septal','Inferior','Lateral','Anterior','Global'};[reshape([fields;repmat({''},1,numel(fields));repmat(repmat({''},1,numel(fields)),2,1)],[],1),repmat({'Endo','Mid','Epi','Myocardium'}',numel(fields),1),tEpeak]];
xlswrite(strcat('Summary_PeakStrain_Iso_',fname),table);

table = [{'Components','PeakStrain','Septal','Inferior','Lateral','Anterior','Global'};[reshape([fields;repmat({''},1,numel(fields));repmat(repmat({''},1,numel(fields)),2,1)],[],1),repmat({'Endo','Mid','Epi','Myocardium'}',numel(fields),1),num2cell(tMean)]];
xlswrite(strcat('Summary_PeakStrain_Iso__Mean',fname),table);

table = [{'Components','PeakStrain','Septal','Inferior','Lateral','Anterior','Global'};[reshape([fields;repmat({''},1,numel(fields));repmat(repmat({''},1,numel(fields)),2,1)],[],1),repmat({'Endo','Mid','Epi','Myocardium'}',numel(fields),1),num2cell(tSEM)]];
xlswrite(strcat('Summary_PeakStrain_Iso__SEM',fname),table);

save(strcat('Var_PeakStrain_Iso_',fname),'Epeak','tEpeak','tMean','tSEM');

end


end