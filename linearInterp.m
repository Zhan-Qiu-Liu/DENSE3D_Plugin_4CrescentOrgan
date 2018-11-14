function linearInterp(dataObj,roiObj)
% Copyright (c) 2016 DENSEanalysis Contributors
% The codes below are improved from the function "linearInterp" in "roitool.m".
% It's much more versatile to all scenarios
% Modified By: Zhan-Qiu Liu (lafeir.lew@gmail.com)
% Last Modified: 22:06 October 28, 2018

	% import plugins.DENSE3D_Plugin_4CrescentOrgan.*
	% handles = guidata(self.hfig(1));
	% dataObj = handles.hdata;
	% roiObj = handles.hdense.hroi;

 	if isempty(roiObj.ROIIndex); return; end
	
	hasNodes = ~cellfun(@isempty, dataObj.roi(roiObj.ROIIndex).Position);
	[nFrames, nLayers] = size(dataObj.roi(roiObj.ROIIndex).Position);
	
   % test if prev not empty:
    % if roiObj.ROIFrame==1 || isempty(dataObj.roi(roiObj.ROIIndex).Position{roiObj.ROIFrame-1,1})
	if isempty(find(hasNodes(1:nFrames))); return; end
	
	f1 = find(hasNodes(1:roiObj.ROIFrame+(nLayers-1)*nFrames-1), 1, 'last');
	[f1, ~] = ind2sub([nFrames nLayers], f1);
	% f1 = roiObj.ROIFrame - 1;
	
	% search down column first then next row:
	f2 = find(hasNodes(roiObj.ROIFrame+1:end), 1, 'first');
	[f2, ~] = ind2sub([nFrames nLayers], roiObj.ROIFrame+f2);
	
	if f2 == f1
		% n=0 -> diff=NaN (inf.)
		nodes = dataObj.roi(roiObj.ROIIndex).Position(f1,:);
	else
		m = mod(roiObj.ROIFrame - f1, nFrames);
		n = mod(f2 - f1, nFrames);
		
		nodes1 = dataObj.roi(roiObj.ROIIndex).Position(f1,:);
		nodes2 = dataObj.roi(roiObj.ROIIndex).Position(f2,:);
		
		num1 = cell2mat(cellfun(@(x)size(x, 1), nodes1, 'uni', 0));
		num2 = cell2mat(cellfun(@(x)size(x, 1), nodes2, 'uni', 0));
		
		if ~isequal(num1, num2)
			errordlg('Contours must have the same number of control points', 'Invalid Input');
		end
		
		diff = cellfun(@minus, nodes2, nodes1, 'uni', 0);
		diff = cellfun(@(x)x*m/n, diff, 'uni', 0); 
		nodes = cellfun(@plus, nodes1, diff, 'uni', 0);
	end
	
	isCurved = dataObj.roi(roiObj.ROIIndex).IsCurved(f1,:);
	isCorner = dataObj.roi(roiObj.ROIIndex).IsCorner(f1,:);
	
	current = struct(...
		'Position', {nodes}, ...
		'IsClosed',{roiObj.cLine.IsClosed}, ...
		'IsCorner',{isCorner}, ...
		'IsCurved',{isCurved});
		
	dataObj.updateROI(roiObj.ROIIndex, roiObj.ROIFrame, current);
	
	redraw(roiObj);
end