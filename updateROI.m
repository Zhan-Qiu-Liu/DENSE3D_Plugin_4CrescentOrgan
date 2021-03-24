function imageryROIfcn(self,handles,flag)
% Copyright (c) 2016 DENSEanalysis Contributors
% The codes below are modified from the function "switchstate" in "DENSEanalysis.m".
% Last Modified: 11:10 September 13, 2019
% Modified By: Zhan-Qiu Liu (lafeir.lew@gmail.com)

	import plugins.dense3D_plugin_4crescentorgan.*
	%{ 
	handles = guidata(obj.Parent);
	% hroi = get(handles.hfig, 'UserData');
	obj = handles.hdata;
	if flag; keyboard; end
	
    if isempty(obj.seq) || isempty(obj.dns)
        return
    end
	%}
	
	didx = handles.hdense.DENSEIndex;
	ridx = handles.hdense.ROIIndex;
	frame = handles.hdense.Frame;

	imdata = handles.hdata.imagedata(didx);
	%% Contours from interpolating points
	cndata = handles.hdata.contourdata(ridx, frame);
	roi = handles.hdata.roi(ridx);

	% initialize any empty phase handles.hdata to zero
	tags = {'Xpha','Ypha'};
	for ti = 1:numel(tags)
		if isempty(imdata.(tags{ti}))
			imdata.(tags{ti}) = zeros(size(imdata.Mag));
			imdata.EncFreq(ti) = 1;
			imdata.Scale(ti)   = 1;
		end
	end

	% ensure DENSE/ROI indices overlap
	sidx = [imdata.MagIndex; imdata.PhaIndex];

	if isempty(intersect(cndata.SeqIndex(:),sidx(:)))
		roibuf = sprintf('%d,',   cndata.SeqIndex);
		dnsbuf = sprintf('%d/%d,',sidx(~isnan(sidx)));
		error(sprintf('%s:indicesDoNotMatch',mfilename),'%s',...
			'The ROI sequence indices [',roibuf(1:end-1),...
			'] do not match the DENSE sequence indices [',...
			dnsbuf(1:end-1),']');
	end

	% reduce contours
	cndata.Contour = cndata.Contour(frame,:);
	cndata.ContourFrame = frame;

	% if the current resident options "spl" is of the same indices,
	% gather those parameters & pass to the DENSESPLINE function
	options = struct;
	tags = {'ResampleMethod','SpatialSmoothing',...
		'TemporalOrder','Xseed','Yseed'};
	if ~isempty(self.mgsoptions)
		opts = self.mgsoptions;
		dlast = handles.hdata.UIDtoIndexDENSE(opts.DENSEUID);
		rlast = handles.hdata.UIDtoIndexROI(opts.ROIUID);
		if isequal(dlast,didx) && isequal(rlast,ridx)
			for ti = 1:numel(tags)
				options.(tags{ti}) = opts.(tags{ti});
			end

			if opts.FramesForAnalysis(1)<=frame && ...
			frame<=opts.FramesForAnalysis(2)
				options.FramesForAnalysis = opts.FramesForAnalysis;
			end
		end
	end

	if any(strcmpi(cndata.ROIType, {'open', 'closed'}))
		options.UnwrapConnectivity = 8;
	end

	options.Method = handles.hdense.Method;

	% initialize waitbar timer
	hwait = self.waitbar(...
		'String', 'Motion Guided Segmentation...', ...
		'WindowStyle', 'modal', ...
		'AllowClose', false);
	hwait.start();

	% Perform analysis
	pos = roi.Position(frame(1),:);
	try
		[newpos,options] = DENSEmgs(imdata,cndata,...
			'Position', pos,...
			'OptionsPanel', true, ...
			'Visualize', false, ...
			options);
	catch ERR
		hwait.stop;
		delete(hwait);
		str = ['Motion Guided Segmentation was unsuccessful. ',...
			'If you attempt MGS again, try selecting different ',...
			'unwrapped pixels, or increase the "smoothness" value ',...
			'to accomodate images with low SNR.  Alternatively, ',...
			'you may manually segment the dataset.'];
		h = warndlg(str,'MGS unsuccessful','modal');
		waitfor(h);
		disp('MGS report (unsuccessful):')
		disp(getReport(ERR))
		return
	end

	% stop/delete the waitbar
	hwait.stop; pause(1e-5);
	delete(hwait);

	% ensure newpos is the same size as pos
	if isempty(newpos)
		return
	elseif ~ismatrix(newpos) || ~all(size(roi.Position)==size(newpos))
		error(sprintf('%s:invalidMGS',mfilename),'%s',...
			'The Motion Guide Segmentation output position size ',...
			'does not match the original position size.');
	end

	% save handles.hdata
	Nfr = size(newpos,1);
	for fr = 1:Nfr
		if all(cellfun(@(p)~isempty(p),newpos(fr,:)))
			r.Position = newpos(fr,:);
			r.IsClosed = roi.IsClosed(frame,:);
			r.IsCurved = roi.IsCurved(frame,:);
			r.IsCorner = roi.IsCorner(frame,:);
			handles.hdata.updateROI(ridx, fr, r);
		end
	end

	self.mgsoptions = options;
	self.mgsoptions.DENSEUID = handles.hdata.dns(didx).UID;
	self.mgsoptions.ROIUID   = handles.hdata.roi(ridx).UID;
end			
