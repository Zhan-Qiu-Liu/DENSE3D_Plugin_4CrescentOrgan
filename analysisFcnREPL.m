function analysisFcnREPL(self,obj,didx,ridx,seedframe,varargin)
% Copyright (c) 2016 DENSEanalysis Contributors
% The codes below are modified from the function "analysisFcn" in "DENSEdata.m".
% Last Modified: 22:28 July 28, 2017
% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)

	import plugins.DENSE3D_Plugin_4CrescentOrgan.*
	
%% ANALYSIS FUNCTION
	% gather data
	imdata = imagedata(obj,didx);
	cndata = contourDataREPL(obj,ridx);

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

	% check seedframe
	rng = cndata.ValidFrames;
	if ~isnumeric(seedframe) || ~isscalar(seedframe) || ...
	   seedframe < rng(1) || rng(2) < seedframe
		error(sprintf('%s:invalidFrame',mfilename),...
			'Frame is invalid.');
	end

	% if the current resident analysis "spl" is of the same indices,
	% gather those parameters & pass to the DENSESPLINE function
	options = struct('FramesForAnalysis',cndata.ValidFrames);
	tagsA = {'ResampleMethod','SpatialSmoothing','TemporalOrder'};
	tagsB = {'Xseed','Yseed','Zseed'};
	if isfield(self.dns(didx),'seed') && seedframe<=size(self.dns(didx).seed,2) && ~isempty(self.dns(didx).seed{1,seedframe})
		for ti = 1:numel(tagsB)
			options.(tagsB{ti}) = self.dns(didx).seed{ti,seedframe};
		end	
	end
	
	if ~isempty(self.spl)
		for ti = 1:numel(tagsA)
			options.(tagsA{ti}) = self.spl.(tagsA{ti});
		end

		dlast = UIDtoIndexDENSE(obj,self.spl.DENSEUID);
		rlast = UIDtoIndexROI(obj,self.spl.ROIUID);
		if isequal(dlast,didx) && isequal(rlast,ridx)
			for ti = 1:numel(tagsB)
				options.(tagsB{ti}) = self.spl.(tagsB{ti});
			end
			if cndata.ValidFrames(1)<=self.spl.frrng(1) && ...
			   self.spl.frrng(2)<=cndata.ValidFrames(2) && ...
			   self.spl.frrng(1)<=seedframe && seedframe<=self.spl.frrng(2)
				options.FramesForAnalysis = self.spl.frrng;
			end
		end
	end
	

	% 'UnwrapConnectivity' option based on ROI type
	if any(strcmpi(cndata.ROIType,{'open','closed'}));
		options.UnwrapConnectivity = 8;
	end

	% initialize waitbar timer
	hwait = waitbartimer;
	cleanupObj = onCleanup(@()delete(hwait(isvalid(hwait))));
	hwait.String = 'Performing Analysis...';
	hwait.WindowStyle = 'modal';
	hwait.AllowClose  = false;
	hwait.start;

	% Perform analysis
	try
		self.spl = DENSEspline(imdata,cndata,options,...
			'SeedFrame',seedframe,'OptionsPanel',true, varargin{:});
	catch ERR
		hwait.stop;
		delete(hwait);
		errstr = sprintf('%s','Analysis error - more information will be ',...
			'printed to the command line.');
		h = errordlg(errstr,'Analysis error!','modal');
		waitfor(h);
		rethrow(ERR);
	end

	% remove waitbar timer
	hwait.stop;
	delete(hwait)

	if ~isempty(self.spl)
		self.spl.DENSEUID = obj.dns(didx).UID;
		self.spl.ROIUID   = obj.roi(ridx).UID;
		self.spl.Mag      = imdata.Mag;
		self.spl.DENSEType = imdata.DENSEType;
		self.spl.ROIType   = cndata.ROIType;

		if ~isempty(self.spl.Xseed)
			self.dns(didx).seed{1,self.spl.Xseed(end)} = self.spl.Xseed;
			self.dns(didx).seed{2,self.spl.Yseed(end)} = self.spl.Yseed;
			self.dns(didx).seed{3,self.spl.Zseed(end)} = self.spl.Zseed;
		end
		
		% notify(obj,'NewState',DENSEEventData('new','spl'));
	end
end