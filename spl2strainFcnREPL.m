function data = spl2strainFcnREPL(self,didx)%,obj
% Copyright (c) 2016 DENSEanalysis Contributors
% The codes below are modified from the function "spl2strainFcn" in "AnalysisViewer.m".

    % waitbar
    hwait = waitbartimer;
    cleanupObj = onCleanup(@()delete(hwait(isvalid(hwait))));
    hwait.String = 'Calculating strain & displacement patterns';
    hwait.WindowStyle = 'modal';
    hwait.AllowClose = false;
    hwait.start;

    % ROI type
    type = self.spl.ROIType;

    % frames for analysis
    Nfr    = size(self.spl.Mag,3);
    frrng  = self.spl.frrng;
    frames = frrng(1):frrng(2);

	
	api = struct(...
		'meshCtrl',            	self.dns(didx).meshCtrl,...
		'res',              	1,...
		'InstallDir',           self.InstallDir,...
		'Type',             	type,...
		'Mag',                  self.spl.Mag,...
		'RestingContour',   	{self.spl.RestingContour},...
		'FramesForAnalysis',	frrng,...
		'spldx',            	self.spl.spldx,...
		'spldy',            	self.spl.spldy,...
		'spldz',            	self.spl.spldx,...
		'MaskFcn',          	self.spl.MaskFcn,...
		'xrng',             	self.spl.jrng,...
		'yrng',             	self.spl.irng,...
		'Resolution',           1,...
		'PositionA',            [],...
		'PositionB',            self.dns(didx).SegPos,...
		'Nmodel',               [],...
		'Nseg',                 [],...
		'Clockwise',            [],...
		'PositionIndices',      [],...
		'SegmentModelPanel',    false);

    if strcmpi(type,'SA')
        % if ~isempty(self.straindata)
            % api.PositionA = self.straindata.PositionA;
            % api.PositionB = self.straindata.PositionB;
            % api.Nmodel    = self.straindata.Nmodel;
            % api.Nseg      = self.straindata.Nseg;
            % api.Clockwise = self.straindata.Clockwise;
        % end
		
		import plugins.dense3D_plugin_4crescentorgan.*
        opts = MeshControl(api);
        drawnow
        if isempty(opts); data = []; return; end
        tags = fieldnames(opts);
        for ti = 1:numel(tags)
            api.(tags{ti}) = opts.(tags{ti});
        end
    elseif strcmpi(type,'LA')
        api.Mag  = self.spl.Mag;
        api.SegmentModelPanel = true;
        % if ~isempty(self.straindata)
            % api.PositionA = self.straindata.PositionA;
            % api.PositionB = self.straindata.PositionB;
            % api.Nmodel    = self.straindata.Nmodel;
            % api.Nseg      = self.straindata.Nseg;
            % api.Clockwise = self.straindata.Clockwise;
        % end
    elseif any(strcmpi(type,{'open','closed'}))
        api.Mag  = self.spl.Mag;
        api.SegmentModelPanel = true;
        if ~isempty(self.straindata)
            api.PositionIndices = self.straindata.PositionIndices;
        end
    end
	
    % strain data
    data = spl2strain(api);
    if isempty(data), return; end

    % append frame range to strain data
    data.frrng = frrng;

    % gather additional data
    strtypes = {'','pix'};
    dsptypes = {'X','Y','Z'};

    for si = 1:numel(strtypes)

        % face/vertex & strain tags
        fvtag = ['fv' strtypes{si}];
        sttag = ['strain' strtypes{si}];

        if ~isfield(data,fvtag), continue; end

        % face centroids (trajectory origins)
        xv = data.(fvtag).vertices(:,1);
        yv = data.(fvtag).vertices(:,2);
        f  = data.(fvtag).faces;

        ori = [mean(xv(f),2),mean(yv(f),2)];
        Nx  = size(ori,1);

        % trajectories
        pos = ori(:,[2 1])';
        for di = 1:numel(dsptypes)
            tagA = ['d' upper(dsptypes{di})];
            tagB = ['spld' lower(dsptypes{di})];

            tmp = zeros(Nx,numel(frames));
            for fridx = 1:numel(frames)
                pos(3,:) = frames(fridx);
                tmp(:,fridx) = fnvalmod(self.spl.(tagB),pos);
            end
            data.(sttag).(tagA) = tmp;
        end

        % polar displacement & twist
        if any(strcmpi(type,{'SA'}))

            % polar origins (bulk corrected)
            X0 = ori(:,1) - mean(ori(:,1));
            Y0 = ori(:,2) - mean(ori(:,2));
            [t0,r0] = cart2pol(X0,Y0);

            % polar trajectories (bulk corrected)
            X = bsxfun(@plus,ori(:,1),data.(sttag).dX);
            X = bsxfun(@minus,X,mean(X,1));
            Y = bsxfun(@plus,ori(:,2),data.(sttag).dY);
            Y = bsxfun(@minus,Y,mean(Y,1));
            [t,r] = cart2pol(X,Y);

            % unwrap trajectory angle
            tall = [t0,t];
            tall = unwrap(tall,[],2);
            t = tall(:,2:end);

            % change in radius/angle
            dR = bsxfun(@minus,r,r0);
            dT = bsxfun(@minus,t,t0);

            % change angle direction if necessary
            if ~data.Clockwise
                dT = -dT;
            end

            % arc length (circumferential displacement)
            dC = bsxfun(@times,dT,r0);

            % save polar information
            data.(sttag).dR = dR;
            data.(sttag).dC = dC;
            data.(sttag).twist = dT*180/pi;
        end

        % extend all strain data to size [Nx x Nfr]
        tags = fieldnames(data.(sttag));
        for ti = 1:numel(tags)
            tmp = zeros([Nx,Nfr]);
            tmp(:,frames) = data.(sttag).(tags{ti});
            data.(sttag).(tags{ti}) = tmp;
        end
    end

    % save to object
    data.meshCtrl = api.meshCtrl;	
    try data.SegDistribution = api.SegDistribution;	end
	data.PositionB = api.PositionB;
end