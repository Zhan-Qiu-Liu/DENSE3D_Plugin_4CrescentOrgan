function [file, out] = exportMatREPL(self,obj,startpath)
% Copyright (c) 2016 DENSEanalysis Contributors
% The codes below are modified from the function "exportMatFcn" in "AnalysisViewer.m".

	import plugins.DENSE3D_Plugin_4CrescentOrgan.*
%{ 
    % determine strain object
    if isempty(self.straindata)
        spl2strainFcnREPL(self,obj);
        if isempty(self.straindata)
            file = [];
            return;
        end
    end
 %}
    % check for startpath
    if nargin < 2 || isempty(startpath)
        startpath = pwd;
    end

    if ~exist(startpath, 'dir')
        [origstartpath,~,ext] = fileparts(startpath);
        file = startpath;
        startpath = origstartpath;
        if isempty(ext)
            % If there is no extension then MATLAB save will add
            % '.mat'. This will allow us to return the correct
            % filename to the user.
            file = [file '.mat'];
        end
    end

    % DENSE & ROI indices
    duid = self.spl.DENSEUID;
    didx = obj.hdata.UIDtoIndexDENSE(duid);

    ruid = self.spl.ROIUID;
    ridx = obj.hdata.UIDtoIndexROI(ruid);

    % If the startpath is not a mat file then make one
    if ~exist('file', 'var')
        % file name
        header = sprintf('%s_%s',...
                         self.dns(didx).Name,...
                         obj.hdata.roi(ridx).Name);

        expr = '[\\\/\?\%\:\*\"\<\>\|]';
        header = regexprep(header,expr,'_');

        file = fullfile(startpath,[header '.mat']);
        cnt = 0;
        while isfile(file)
            cnt = cnt+1;
            file = fullfile(startpath,sprintf('%s (%d).mat',header,cnt));
        end

        % allow user to change file name
        [uifile,uipath] = uiputfile('*.mat',[],file);
        if isequal(uifile,0)
            file = [];
            return;
        end

        % check extension
        file = fullfile(uipath,uifile);
        [~,~,e] = fileparts(file);
        if ~isequal(e,'.mat')
            file = [file, '.mat'];
        end
    end

    % start waitbartimer
    hwait = waitbartimer;
    cleanupObj = onCleanup(@()delete(hwait(isvalid(hwait))));
    hwait.String = 'Saving MAT file...';
    hwait.WindowStyle = 'modal';
    hwait.AllowClose  = false;
    start(hwait);
    drawnow

    % magnitude/phase indices
    midx = self.dns(didx).MagIndex;
    pidx = self.dns(didx).PhaIndex;
    sidx = [midx; pidx];

    % image information
    tags = {'DENSEType','Multipliers','Mag'};
    tfsingle = logical([0 0 0]);
    if ~isnan(pidx(1))
        tags = cat(2,tags,{'Xwrap','Xunwrap'});
        tfsingle = [tfsingle, true, true];
    end
    if ~isnan(pidx(2))
        tags = cat(2,tags,{'Ywrap','Yunwrap'});
        tfsingle = [tfsingle, true, true];
    end
    if ~isnan(pidx(3))
        tags = cat(2,tags,{'Zwrap','Zunwrap'});
        tfsingle = [tfsingle, true, true];
    end

    ImageInfo = struct;
    for ti = 1:numel(tags)
        tag = tags{ti};
        if tfsingle(ti)
            ImageInfo.(tag) = single(self.spl.(tag));
        else
            ImageInfo.(tag) = self.spl.(tag);
        end
    end

    % ROI information
    tags = {'ROIType','RestingContour','Contour'};
    ROIInfo = struct;
    for ti = 1:numel(tags)
        ROIInfo.(tags{ti}) = self.spl.(tags{ti});
    end

    % analysis information
    tags = {'ResampleMethod','ResampleDistance','SpatialSmoothing',...
        'TemporalOrder','Xseed','Yseed','Zseed'};
    AnalysisInfo = struct;
    for ti = 1:numel(tags)
        AnalysisInfo.(tags{ti}) = self.spl.(tags{ti});
    end

    % frame range analyized
    frrng  = self.spl.frrng;
    frames = frrng(1):frrng(2);
    AnalysisInfo.FramesForAnalysis = frrng;

    % segment model & orientation
    if any(strcmpi(self.spl.ROIType,{'sa','la'}))
		% save to *.mat
        AnalysisInfo.Nmodel    = self.straindata.Nmodel;
        AnalysisInfo.Nseg      = self.straindata.Nseg;
        try AnalysisInfo.SegDistribution      = self.straindata.SegDistribution; end
        AnalysisInfo.PositionA = self.straindata.PositionA;
        AnalysisInfo.PositionB = self.straindata.PositionB;
        AnalysisInfo.Clockwise = self.straindata.Clockwise;
    end

    % DENSE group information
    DENSEInfo = self.dns(didx);

    % original sequence header information
    SequenceInfo = repmat(struct,[2 3]);
    for k = 1:6
        if ~isnan(sidx(k))
            tags = fieldnames(obj.hdata.seq(sidx(k)));
            for ti = 1:numel(tags)
                tag = tags{ti};
                SequenceInfo(k).(tag) = ...
                    obj.hdata.seq(sidx(k)).(tag);
            end
        end
    end

    % Displacement Info
    Isz = size(self.spl.Xwrap(:,:,1));
    Nfr = size(self.spl.Xwrap,3);

    x = 1:Isz(2);
    y = 1:Isz(1);
    [X,Y] = meshgrid(x,y,0);

    mask0 = self.spl.MaskFcn(...
        X,Y,self.spl.RestingContour);
    Npts = sum(mask0(:));

    DisplacementInfo = struct(...
        'X',    X(mask0),...
        'Y',    Y(mask0),...
        'dX',   NaN([Npts,Nfr]),...
        'dY',   NaN([Npts,Nfr]),...
        'dZ',   NaN([Npts,Nfr]));

    pts = [Y(:),X(:),zeros(size(X(:)))];
    pts = pts(mask0,:)';

    for fr = frames
        pts(3,:) = fr;
        DisplacementInfo.dX(:,fr) = fnvalmod(self.spl.spldx,pts);
        DisplacementInfo.dY(:,fr) = fnvalmod(self.spl.spldy,pts);
        DisplacementInfo.dZ(:,fr) = fnvalmod(self.spl.spldz,pts);
    end

    % short-axis angle
    % for 6-segment model, [0,60)=anterior, [60,120)=anteroseptal, etc.
    % for 4-segment model, [0,90)=anterior, [90 180)=septal, etc.
    if strcmpi(self.spl.ROIType,'sa')
        origin  = self.straindata.PositionA;
        posB    = self.straindata.PositionB;
        flag_clockwise = self.straindata.Clockwise;

        theta0 = atan2(posB(2)-origin(2),posB(1)-origin(1));
        theta  = atan2(Y(mask0)-origin(2),X(mask0)-origin(1)) - theta0;
        if ~flag_clockwise, theta = -theta; end

        theta(theta<0) = theta(theta<0) + 2*pi;
        theta(theta>=2*pi) = 0;

        DisplacementInfo.Angle = theta(:);
    end

    % fv & strain tags
    if any(strcmpi(self.spl.ROIType,{'open','closed'}))
        fvtag = 'fv';
        sttag = 'strain';
    else
        fvtag = 'fvpix';
        sttag = 'strainpix';
    end

    % Strain data
    StrainInfo = struct(...
        'X',        X,...
        'Y',        Y,...
        'Faces',    self.straindata.(fvtag).faces,...
        'Vertices', self.straindata.(fvtag).vertices);

    % strain orientation
    if any(strcmpi(self.spl.ROIType,{'sa','la'}))
        StrainInfo.PolarOrientation = self.straindata.(fvtag).orientation;
    end

    % expand the maskimage
    if ~any(strcmpi(self.spl.ROIType,{'open','closed'}))
        irng = self.spl.irng;
        jrng = self.spl.jrng;

        mask = false(Isz);
        mask(irng(1):irng(2),jrng(1):jrng(2)) = ...
            self.straindata.(fvtag).maskimage;
        StrainInfo.Mask = mask;
    end

    % strain tag labels
    switch lower(self.spl.ROIType)
        case 'sa'
            tagA = {'XX','YY','XY','YX','RR','CC','RC','CR','p1','p2','p1or'};
            tagB = tagA;
        case 'la'
            tagA = {'XX','YY','XY','YX','RR','CC','RC','CR','p1','p2','p1or'};
            tagB = {'XX','YY','XY','YX','RR','LL','RL','LR','p1','p2','p1or'};
        case {'open','closed'}
            tagA = {'SS'};
            tagB = tagA;
        otherwise
            tagA = {'XX','YY','XY','YX','p1','p2','p1or'};
            tagB = tagA;
    end

    % copy relevant strain information
    for ti = 1:numel(tagA)
        StrainInfo.(tagB{ti}) = self.straindata.(sttag).(tagA{ti});
    end

    % export file
    AnalysisInstanceUID = dicomuid;

    out.ImageInfo = ImageInfo;
    out.ROIInfo = ROIInfo;
    out.DisplacementInfo = DisplacementInfo;
    out.AnalysisInfo = AnalysisInfo;
    out.DENSEInfo = DENSEInfo;
    out.SequenceInfo = SequenceInfo;
    out.StrainInfo = StrainInfo;
    out.TransmuralStrainInfo = computeTransmuralDataREPL(self);
    out.AnalysisInstanceUID = AnalysisInstanceUID;

    save(file, '-struct', 'out')
end