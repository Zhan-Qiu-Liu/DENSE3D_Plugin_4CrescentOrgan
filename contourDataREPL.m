function cndata = contourDataREPL(obj,ridx,frames,checkonlyflag)
% Copyright (c) 2016 DENSEanalysis Contributors
% The codes below are modified from the function "contourData" in "DENSEdata.m".

    if nargin < 4 || isempty(checkonlyflag)
        checkonlyflag = false;
    end

    % check for empty object
    if isempty(obj.seq) || isempty(obj.roi)
        error(sprintf('%s:invalidInput',mfilename),'%s',...
            'One or more object fields are empty.');
    end

    % check indices
    rng = [1 numel(obj.roi)];
    if ~isnumeric(ridx) || ~isscalar(ridx) || ...
       ridx < rng(1) || rng(2) < ridx
        error(sprintf('%s:invalidROIIndex',mfilename),...
            'ROI Index is invalid.');
    end

    % gather object
    roi  = obj.roi(ridx);
    pos  = obj.roi(ridx).Position;
    Nfr  = size(pos,1);

    % empty ROI frames positions
    tf = cellfun(@isempty,pos);
    emptypos = any(tf,2);

    % check for ROI continuity
    idx = find(~emptypos);
    if any(emptypos(idx(1):idx(end)))
        error(sprintf('%s:invalidPosition',mfilename),...
            'The ROI is not defined on a continuous frame range.');
    end

    % default inputs
    if nargin < 3 || isempty(frames), frames = idx; end

    % check frames
    if ~isnumeric(frames) || any(mod(frames,1)~=0) || ...
       any(frames < 1) || any(Nfr < frames)
        error(sprintf('%s:invalidFrame',mfilename),...
            'Frame Index is invalid.');
    end

    % check for valid positions
    if any(emptypos(frames))
        error(sprintf('%s:invalidFrame',mfilename),...
            'The ROI has not been defined on one or more frames.');
    end

    if checkonlyflag
        maskfcn = [];
        C = [];
    else
        ROI = findobj(obj.roitypes, 'Type', roi.Type);
        maskfcn = @(X,Y,C)ROI.mask(X,Y,C);

        % contour retrieval
        C = repmat({zeros(0,2)},size(pos));
        for fr = frames(:)'
            for k = 1:size(C,2)
                seg = clinesegments(pos{fr,k},...
                    roi.IsClosed{fr,k},roi.IsCurved{fr,k},...
                    roi.IsCorner{fr,k},0.5);
                crv = cat(1,seg{:});
                tf = [true; all(crv(2:end,:)~=crv(1:end-1,:),2)];
                C{fr,k} = crv(tf,:);
            end
        end
    end

    % save contour data
    cndata = struct(...
        'ROIType',      roi.Type,...
        'SeqIndex',     roi.SeqIndex,...
        'Contour',      {C},...
        'MaskFcn',      maskfcn,...
        'ValidFrames',  frames([1 end]));
end