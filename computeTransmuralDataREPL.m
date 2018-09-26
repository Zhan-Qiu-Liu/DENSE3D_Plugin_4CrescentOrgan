function tdata = computeTransmuralDataREPL(self)
% Copyright (c) 2016 DENSEanalysis Contributors
% The codes below are modified from the function "computeTransmuralData" in "AnalysisViewer.m".

    % Compute transmural strain values
    straintags  = {'RR','CC','p1','p2'};
    strainnames = {'Err','Ecc','E1','E2'};

    switch lower(self.spl.ROIType)
        case 'la'
            strainnames{2} = 'Ell';
        case 'sa'
            straintags{end+1} = 'twist';
            strainnames{end+1} = 'Twist';
    end

    % layers to report
    layernames = {'mid','subepi','subendo','ave'};
    layerids = {3, 1, 5, 1:5};

    % Initialize the output structure
    fields = cat(1,layernames, cell(size(layernames)));
    substruct = struct(fields{:});

    fields = cat(1,strainnames,repmat({substruct},size(strainnames)));
    tdata = struct(fields{:});

    [i,j]  = ndgrid(1:numel(strainnames),1:numel(layernames));
    tags  = straintags(i);
    lids  = layerids(j);
    types = strainnames(i);
    trans = layernames(j);

    Nsec   = self.straindata.Nseg;
    Nfr    = size(self.straindata.strain.RR,2);
    frrng  = self.straindata.frrng;
    frames = frrng(1):frrng(2);

    layer  = self.straindata.fv.layerid;
    sector = self.straindata.fv.sectorid;

    for k = numel(tags):-1:1
        lid  = lids{k};
        tag  = tags{k};
        data = zeros(Nfr,Nsec);

        ltf = ismember(layer,lid);

        for sid = 1:Nsec
            stf = (sector == sid);
            tf = ltf & stf;

            for fr = frames
                data(fr,sid) = mean(self.straindata.strain.(tag)(tf,fr));
            end

        end

        tdata.(types{k}).(trans{k}) = data(frames,:);
    end

    % CURE/RURE output
    if isequal(lower(self.spl.ROIType),'sa')

        RR = zeros(Nfr,Nsec);
        CC = RR;
        for sid = 1:Nsec
            tf = (sector == sid);
            for fr = frames
                RR(fr,sid) = mean(self.straindata.strain.RR(tf,fr));
                CC(fr,sid) = mean(self.straindata.strain.CC(tf,fr));
            end
        end

        valcure = CURE(CC(frames,:));
        valrure = CURE(RR(frames,:));

        tdata.CURE = valcure;
        tdata.RURE = valrure;
    end
end