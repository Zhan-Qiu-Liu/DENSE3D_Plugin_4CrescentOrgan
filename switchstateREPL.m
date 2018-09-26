function switchstateREPL(hfig)
% Copyright (c) 2016 DENSEanalysis Contributors
% The codes below are modified from the function "switchstate" in "DENSEanalysis.m".

%% SWITCH TAB (DICOM/DENSE/ANALYSIS)

    handles = guidata(hfig);

    % Get the indices of the sidebar tabs
    sidebarinds = handles.hsidebar.find([
        handles.dicom_hpanel
        handles.dense_hpanel
        handles.analysis_hpanel]);

    % current tab
    tabidx = handles.hsidebar.ActiveTab;
    handles.renderer{handles.LastTab} = get(handles.hfig, 'renderer');

    if ~ismember(tabidx, sidebarinds)
        return
    end

    set(handles.hfig, 'Colormap', gray)

    wild = {'*zoom*','*pan','*rotate*','tool*','*save*','*contrast*'};
    tf = cellfun(@(tag)any(strwcmpi(tag,wild)),get(handles.htools,'tag'));
    set(handles.htools(tf),'Enable','on');

    handles.hdicom.ROIEdit = 'off';
    handles.hdense.ROIEdit = 'off';
    set(handles.tool_roi,'State','off','enable','on');

    popupinds = handles.hpopup.find([
        handles.popup_dicom
        handles.popup_arial
        handles.popup_dense
        handles.popup_slice
        handles.popup_analysis
    ]);

    switch tabidx

        case 1

            tf = logical([1 0 1 1 0]);
            handles.hpopup.Visible(popupinds(tf))  = {'on'};
            handles.hpopup.Visible(popupinds(~tf)) = {'off'};

            handles.hanalysis.Enable = 'off';

            % transfer slice/arial to DICOM viewer
            handles.hdense.SliceViewer = [];
            handles.hdense.ArialViewer = [];
            handles.hdicom.SliceViewer = handles.hslice;
            handles.hdicom.ArialViewer = handles.harial;
            redraw(handles.hdicom);

        case 2
            tf = logical([0 1 1 1 0]);
            handles.hpopup.Visible(popupinds(tf))  = {'on'};
            handles.hpopup.Visible(popupinds(~tf)) = {'off'};

            handles.hanalysis.Enable = 'off';

            % transfer slice/arial to DENSE viewer
            handles.hdicom.SliceViewer = [];
            handles.hdicom.ArialViewer = [];
            handles.hdense.SliceViewer = handles.hslice;
            handles.hdense.ArialViewer = handles.harial;

            redraw(handles.hdense);

        case 3

            tf = logical([0 0 0 0 1]);
            handles.hpopup.Visible(popupinds(tf))  = {'on'};
            handles.hpopup.Visible(popupinds(~tf)) = {'off'};
            handles.hanalysis.Enable = 'on';
            set(handles.tool_roi,'enable','off');
            redraw(handles.hanalysis);

    end


    set(handles.hfig,'renderer',handles.renderer{tabidx});
    handles.LastTab = tabidx;
    guidata(hfig,handles);

end
