function handles = resetFcnREPL(hfig)
% Copyright (c) 2016 DENSEanalysis Contributors
% The codes below are modified from the function "resetFcn" in "DENSEanalysis.m".
% this function ensures the state of the GUI matches the data loaded into
% the DENSEdata object (handles.hdata)

	import plugins.DENSE3D_Plugin_4CrescentOrgan.*
	
    % gather gui data
    handles = guidata(hfig);

    % disable zoom/pan/rot
    handles.hzoom.Enable = 'off';
    handles.hpan.Enable  = 'off';
    handles.hrot.Enable  = 'off';

    % reset tabs to initial state: got to be #1 since no switchstate func
    handles.hsidebar.ActiveTab = 2;

    handles.hsidebar.Enable();

    % inds = handles.hsidebar.find([handles.dense_hpanel,handles.analysis_hpanel]);
    % handles.hsidebar.Enable(inds)   = {'off'};

    H = [handles.popup_dicom, handles.popup_dense, handles.popup_arial, ...
         handles.popup_slice, handles.popup_analysis];
    inds = handles.hpopup.find(H);

    handles.hpopup.Visible(:) = {'off'};
    handles.hpopup.Visible(inds) = {'on'};
    handles.hpopup.Enable(:) = {'off'};

    % reset renderer
    handles.renderer(:) = {'painters'};
    handles.LastTab = 2;

    % disable toolbar tools
    set(handles.htools,'Enable','off');
    wild = {'*open','*new'};
    tf = cellfun(@(tag)any(strwcmpi(tag,wild)),get(handles.htools,'tag'));
    set(handles.htools(tf),'Enable','on');


    % quit if no data to display
    if numel(handles.hdata.seq) == 0
        handles.hsidebar.redraw();
        handles.hpopup.redraw();
        return
    end

    % enable/open popups
    handles.hpopup.Enable(inds) = {'on'};
    handles.hpopup.IsOpen(inds) = 1;

    % enable sidetabs
    inds = handles.hsidebar.find([handles.dense_hpanel,handles.analysis_hpanel]);
    handles.hsidebar.Enable(inds) = {'on'};

    % activate/deactivate stuff
    switchstateREPL(handles.hfig)

    % trigger resize functions
    % (panels are resized via "hsidebar" object)
    handles.hsidebar.redraw();
    handles.hpopup.redraw();
end
