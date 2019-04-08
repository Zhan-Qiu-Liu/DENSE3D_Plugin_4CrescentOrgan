%
% File: Adjust_Mesh.m
% Author: Oskar Skrinjar
% Date: May 2014
% ----------------------------------------------------------------------
% This function adjusts the mesh input mesh to better match the points.
%

function mesh = Adjust_Mesh(mesh, points, APPEARANCE, name)

    % method parameters
    initial_NOR_R = 10; % [mm]
    final_NOR_R = 5; % [mm]
    initial_TAN_R = 5; % [mm]
    final_TAN_R = 5; % [mm]
    initial_smoothness_weight = .9;
    final_smoothness_weight = .9;
    initial_points_weight = .1;
    final_points_weight = .2;
    number_of_iterations = 100;

    % COLOR PARAMETERS
    STATE.SURFACE_COLOR = 'y';
    STATE.EDGE_COLOR = 'k';
    STATE.SMOOTH_DISP_COLOR = 'm';
    STATE.POINT_DISP_COLOR = 'b';
    POINTS_COLOR = 'r.';
    INITIAL_SURFACE_ALPHA = .9;
    
    
    % hidden figure
    STATE.hidden_figure = figure('Visible', 'off'); % its purpose is to pass mesh structure (in its 'UserData' field) once the main figure is closed
    
    % save mesh, points, and name
    STATE.initial_mesh = mesh;
    STATE.points = points;
    STATE.name = name;
    
    % initialize handles
    STATE.surface_handle = [];
    STATE.smooth_disp_handle = [];
    STATE.point_disp_handle = [];
    

    % pre-compute node neighborhoods
    if ~isfield(mesh,'node_neighborhood'); mesh.node_neighborhood = Node_Neighborhoods(mesh); end
	STATE.nodes = mesh.node_neighborhood;
    
    
    %
    % Main figure
    %  
    STATE.fig = figure('Units', 'normalized', 'Position', [.05 .1 .9 .8], 'Name', sprintf('Generate %s mesh', STATE.name), 'NumberTitle', 'off', 'Menubar', 'none', 'Toolbar', 'none', 'Color', APPEARANCE.FIGURE_BACKGROUND_COLOR);

    
    % control panel
    PANEL_WIDTH = .15;
    panel_left = 1-PANEL_WIDTH+.01;
    panel_width = PANEL_WIDTH-.02;
    
    %
    % "GUI" panel
    %
    GUI_panel = uipanel('Parent', STATE.fig, 'Units', 'normalized', 'Position', [panel_left .01 panel_width .98], 'Title', '', 'BackgroundColor', APPEARANCE.PANEL_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);
    

    %
    % "Method" panel
    %
    method_panel = uipanel('Parent', GUI_panel, 'Units', 'normalized', 'Position', [.01 .31 .98 .68], 'Title', 'Method', 'BackgroundColor', APPEARANCE.PANEL_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);
    
    %
    % "Parameters" panel
    %
    parameters_panel = uipanel('Parent', method_panel, 'Units', 'normalized', 'Position', [.01 .43 .98 .56], 'Title', 'Parameters', 'BackgroundColor', APPEARANCE.PANEL_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);
    
    % nor R label and edits
    axes('Parent', parameters_panel, 'Position', [.01 .81 .48 .18], 'Color', 'none', 'Visible', 'off');
    text(1,.5, 'nor R [mm]:', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    STATE.nor_R_1_edit = uicontrol('Parent', parameters_panel, 'Style', 'edit', 'Units', 'normalized', 'position', [.51 .81 .23 .18], 'String', sprintf('%g', initial_NOR_R));
    STATE.nor_R_2_edit = uicontrol('Parent', parameters_panel, 'Style', 'edit', 'Units', 'normalized', 'position', [.76 .81 .23 .18], 'String', sprintf('%g', final_NOR_R));
    
    % tan R label and edits
    axes('Parent', parameters_panel, 'Position', [.01 .61 .48 .18], 'Color', 'none', 'Visible', 'off');
    text(1,.5, 'tan R [mm]:', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    STATE.tan_R_1_edit = uicontrol('Parent', parameters_panel, 'Style', 'edit', 'Units', 'normalized', 'position', [.51 .61 .23 .18], 'String', sprintf('%g', initial_TAN_R));
    STATE.tan_R_2_edit = uicontrol('Parent', parameters_panel, 'Style', 'edit', 'Units', 'normalized', 'position', [.76 .61 .23 .18], 'String', sprintf('%g', final_TAN_R));
    
    % smooth label and edits
    axes('Parent', parameters_panel, 'Position', [.01 .41 .48 .18], 'Color', 'none', 'Visible', 'off');
    text(1,.5, 'smooth:', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    STATE.smooth_1_edit = uicontrol('Parent', parameters_panel, 'Style', 'edit', 'Units', 'normalized', 'position', [.51 .41 .23 .18], 'String', sprintf('%g', initial_smoothness_weight));
    STATE.smooth_2_edit = uicontrol('Parent', parameters_panel, 'Style', 'edit', 'Units', 'normalized', 'position', [.76 .41 .23 .18], 'String', sprintf('%g', final_smoothness_weight));
    
    % points label and edits
    axes('Parent', parameters_panel, 'Position', [.01 .21 .48 .18], 'Color', 'none', 'Visible', 'off');
    text(1,.5, 'points:', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    STATE.points_1_edit = uicontrol('Parent', parameters_panel, 'Style', 'edit', 'Units', 'normalized', 'position', [.51 .21 .23 .18], 'String', sprintf('%g', initial_points_weight));
    STATE.points_2_edit = uicontrol('Parent', parameters_panel, 'Style', 'edit', 'Units', 'normalized', 'position', [.76 .21 .23 .18], 'String', sprintf('%g', final_points_weight));
    
    %  iterations label and edit
    axes('Parent', parameters_panel, 'Position', [.01 .01 .68 .18], 'Color', 'none', 'Visible', 'off');
    text(1,.5, 'iterations:', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    STATE.iterations_edit = uicontrol('Parent', parameters_panel, 'Style', 'edit', 'Units', 'normalized', 'position', [.71 .01 .28 .18], 'String', sprintf('%d', number_of_iterations));
    
    
    
    %
    % "Status" panel
    %
    status_panel = uipanel('Parent', method_panel, 'Units', 'normalized', 'Position', [.01 .29 .98 .12], 'Title', 'Status', 'BackgroundColor', APPEARANCE.PANEL_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);
    
    % "iteration label"
    axes('Parent', status_panel, 'Position', [.05 .2 .9 .6], 'Color', 'none', 'Visible', 'off');
    STATE.state_label = text(0,.5, 'iteration:', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    
    %
    % "Controls" panel
    %
    controls_panel = uipanel('Parent', method_panel, 'Units', 'normalized', 'Position', [.01 .01 .98 .26], 'Title', 'Controls', 'BackgroundColor', APPEARANCE.PANEL_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);
    
    % reset button
    reset_pushbutton = uicontrol('Parent', controls_panel, 'Style', 'pushbutton', 'Units', 'normalized', 'position', [.05 .6 .4 .3], 'String', 'Reset', 'BackgroundColor', APPEARANCE.GUI_BACKGROUND_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);

    % run button
    run_pushbutton = uicontrol('Parent', controls_panel, 'Style', 'pushbutton', 'Units', 'normalized', 'position', [.55 .6 .4 .3], 'String', 'Run', 'BackgroundColor', APPEARANCE.GUI_BACKGROUND_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);

    % update displacements button
    update_displacements_pushbutton = uicontrol('Parent', controls_panel, 'Style', 'pushbutton', 'Units', 'normalized', 'position', [.05 .1 .9 .3], 'String', 'Update Displacements', 'BackgroundColor', APPEARANCE.GUI_BACKGROUND_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);

    
    
    %
    % "Display" panel
    %
    display_panel = uipanel('Parent', GUI_panel, 'Units', 'normalized', 'Position', [.01 .11 .98 .18], 'Title', 'Visibility', 'BackgroundColor', APPEARANCE.PANEL_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);
    
    % surface radiobutton
    STATE.surface_checkbox = uicontrol('Parent', display_panel, 'Style', 'checkbox', 'Units', 'normalized', 'position', [.05 .7 .4 .2], 'String', 'Surface', 'Value', 1, 'BackgroundColor', APPEARANCE.PANEL_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);
    
    % surface transparency slider
    STATE.surface_transparency_slider = uicontrol('Parent', display_panel, 'Style', 'slider', 'Units', 'normalized', 'position', [.51 .7 .48 .2], 'Min', 0, 'Max', 1, 'Value', INITIAL_SURFACE_ALPHA, 'BackgroundColor', APPEARANCE.PANEL_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);
    
    % smooth displacements radiobutton
    STATE.smooth_disp_checkbox = uicontrol('Parent', display_panel, 'Style', 'checkbox', 'Units', 'normalized', 'position', [.05 .4 .9 .2], 'String', 'Smooth displacements', 'Value', 1, 'BackgroundColor', APPEARANCE.PANEL_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);
    
    % surface radiobutton
    STATE.point_disp_checkbox = uicontrol('Parent', display_panel, 'Style', 'checkbox', 'Units', 'normalized', 'position', [.05 .1 .9 .2], 'String', 'Point displacements', 'Value', 1, 'BackgroundColor', APPEARANCE.PANEL_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);
    
    
    
    %
    % "Misc" panel
    %
    misc_panel = uipanel('Parent', GUI_panel, 'Units', 'normalized', 'Position', [.01 .01 .98 .08], 'Title', 'Misc', 'BackgroundColor', APPEARANCE.PANEL_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);
    
    % accept button
    accept_pushbutton = uicontrol('Parent', misc_panel, 'Style', 'pushbutton', 'Units', 'normalized', 'position', [.05 .2 .4 .6], 'String', 'Accept', 'BackgroundColor', APPEARANCE.GUI_BACKGROUND_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);

    % help button
    help_pushbutton = uicontrol('Parent', misc_panel, 'Style', 'pushbutton', 'Units', 'normalized', 'position', [.55 .2 .4 .6], 'String', 'Help', 'BackgroundColor', APPEARANCE.GUI_BACKGROUND_COLOR, 'ForegroundColor', APPEARANCE.TEXT_COLOR);

    
    %
    % "3D" panel
    %
    TD_panel = uipanel('Parent', STATE.fig, 'Units', 'normalized', 'Position', [.01 .01 panel_left-.02 .98], 'Title', '');
    
    
    %
    % 3D Axes
    %
    STATE.main_axes = axes('Parent', TD_panel, 'Position', [0 0 1 1], 'Color', 'none', 'Visible', 'off');
    hold on;
    
    
    % display points
    plot3(points(:,1), points(:,2), points(:,3), POINTS_COLOR);
    
    
    % set lighting
    light;
    
    % set the axis
    hold on;
    axis equal;
    axis vis3d;
    set(STATE.main_axes, 'Projection', 'perspective');
    axis off;
    
    % set the cameratoolbar
    cameratoolbar(STATE.fig);
    cameratoolbar('SetMode', 'orbit');
    
    % set the renderer
    set(STATE.fig, 'Renderer', 'OpenGL');
    
    
    % help widnow
    help_str = { sprintf('This windows allows the user to fit an %s surface mesh to the slice boundary (red) points.', STATE.name) };
    help_str = [ help_str ' ' '"Parameters" panel allows the user to set the starting and ending parameter values that are lineraly varied (from the starting to the ending values) over the iterations. The parameters are "normal R", "tagent R", "smoothing weight", "points weight" and the number of interations.' ];
    help_str = [ help_str ' ' '"Status" panel shows the current iteration.' ];
    help_str = [ help_str ' ' '"Controls" panel allows the user to reset the surface to the initial shape, to run the method with the selected parameters or to update the displacements with the starting parameters.' ];
    help_str = [ help_str ' ' '"Visibility" panel controls the visibility of the mesh, smooth and point displacements. In addition, the slider controls the transparency of the surface mesh.' ];
    help_str = [ help_str ' ' 'If you are satisfied with the generated mesh you can accepte it by pressing "Accept" button. Otherwise keep adjusting the mesh or terminate the application by pressing "x" on the window bar.' ];
    help_position = [ .2 .35 .4 .3 ];
    APPEARANCE.FIGURE_BACKGROUND_COLOR = [ 231 234 241 ] / 255; % R,G,B values of the figure background color
	APPEARANCE.TEXT_COLOR = [ 0 0 0 ] / 255; % R,G,B values of the text color;
	try 
		[STATE.help_fig,~] = Help_Window(help_str, 'Help', APPEARANCE, help_position);
    catch
		STATE.help_fig = msgbox(sprintf('%s\n%s\n%s\n%s\n%s',help_str{:}));
	end
    
    % save STATE
    set(STATE.fig, 'UserData', STATE);
    
    % callbacks
    set(accept_pushbutton, 'Callback', {@Accept_CB, STATE.fig});
    set(help_pushbutton, 'Callback', {@Help_CB, STATE.fig});
    set(reset_pushbutton, 'Callback', {@Reset_CB, STATE.fig});
    set(run_pushbutton, 'Callback', {@Run_CB, STATE.fig});
    set(STATE.surface_checkbox, 'Callback', {@Visibility_CB, STATE.fig});
    set(STATE.smooth_disp_checkbox, 'Callback', {@Visibility_CB, STATE.fig});
    set(STATE.point_disp_checkbox, 'Callback', {@Visibility_CB, STATE.fig});
    set(update_displacements_pushbutton, 'Callback', {@Update_Displacements_CB, STATE.fig});
    set(STATE.surface_transparency_slider, 'Callback', {@Surface_Transparency_Slider_CB, STATE.fig});
    set(STATE.nor_R_1_edit, 'Callback', {@Nor_R_1_CB, STATE.fig});
    set(STATE.nor_R_2_edit, 'Callback', {@Nor_R_2_CB, STATE.fig});
    set(STATE.tan_R_1_edit, 'Callback', {@Tan_R_1_CB, STATE.fig});
    set(STATE.tan_R_2_edit, 'Callback', {@Tan_R_2_CB, STATE.fig});
    set(STATE.smooth_1_edit, 'Callback', {@Smooth_1_CB, STATE.fig});
    set(STATE.smooth_2_edit, 'Callback', {@Smooth_2_CB, STATE.fig});
    set(STATE.points_1_edit, 'Callback', {@Points_1_CB, STATE.fig});
    set(STATE.points_2_edit, 'Callback', {@Points_2_CB, STATE.fig});
    set(STATE.iterations_edit, 'Callback', {@Iterations_CB, STATE.fig});
    set(STATE.fig, 'CloseRequestFcn', {@Terminate_CB, STATE.fig});
    
    % start the process
    Reset_CB([], [], STATE.fig);
    
    % wait for the user to accept the mesh or cancel
    uiwait(STATE.fig);

    % set the output argument
    mesh = get(STATE.hidden_figure, 'UserData');
    
    % delete the hidden figure
    delete(STATE.hidden_figure);

end % Adjust_Mesh


%
% Function: Delete_Figures
%
function Delete_Figures(STATE)

    % delete the help figure
    try delete(STATE.help_fig); end
    
    % delete the main figure
    delete(STATE.fig);
    
end % Delete_Figures


%
% Function: Nor_R_2_CB
%
function Nor_R_2_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % check if the value is proper
    Read_Nor_R_2(STATE);
    
end % Nor_R_2_CB


%
% Function: Read_Nor_R_2
%
function nor_R_2 = Read_Nor_R_2(STATE)

    % read nor_R_2
    nor_R_2 = str2double(get(STATE.nor_R_2_edit, 'String'));

    % validity check
    if isnan(nor_R_2)
        uiwait(warndlg('Invalid final nor R!', 'Warning', 'modal'));
        nor_R_2 = [];
        return;
    end
    
    % validity check
    if (nor_R_2 < 0)
        uiwait(warndlg('Invalid final nor R: it needs to be a non-negative number!', 'Warning', 'modal'));
        nor_R_2 = [];
        return;
    end

end % Read_Nor_R_2


%
% Function: Tan_R_2_CB
%
function Tan_R_2_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % check if the value is proper
    Read_Tan_R_2(STATE);
    
end % Tan_R_2_CB


%
% Function: Read_Tan_R_2
%
function tan_R_2 = Read_Tan_R_2(STATE)

    % read tan_R_2
    tan_R_2 = str2double(get(STATE.tan_R_2_edit, 'String'));

    % validity check
    if isnan(tan_R_2)
        uiwait(warndlg('Invalid final tan R!', 'Warning', 'modal'));
        tan_R_2 = [];
        return;
    end
    
    % validity check
    if (tan_R_2 < 0)
        uiwait(warndlg('Invalid final tan R: it needs to be a non-negative number!', 'Warning', 'modal'));
        tan_R_2 = [];
        return;
    end

end % Read_Tan_R_2


%
% Function: Smooth_2_CB
%
function Smooth_2_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % check if the value is proper
    Read_Smooth_2(STATE);
    
end % Smooth_2_CB


%
% Function: Read_Smooth_2
%
function smooth_2 = Read_Smooth_2(STATE)

    % read smooth_1
    smooth_2 = str2double(get(STATE.smooth_2_edit, 'String'));

    % validity check
    if isnan(smooth_2)
        uiwait(warndlg('Invalid final smooth weight!', 'Warning', 'modal'));
        smooth_2 = [];
        return;
    end
    
    % validity check
    if ( (smooth_2 < 0) || (smooth_2 > 1) )
        uiwait(warndlg('Invalid final smooth weight: it needs to be a value from [0,1]!', 'Warning', 'modal'));
        smooth_2 = [];
        return;
    end

end % Read_Smooth_2


%
% Function: Points_2_CB
%
function Points_2_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % check if the value is proper
    Read_Points_2(STATE);
    
end % Points_2_CB


%
% Function: Read_Points_2
%
function points_2 = Read_Points_2(STATE)

    % read points_2
    points_2 = str2double(get(STATE.points_2_edit, 'String'));

    % validity check
    if isnan(points_2)
        uiwait(warndlg('Invalid final points weight!', 'Warning', 'modal'));
        points_2 = [];
        return;
    end
    
    % validity check
    if ( (points_2 < 0) || (points_2 > 1) )
        uiwait(warndlg('Invalid final points weight: it needs to be a value from [0,1]!', 'Warning', 'modal'));
        points_2 = [];
        return;
    end

end % Read_Points_2


%
% Function: Iterations_CB
%
function Iterations_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % check if the value is proper
    Read_Iterations(STATE);
    
end % Iterations_CB


%
% Function: Read_Iterations
%
function iterations = Read_Iterations(STATE)

    % read iterations
    iterations = str2double(get(STATE.iterations_edit, 'String'));

    % validity check
    if isnan(iterations)
        uiwait(warndlg('Invalid number of iterations!', 'Warning', 'modal'));
        iterations = [];
        return;
    end
    
    % validity check
    if ( (iterations < 1) || (round(iterations) ~= iterations) )
        uiwait(warndlg('Invalid number of iterations: it needs to be a positive integer!', 'Warning', 'modal'));
        iterations = [];
        return;
    end

end % Read_Iterations


%
% Function: Points_1_CB
%
function Points_1_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % check if the value is proper
    Read_Points_1(STATE);
    
end % Points_1_CB


%
% Function: Read_Points_1
%
function points_1 = Read_Points_1(STATE)

    % read points_1
    points_1 = str2double(get(STATE.points_1_edit, 'String'));

    % validity check
    if isnan(points_1)
        uiwait(warndlg('Invalid initial points weight!', 'Warning', 'modal'));
        points_1 = [];
        return;
    end
    
    % validity check
    if ( (points_1 < 0) || (points_1 > 1) )
        uiwait(warndlg('Invalid initial points weight: it needs to be a value from [0,1]!', 'Warning', 'modal'));
        points_1 = [];
        return;
    end

end % Read_Points_1


%
% Function: Smooth_1_CB
%
function Smooth_1_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % check if the value is proper
    Read_Smooth_1(STATE);
    
end % Smooth_1_CB


%
% Function: Read_Smooth_1
%
function smooth_1 = Read_Smooth_1(STATE)

    % read smooth_1
    smooth_1 = str2double(get(STATE.smooth_1_edit, 'String'));

    % validity check
    if isnan(smooth_1)
        uiwait(warndlg('Invalid initial smooth weight!', 'Warning', 'modal'));
        smooth_1 = [];
        return;
    end
    
    % validity check
    if ( (smooth_1 < 0) || (smooth_1 > 1) )
        uiwait(warndlg('Invalid initial smooth weight: it needs to be a value from [0,1]!', 'Warning', 'modal'));
        smooth_1 = [];
        return;
    end

end % Read_Smooth_1


%
% Function: Nor_R_1_CB
%
function Nor_R_1_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % check if the value is proper
    Read_Nor_R_1(STATE);
    
end % Nor_R_1_CB


%
% Function: Read_Nor_R_1
%
function nor_R_1 = Read_Nor_R_1(STATE)

    % read nor_R_1
    nor_R_1 = str2double(get(STATE.nor_R_1_edit, 'String'));

    % validity check
    if isnan(nor_R_1)
        uiwait(warndlg('Invalid initial nor R!', 'Warning', 'modal'));
        nor_R_1 = [];
        return;
    end
    
    % validity check
    if (nor_R_1 < 0)
        uiwait(warndlg('Invalid initial nor R: it needs to be a non-negative number!', 'Warning', 'modal'));
        nor_R_1 = [];
        return;
    end

end % Read_Nor_R_1


%
% Function: Tan_R_1_CB
%
function Tan_R_1_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % check if the value is proper
    Read_Tan_R_1(STATE);
    
end % Tan_R_1_CB


%
% Function: Read_Tan_R_1
%
function tan_R_1 = Read_Tan_R_1(STATE)

    % read tan_R_1
    tan_R_1 = str2double(get(STATE.tan_R_1_edit, 'String'));

    % validity check
    if isnan(tan_R_1)
        uiwait(warndlg('Invalid initial tan R!', 'Warning', 'modal'));
        tan_R_1 = [];
        return;
    end
    
    % validity check
    if (tan_R_1 < 0)
        uiwait(warndlg('Invalid initial tan R: it needs to be a non-negative number!', 'Warning', 'modal'));
        tan_R_1 = [];
        return;
    end

end % Read_Tan_R_1



%
% Function: Reset_CB
%
function Reset_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % set the current mesh
    STATE.current_mesh = STATE.initial_mesh;
    
    % set the iteration
    STATE.iteration = 0;
    
    % display the surface
    STATE = Display_Surface(STATE);
    
    % save STATE
    set(fig, 'UserData', STATE);
    
    % update displacements
    Update_Displacements_CB([], [], STATE.fig);
    
end % Reset_CB


%
% Function: Display_Displacements
%
function STATE = Display_Displacements(STATE)

    % display smooth displacements
    disp_x = zeros(1,STATE.current_mesh.node_num);
    disp_y = zeros(1,STATE.current_mesh.node_num);
    disp_z = zeros(1,STATE.current_mesh.node_num);
    for n = 1:STATE.current_mesh.node_num,
        disp_x(n) = STATE.nodes(n).smooth_disp(1);
        disp_y(n) = STATE.nodes(n).smooth_disp(2);
        disp_z(n) = STATE.nodes(n).smooth_disp(3);
    end
    delete(STATE.smooth_disp_handle);
    vis = 'off';
    if (get(STATE.smooth_disp_checkbox,'Value') == 1)
        vis = 'on';
    end
    STATE.smooth_disp_handle = quiver3(STATE.current_mesh.node_x, STATE.current_mesh.node_y, STATE.current_mesh.node_z, disp_x, disp_y, disp_z, 0, STATE.SMOOTH_DISP_COLOR, 'Visible', vis);
    
    % display point displacements
    for n = 1:STATE.current_mesh.node_num,
        disp_x(n) = STATE.nodes(n).point_disp(1);
        disp_y(n) = STATE.nodes(n).point_disp(2);
        disp_z(n) = STATE.nodes(n).point_disp(3);
    end
    delete(STATE.point_disp_handle);
    vis = 'off';
    if (get(STATE.point_disp_checkbox,'Value') == 1)
        vis = 'on';
    end
    STATE.point_disp_handle = quiver3(STATE.current_mesh.node_x, STATE.current_mesh.node_y, STATE.current_mesh.node_z, disp_x, disp_y, disp_z, 0, STATE.POINT_DISP_COLOR, 'Visible', vis);

end % Display_Displacements


%
% Function: Display_State
%
function STATE = Display_Surface(STATE)

    % update the state label
    set(STATE.state_label, 'String', sprintf('iteration: %d', STATE.iteration));
    
    % display the surface
    delete(STATE.surface_handle);
    vis = 'off';
    if (get(STATE.surface_checkbox,'Value') == 1)
        vis = 'on';
    end
    axes(STATE.main_axes);
    STATE.surface_handle = Display_Mesh(STATE.current_mesh, STATE.SURFACE_COLOR, get(STATE.surface_transparency_slider,'Value'), STATE.EDGE_COLOR, vis);
    
end % Display_State


%
% Function: Display_Mesh
%
function surface_handle = Display_Mesh(mesh, face_color, face_alpha, edge_color, vis)

    surface_handle = trisurf([mesh.tri_n1' mesh.tri_n2' mesh.tri_n3'], mesh.node_x, mesh.node_y, mesh.node_z, 'FaceColor', face_color, 'FaceAlpha', face_alpha, 'EdgeColor', edge_color, 'FaceLighting', 'phong', 'Visible', vis);
    
end % Display_Mesh


%
% Function: Visibility_CB
%
function Visibility_CB(obj, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % get the handle
    if (obj == STATE.surface_checkbox)
        handle = STATE.surface_handle;
    end
    if (obj == STATE.smooth_disp_checkbox)
        handle = STATE.smooth_disp_handle;
    end
    if (obj == STATE.point_disp_checkbox)
        handle = STATE.point_disp_handle;
    end
    
    % set the visibility
    vis = 'off';
    if (get(obj,'Value') == 1)
        vis = 'on';
    end
    set(handle, 'Visible', vis);
    
end % Visibility_CB


%
% Function: Surface_Transparency_Slider_CB
%
function Surface_Transparency_Slider_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % set the surface transparency
    set(STATE.surface_handle, 'FaceAlpha', get(STATE.surface_transparency_slider,'Value'));
    
end % Surface_Transparency_Slider_CB


%
% Function: Update_Displacements_CB
%
function Update_Displacements_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % get nor_R_1
    nor_R_1 = Read_Nor_R_1(STATE);
    if isempty(nor_R_1)
        return;
    end
    
    % get tan_R_1
    tan_R_1 = Read_Tan_R_1(STATE);
    if isempty(tan_R_1)
        return;
    end
    
    % compute node normals
    STATE.nodes = Surface_Normals_at_Nodes(STATE.nodes, STATE.current_mesh);
    
    % compute displacements
    STATE = Compute_Displacements(STATE, nor_R_1, tan_R_1);
    
    % display the surface
    STATE = Display_Displacements(STATE);
    
    % save STATE
    set(fig, 'UserData', STATE);
    
end % Update_Displacements_CB


%
% Function: Surface_Normals_at_Nodes
%
function nodes = Surface_Normals_at_Nodes(nodes, mesh)
   
    % process node by node
    for n = 1:mesh.node_num,
        
        % compute the center
        ind = nodes(n).neighbors(:,1)';
        mean_x = mean(mesh.node_x(ind));
        mean_y = mean(mesh.node_y(ind));
        mean_z = mean(mesh.node_z(ind));

        % node normal
        X = [ mesh.node_x(ind)'-mean_x mesh.node_y(ind)'-mean_y mesh.node_z(ind)'-mean_z ];
        
        % compute the eigenvalues of X'*X
        [V,D] = eig(X'*X);
        
        % the first eigenvalue should be the smallest - check
        m = min([ D(1,1) D(2,2) D(3,3) ]);
        if (m ~= D(1,1))
            error('Adjust_Mesh::Surface_Normals_at_Nodes: the first eignevalue is not the smallest, which is assumed to be the case!');
        end

        % node normal
        nodes(n).normal = V(:,1);
        
        % set the normal orientation
        n1 = mesh.tri_n1(nodes(n).tri_ind);
        n2 = mesh.tri_n2(nodes(n).tri_ind);
        n3 = mesh.tri_n3(nodes(n).tri_ind);
        e12 = [ mesh.node_x(n2)-mesh.node_x(n1) mesh.node_y(n2)-mesh.node_y(n1) mesh.node_z(n2)-mesh.node_z(n1) ];
        e13 = [ mesh.node_x(n3)-mesh.node_x(n1) mesh.node_y(n3)-mesh.node_y(n1) mesh.node_z(n3)-mesh.node_z(n1) ];
        if (dot(nodes(n).normal,cross(e12,e13)) > 0)
            nodes(n).normal = -nodes(n).normal;
        end
   
    end % n loop

end % Surface_Normals_at_Nodes


%
% Function: Run_CB
%
function Run_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    
    % read the parameters
    nor_R_1 = Read_Nor_R_1(STATE);
    if isempty(nor_R_1)
        return;
    end
    nor_R_2 = Read_Nor_R_2(STATE);
    if isempty(nor_R_2)
        return;
    end
    tan_R_1 = Read_Tan_R_1(STATE);
    if isempty(tan_R_1)
        return;
    end
    tan_R_2 = Read_Tan_R_2(STATE);
    if isempty(tan_R_2)
        return;
    end
    smooth_1 = Read_Smooth_1(STATE);
    if isempty(smooth_1)
        return;
    end
    smooth_2 = Read_Smooth_2(STATE);
    if isempty(smooth_2)
        return;
    end
    points_1 = Read_Points_1(STATE);
    if isempty(points_1)
        return;
    end
    points_2 = Read_Points_2(STATE);
    if isempty(points_2)
        return;
    end
    iterations = Read_Iterations(STATE);
    if isempty(iterations)
        return;
    end

    
    % iterations
    for n = 1:iterations,
        
        k = (n-1)/(iterations-1);
        nor_R = nor_R_1*(1-k) + nor_R_2*k;
        tan_R = tan_R_1*(1-k) + tan_R_2*k;
        smoothness_weight = smooth_1*(1-k) + smooth_2*k;
        points_weight = points_1*(1-k) + points_2*k;

        fprintf('Adjust_Mesh: iteration %d (nor_R=%g   tan_R=%g   smoothness=%g   points=%g)\n', n, nor_R, tan_R, smoothness_weight, points_weight);
        STATE.current_mesh = Next(STATE.nodes, STATE.current_mesh, STATE.points, nor_R, tan_R, smoothness_weight, points_weight);
        
        
        % update the iteration
        STATE.iteration = STATE.iteration + 1;
    
        % display the surface
        STATE = Display_Surface(STATE);
        
        % display the displacements
        STATE = Display_Displacements(STATE);
        
        % force the rendering
        drawnow;

    end % m loop
    
    
    % save STATE
    set(fig, 'UserData', STATE);
    
end % Run_CB


%
% Function: Terminate_CB
%
function Terminate_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');

    % ask the user to confirm
    ButtonName = questdlg('This will terminate the application. To accept the mesh press "Accept". Still terminate?', 'Terminate?', 'Terminate', 'Cancel', 'Cancel');

    % if the user wants to stay
    if (strcmp(ButtonName, 'Cancel') == 1)
        return;
    end
    
    % empty output
    set(STATE.hidden_figure, 'UserData', []);
    
    % delete the figures
    Delete_Figures(STATE);
    
end % Terminate_CB


%
% Function: Accept_CB
%
function Accept_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % ask the user to confirm
    ButtonName = questdlg(sprintf('Accept the %s mesh?', STATE.name), 'Accept?', 'Accept', 'Cancel', 'Cancel');

    % if the user wants to stay
    if (strcmp(ButtonName, 'Cancel') == 1)
        return;
    end

    % store the mesh for the output
    set(STATE.hidden_figure, 'UserData', STATE.current_mesh);

    % delete the figures
    Delete_Figures(STATE);

    
end % Accept_CB


%
% Function: Help_CB
%
function Help_CB(~, ~, fig)

    % get STATE
    STATE = get(fig, 'UserData');
    
    % show the help window
    set(STATE.help_fig, 'Visible', 'on');
    
end % Help_CB



%
% Function: Next
%
function mesh = Next(nodes, mesh, points, NOR_R, TAN_R, smoothness_weight, points_weight)

    % compute node normals
    nodes = Surface_Normals_at_Nodes(nodes, mesh);

    % node displacements from the smoothness
    nodes = Node_Displacements_from_Smoothness(nodes, mesh);
    
    % node displacements from points
    nodes = Node_Displacements_from_Points(nodes, mesh, points, NOR_R, TAN_R);
    
    % update the mesh
    for n = 1:mesh.node_num,
        mesh.node_x(n) = mesh.node_x(n) + smoothness_weight*nodes(n).smooth_disp(1) + points_weight*nodes(n).point_disp(1);
        mesh.node_y(n) = mesh.node_y(n) + smoothness_weight*nodes(n).smooth_disp(2) + points_weight*nodes(n).point_disp(2);
        mesh.node_z(n) = mesh.node_z(n) + smoothness_weight*nodes(n).smooth_disp(3) + points_weight*nodes(n).point_disp(3);
    end

end % Next


%
% Function: Compute_Displacements
%
function STATE = Compute_Displacements(STATE, nor_R, tan_R)

    % node displacements from the smoothness
    STATE.nodes = Node_Displacements_from_Smoothness(STATE.nodes, STATE.current_mesh);
    
    % node displacements from points
    STATE.nodes = Node_Displacements_from_Points(STATE.nodes, STATE.current_mesh, STATE.points, nor_R, tan_R);

end % Compute_Displacements


%
% Function: Node_Displacements_from_Smoothness
% -----------------------------------------------------------------------
%
function nodes = Node_Displacements_from_Smoothness(nodes, mesh)

    % loop over nodes
    for n = 1:mesh.node_num,
        
        % compute the center
        ind = nodes(n).neighbors(:,1)';
        ind2 = [ ind(2:end) ind(1) ];
        dx = mesh.node_x(ind2) - mesh.node_x(ind);
        dy = mesh.node_y(ind2) - mesh.node_y(ind);
        dz = mesh.node_z(ind2) - mesh.node_z(ind);
        d = sqrt( dx.^2 + dy.^2 + dz.^2 );
        mx = .5 * ( mesh.node_x(ind2) + mesh.node_x(ind) );
        my = .5 * ( mesh.node_y(ind2) + mesh.node_y(ind) );
        mz = .5 * ( mesh.node_z(ind2) + mesh.node_z(ind) );
        center_x = sum(mx.*d) / sum(d);
        center_y = sum(my.*d) / sum(d);
        center_z = sum(mz.*d) / sum(d);
        
        % compute the average normal distance to triangle sides
        d21x = mesh.node_x(ind2) - mesh.node_x(ind);
        d21y = mesh.node_y(ind2) - mesh.node_y(ind);
        d21z = mesh.node_z(ind2) - mesh.node_z(ind);
        d10x = mesh.node_x(ind) - center_x;
        d10y = mesh.node_y(ind) - center_y;
        d10z = mesh.node_z(ind) - center_z;
        d = sqrt( (d10x.^2+d10y.^2+d10z.^2).*(d21x.^2+d21y.^2+d21z.^2) - (d10x.*d21x+d10y.*d21y+d10z.*d21z).^2 ) ./ sqrt(d21x.^2+d21y.^2+d21z.^2);
        d_avg = mean(d);  
        
        % node normal
        node_normal = nodes(n).normal;
        
        
        % compute the intersection of second triangle planes with the node
        % normal
        cp = [ 0 0 0 ]';
        l0 = [ center_x center_y center_z ]';
        M = size(nodes(n).neighbors,1);
        for m = 1:M,
            mn = m+1;
            if (m == M)
                mn = 1;
            end
            n1 = nodes(n).neighbors(m,1);
            n2 = nodes(n).neighbors(mn,1);
            n3 = nodes(n).neighbors(m,2);
            p1 = [ mesh.node_x(n1) mesh.node_y(n1) mesh.node_z(n1) ]';
            p2 = [ mesh.node_x(n2) mesh.node_y(n2) mesh.node_z(n2) ]';
            p3 = [ mesh.node_x(n3) mesh.node_y(n3) mesh.node_z(n3) ]';
            p = Line_Plane_Intersection(l0, node_normal, p1, p2, p3);
            cp = cp + p;
            % tmp
            % plot3([p1(1) p(1)], [p1(2) p(2)], [p1(3) p(3)], 'm');
            % plot3([p2(1) p(1)], [p2(2) p(2)], [p2(3) p(3)], 'm');
            % plot3([p3(1) p(1)], [p3(2) p(2)], [p3(3) p(3)], 'm');
        end % m loop
        cp = cp / M;
  
        % compute the new point
        H = norm(cp-l0);
        h = d_avg * tan(atan(H/d_avg)/3);
        np = l0 + h * (cp-l0)/H;
        
        % store the displacement due to smoothness
        nodes(n).smooth_disp = [ np(1)-mesh.node_x(n) np(2)-mesh.node_y(n) np(3)-mesh.node_z(n) ];
        
    end % n loop

end % Node_Displacements_from_Smoothness


%
% Function: Line_Plane_Intersection
% -----------------------------------------------------------------------
% This function computes the intersection of a line (defined by point l0
% and direction l) and a plane (defined by three points p1, p2, p3). The
% function returns the point of intersection p. If there is no intersection
% p = [].
%
function p = Line_Plane_Intersection(l0, l, p1, p2, p3)

    % plane normal
    n = cross(p2-p1,p3-p1);
    nn = norm(n);
    if (nn <= eps)
        error('Adjust_Mesh::Line_Plane_Intersection: The three points do not define a plane!');
    end
    n = n / nn;
    
    num = dot(l,n);
    if (abs(num) <= eps)
        p = [];
    else
        d = dot((p1-l0),n) / num;
        p = l0 + d*l;
    end

end % Line_Plane_Intersection


%
% Function: Node_Neighborhoods
% -----------------------------------------------------------------------
% It computes the neighbood for each node. Each node will have a two-column
% matrix "neighbors". The first column contains the index of the
% neighboring node, while the second column contines the index of the "star
% tip second neighbor". Typically a node has six neighbors and therefore
% typically "neighbors" is a 6x2 matrix, representing neighbor info around
% the node.
%
function nodes = Node_Neighborhoods(mesh)

    % preallocate and initialize the nodes array
    nodes(mesh.node_num).count = [];
    for n = 1:mesh.node_num,
        nodes(n).count = 0;
        nodes(n).point_disp = [ 0 0 0 ];
    end
    
    % loop over nodes
    fprintf('Adjust_Mesh: Computing node neighbors: ');
    previous_percent = -1;
    for n = 1:mesh.node_num,
        
        
        % print progress
        percent = 10*floor(10*n/mesh.node_num);
        if (percent ~= previous_percent)
            if (percent == 10)
                fprintf('\b\b');
            end
            if (percent > 10)
                fprintf('\b\b\b');
            end
            fprintf('%d%%', percent);
            previous_percent = percent;
        end
        
        % find a triangle that contains node n
        ind = find(mesh.tri_n1 == n, 1);
        if isempty(ind)
            ind = find(mesh.tri_n2 == n, 1);
            if isempty(ind)
                ind = find(mesh.tri_n3 == n, 1);
                if isempty(ind)
                    error('Adjust_Mesh::Node_Neighborhoods: Invalid situation!');
                else
                    n1 = mesh.tri_n1(ind);
                    n2 = mesh.tri_n2(ind);
                end
            else
                n1 = mesh.tri_n3(ind);
                n2 = mesh.tri_n1(ind);
            end
        else
            n1 = mesh.tri_n2(ind);
            n2 = mesh.tri_n3(ind);
        end
        
        % store the triangle index
        nodes(n).tri_ind = ind;
        
        % find the forth point
        [n3, ~] = Find_Triangle(mesh, n1, n2, n);
        
        % store the neighbors
        nodes(n).neighbors = [ n1 n3 ];
        
        % go around the node
        ind = 1;
        while(ind)
            
            % find the next triangle
            [n4, ~] = Find_Triangle(mesh, n, n2, n1);
            
            % rename the points
            n1 = n2;
            n2 = n4;
            
            % find the forth point
            [n3, ~] = Find_Triangle(mesh, n1, n2, n);
            
            % store the neighbors
            nodes(n).neighbors = [ nodes(n).neighbors ; n1 n3 ];
            
            % check if the circle is complete
            if (n4 == nodes(n).neighbors(1,1))
                ind = 0;
            end
            
        end % while loop
        
    end % n loop
    fprintf('\n');
    
end % Node_Neighborhoods


%
% Function: Find_Triangle
% -----------------------------------------------------------------------
% This function finds the triangle, i.e. the node n3 that together with
% nodes n1 and n2 forms a triangle, but it does not contain node n.
%
function [n3, ind] = Find_Triangle(mesh, n1, n2, n)

    % n1 - n2 - n
    ind = find( (mesh.tri_n1 == n1) & (mesh.tri_n2 == n2) & (mesh.tri_n3 ~= n) );
    if (isempty(ind) == 0)
        if (length(ind) > 1)
            error('Adjust_Mesh::Find_Triangle: Invalid situation (n1-n2-n)!');
        end
        n3 = mesh.tri_n3(ind);
        return;
    end
    
    % n1 - n - n2
    ind = find( (mesh.tri_n1 == n1) & (mesh.tri_n2 ~= n) & (mesh.tri_n3 == n2) );
    if (isempty(ind) == 0)
        if (length(ind) > 1)
            error('Adjust_Mesh::Find_Triangle: Invalid situation (n1-n-n2)!');
        end
        n3 = mesh.tri_n2(ind);
        return;
    end
    
    % n2 - n1 - n
    ind = find( (mesh.tri_n1 == n2) & (mesh.tri_n2 == n1) & (mesh.tri_n3 ~= n) );
    if (isempty(ind) == 0)
        if (length(ind) > 1)
            error('Adjust_Mesh::Find_Triangle: Invalid situation (n2-n1-n)!');
        end
        n3 = mesh.tri_n3(ind);
        return;
    end
    
    % n2 - n - n1
    ind = find( (mesh.tri_n1 == n2) & (mesh.tri_n2 ~= n) & (mesh.tri_n3 == n1) );
    if (isempty(ind) == 0)
        if (length(ind) > 1)
            error('Adjust_Mesh::Find_Triangle: Invalid situation (n2-n-n1)!');
        end
        n3 = mesh.tri_n2(ind);
        return;
    end
    
    % n - n1 - n2
    ind = find( (mesh.tri_n1 ~= n) & (mesh.tri_n2 == n1) & (mesh.tri_n3 == n2) );
    if (isempty(ind) == 0)
        if (length(ind) > 1)
            error('Adjust_Mesh::Find_Triangle: Invalid situation (n-n1-n2)!');
        end
        n3 = mesh.tri_n1(ind);
        return;
    end
    
    % n - n2 - n1
    ind = find( (mesh.tri_n1 ~= n) & (mesh.tri_n2 == n2) & (mesh.tri_n3 == n1) );
    if (isempty(ind) == 0)
        if (length(ind) > 1)
            error('Adjust_Mesh::Find_Triangle: Invalid situation (n-n2-n1)!');
        end
        n3 = mesh.tri_n1(ind);
        return;
    end
    
    error('Adjust_Mesh::Find_Triangle: Invalid situation (no triangle match)!');

end % Find_Triangle


%
% Function: Node_Displacements_from_Points
% -----------------------------------------------------------------------
%
function nodes = Node_Displacements_from_Points(nodes, mesh, points, NORMAL_R, TANGENT_R)
    
    % process node by node
    for n = 1:mesh.node_num,
        
        % node normal
        node_normal = nodes(n).normal;
        
        % compute the normal and tangent components
        dx = points(:,1) - mesh.node_x(n);
        dy = points(:,2) - mesh.node_y(n);
        dz = points(:,3) - mesh.node_z(n);
        nor = dx*node_normal(1) + dy*node_normal(2) + dz*node_normal(3);
        tan = sqrt( (dx-nor*node_normal(1)).^2 + (dy-nor*node_normal(2)).^2 + (dz-nor*node_normal(3)).^2 );
        
        eds = (nor/NORMAL_R).^2 + (tan/TANGENT_R).^2;
        ind = find(eds <= 1);  
        
        nodes(n).point_disp = [0 0 0];
        nodes(n).point_disp_nc = 0; % point displacement normal component
        c = 0;
        nodes(n).point_disp_computed = 0;
        if (isempty(ind) == 0)
            for m = 1:length(ind),
                im = ind(m);
                w = 1 - eds(im);
                nm = nor(im);
                nodes(n).point_disp = nodes(n).point_disp + w * nm * node_normal';
                nodes(n).point_disp_nc = nodes(n).point_disp_nc + w * nm;
                c = c + w;
            end % m loop
            nodes(n).point_disp = nodes(n).point_disp / c;
            nodes(n).point_disp_nc = nodes(n).point_disp_nc / c;
            nodes(n).point_disp_computed = 1;
        end
        
    end % n loop
    
    
    %
    % Spread the point displacements using "Laplacian smoothing"
    %
    
    % point_disp_normal
    point_disp_nc = zeros(1, mesh.node_num);
    for n = 1:mesh.node_num,
        if (nodes(n).point_disp_computed == 1)
            point_disp_nc(n) = nodes(n).point_disp_nc;
        end
    end
    
    
    M = 100; % number of iterations (practice shows that this is a reasonable value)
    for m = 1:M,
        
        % tmp storage
        nc = zeros(1, mesh.node_num);
        
        % process node by node
        for n = 1:mesh.node_num,
            
            % if the node point displacement needs to be adjusted
            if (nodes(n).point_disp_computed == 0)
            
                % node neighbors
                ind = nodes(n).neighbors(:,1)';
                
                % the average of neighbors
                nc(n) = mean(point_disp_nc(ind));
                      
            end % if ...
            
        end % n loop
        
        % copy the results
        for n = 1:mesh.node_num,
            if (nodes(n).point_disp_computed == 0)
                point_disp_nc(n) = nc(n);
            end
        end % n loop
        
    end % m loop
    
    
    % get the final point displacements
    for n = 1:mesh.node_num,
        
        nodes(n).point_disp = point_disp_nc(n) * nodes(n).normal';
        
    end % n loop

end % Node_Displacements_from_Points
