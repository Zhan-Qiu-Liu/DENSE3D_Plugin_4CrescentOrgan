classdef animatedgif < hgsetget
% Copyright (c) 2016 DENSEanalysis Contributors
% The codes below are modified from the function "switchstate" in "DENSEanalysis.m".

    % animatedgif - Class for creating animated GIFs from images
    %
    %   Similar to the MATLAB builtin class for generating movies. This
    %   class allows the user to dynamically create an animated GIF frame
    %   by frame.
    %
    %   Typical usage would be to create an animatedgif object and then use
    %   the addFrame method to add all the frames to the object:
    %
    %       A = animatedgif('test.gif');
    %       A.addFrame(rand(100));
    %
    %   Properties that control the display of the animated GIF are
    %   provided as properties of the class and can be set either in the
    %   constructor, through the SET method, or directly:
    %
    %       A = animatedgif('test.gif', 'LoopCount', 10);
    %       A.DelayTime = 0.1;
    %       set(A, 'Comment', 'This is an animated GIF')
    %
    % USAGE:
    %   A = animatedgif(filename);
    %
    % INPUTS:
    %   filename:   String, Filename to store the resulting animated GIF
    %
    % OUTPUTS:
    %   A:          Object, Object handle that can be used to add
    %               successive frames to the animated GIF.

    properties
        DelayTime           = 0.5;  % Time (in sec) between frames
        LoopCount           = inf;  % Number of times to loop through GIF
        Comment             = '';   % Comment to be supplied with GIF
        BackgroundColor     = [];   % Background color for optimization
        DisposalMethod      = 'doNotSpecify';   % Animation Method
        TransparentColor    = [];   % Index of desired transparent color
    end

    properties (Access = private)
        nFrames     = 0;    % Number of frames added to the GIF
        Filename    = '';   % Path to which the GIF is saved
    end

    properties (Dependent, Hidden)
        Inputs  % Cell Array for passing parameters to imwrite
    end

    methods
        %--- Get Functions for Dependent Properties ---%
        function res = get.Inputs(self)
            % Dynamically create property list to send to imwrite

            props       = get(self);

            % Remove all empty and private fields
            toremove    = structfun(@isempty, props);
            toremove(end-3:end) = 1;

            % Convert struct into cellular input
            vals    = struct2cell(props);
            names   = fieldnames(props);

            res = cat(1, names(~toremove)', vals(~toremove)');

            % If this isn't the first frame then we need to append
            if self.nFrames
                [~,ind] = ismember('LoopCount', names);
                res(:,ind) = {'WriteMode', 'Append'};
            end
        end

        function self = animatedgif(filename, varargin)
            % animatedgif - Constructor for animated GIFs
            %
            % USAGE:
            %   A = animatedgif(filename, params)
            %
            % INPUTS:
            %   filename:   String, Path to the GIF file to be created
            %   params:     Param/Value pairs, Parameter value inputs in
            %               the form ('Parameter', value) where the
            %               parameters are properties of this class such as
            %               "DelayTime" and "LoopCount". Type
            %               "properties('animatedgif')" in Matlab to see
            %               all available properties
            %
            % OUTPUTS:
            %   A:          animatedgif object, Instance of an animated gif
            %               object which can be used to append new frames
            %               to an image.

            if ~exist('filename', 'var')
                error(sprintf('%s:invalidInput', mfilename),...
                    'You must provide a filename');
            end

            if ~ischar(filename)
                % Assume that the first input is the image data (old style)
                data = filename;

                if numel(varargin) < 1 || ~ischar(varargin{1})
                    error(sprintf('%s:invalidInput', mfilename),...
                        'You must provide a filename');
                else
                    filename = varargin{1};
                    varargin(1) = [];
                end
            end

            self.Filename = filename;

            % Touch the file to go ahead and create it
            fclose(fopen(self.Filename, 'w'));

            if numel(varargin); set(self, varargin{:}); end

            % If the user supplied data, then go ahead and append it
            if exist('data', 'var')
                self.addFrame(data);
            end
        end

        function delete(self)
            % delete - Destructor for animatedgif object
            %
            %   Performs cleanup when the object instance is deleted
            %
            % USAGE:
            %   self.delete()

            % If there were no files written to the file then remove it
            if self.nFrames == 0
                delete(self.Filename);
            end
        end

        function addFrame(self, img)
            % addFrame - Adds a frame to the current animated GIF
            %
            %   This function can be called to append an image to the
            %   animated GIF object. This method call can be placed inside
            %   of a for loop with other functions such as GETFRAME to
            %   generate an animated GIF from an axes.
            %
            % USAGE:
            %   self.addFrame(mov)
            %
            % INPUTS:
            %   mov:    Struct or Matrix, The input can be either a
            %           structure similar to the output of GETFRAME or a
            %           matrix representing either a grayscale or RGB
            %           image.

            if isstruct(img)
                if numel(img) > 1
                    arrayfun(@(x)self.addFrame(x), img);
                    return;
                else
                    img = frame2im(img);
                end
            elseif ismatrix(img)
                img = repmat(img, [1 1 3]);
            elseif size(img, 3) ~= 3
                error(sprintf('%s:invaludInput', mfilename),...
                    'Input must be either gray, RGB, or GETFRAME output');
            end

            % Convert to an indexed image with a colormap
            [IND, map] = rgb2ind(img, 256);

            % Actually write the image to file using the specified inputs
            imwrite(IND, map, self.Filename, 'gif', self.Inputs{:});

            % Increment number of frames used so that we can ensure that
            % different parameters are passed on successive calls
            self.nFrames = self.nFrames + 1;
        end
    end
end
