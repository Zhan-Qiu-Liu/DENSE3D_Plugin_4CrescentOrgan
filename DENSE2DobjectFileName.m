% DENSE2DobjectFileName  - Read the File Name of DENSE2Dobject
%
%
%	Set current directory in Matlab as the location where Inputs are
%
%
% USAGE:
%   files = DENSE2DobjectFileName()
%
% INPUTS:
%	none
%
% OUTPUTS:
%	Path of files 
%
%
% See also:
%   
%
% Last Modified: 1:06 PM Thursday, October 15, 2015
% Modified By: Zhanqiu Liu (lafeir.lew@gmail.com)

% function [files, tmp] = DENSE2DobjectFileName()
function files = DENSE2DobjectFileName(directory,type)
	%set current Dir.
	if exist('directory', 'var')
		tmp = dir(directory); path = directory;
	else
		path = uigetdir(pwd, 'Select the directory storing DENSE2D Outputs(*.mat)');
		tmp = dir(path);
	end
	
	tmp = {tmp.name};

	%find  structure of garbage 
	temp2 = strncmpi(tmp, 'auto.', 5);	
	jj = 0;
	for ii = 1:size(tmp,2)
		if temp2(ii)
			jj = jj + 1;
			files(jj) = fullfile(path,tmp(ii));%load names
		end
	end

	if exist('type', 'var')
		temp2 = strncmpi({'SA','LA'}, type, 5);
		if sum(temp2) == 0
			error('Incorrect 2nd input parameter--"Type":it should be "SA" or "LA"!');
			return
		end
		for ii = jj:-1:1
			load(char(files(ii)));
			if ~strcmpi(ROIInfo.ROIType,type)
				files(ii) = [];
			end
		end
	end

end


