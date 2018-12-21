function [] = buildCalibrationModels(input,varargin)
%buildCalibrationModels - Generate calibration models for converting G-values to OD equivalents
%
% buildCalibrationModels([])
% buildCalibrationModels(file)
% buildCalibrationModels(..., plot)
%
% Calling the function without inputs will open a window for selecting an
% Excel-file containing calibration data. The function can also be called
% with a string indicating the path and name of a file.
%
% Optional input
% plot      'on' (default) or 'off'. Displays a plot of each calibration model.
% 
% The output is a mat-file which will be used when calling the 'convertGrowthProfilerGvalues' function.

%% HANLDE INPUT
if isempty(input)
    [filename, pathname] = uigetfile( {'*.xls;*.xlsx', 'EXCEL Files (*.xls, *.xlsx)';'*.*','All Files (*.*)'},'Open Excel data file');
    file = fullfile(pathname, filename);
end


if nargin == 1
    plotFlag = 'on';
else
    plotFlag = varargin{1};
end
    
%% GENERATE MODELS
calModels = cell(2,3);
plates = {'24WRO','24GSQ','96WSQ'};

for scanner = 1:2
    for  p = 1:3
        [~, ~, raw] = xlsread(file,plates{p});
        
        if size(raw,1) > 9
            
            switch scanner
                case 1, calData = raw(9:end,1:3);
                case 2, calData = raw(9:end,5:7);
            end
            
            calData = reshape([calData{:}],size(calData));
            
            switch plates{p}
                case '96WSQ', numKnots = 8;
                case {'24WRO','24GSQ'}, numKnots = 10;
            end
            
            model = slmengine(calData(:,2),calData(:,1),'plot',plotFlag,'concaveup','on','incr','on','minvalue',0,'knots',numKnots,'verbosity',1);
            calModels{scanner,p} = model;
        end
    end    
end

[FileName,PathName] = uiputfile({'*.mat', 'MAT File (*.mat)';'*.*','All Files (*.*)'},'Save as mat-file');
if isequal(filename,0) || isequal(pathname,0)
    error('No valid input given. Results not saved.')
else
    save(fullfile(PathName,FileName),'calModels')    
end