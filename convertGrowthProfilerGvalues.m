function [] = convertGrowthProfilerGvalues(calModels,varargin)
%convertGrowthProfilerGvalues - Converts G-values to OD equivalents using
%calibration data
%
% convertGrowthProfilerGvalues(calModels)
% convertGrowthProfilerGvalues(calModels,datafiles)
%
% Required input:
% calModels         Matlab cell array containing calibration models
%                   generated using the 'buildCalibrationModels' function
%
% Calling the function with a single input will open a window for selecting an
% Excel-file containing G-values.
%
% Optional input:
% datafiles         Cell aray of strings indicating full path and filename
%                   of Excel-files with G-values to be batch processes.
%
% The output is an Excel-file in the specified directory or one file in each
% directory in case of batch processing with quanitifed OD equivalents.
%

%% HANLDE INPUT
if nargin == 1
    [filename, pathname] = uigetfile( {'*.xls;*.xlsx', 'EXCEL Files (*.xls, *.xlsx)';'*.*','All Files (*.*)'},'Open Excel data file','MultiSelect','on');        
    if iscell(filename)
        datafile = strcat({pathname},filename');
    else
        datafile = {fullfile(pathname, filename)};
    end
    
    if isequal(filename,0) || isequal(pathname,0)
        error('No valid file selected. Program ended.')
    end
    
else
    if iscell(varargin{1})
        datafile = varargin{1};
    else
        error('Second input argument should be a cell array of strings.')
    end
end

%% CONVERT G-VALUES TO OD EQUIVALENTS
warning('off','MATLAB:xlswrite:AddSheet');
numFiles = length(datafile);
plates = {'24WRO','24GSQ','96WSQ'};

for f = 1:numFiles    
    [pathstr,filename,ext] = fileparts(datafile{f});
    [~,~,expInfo] = xlsread(datafile{f},'Info');
    projectDate = expInfo{2,2};
    projectStart = datestr(expInfo{3,2},'HH:MM:SS');
    expInfo{5,4} = 'Number of time points';    
    
    outputfile = fullfile(pathstr, sprintf('%s - OD %s%s', filename,datestr(date, 29),ext));    
    xlswrite(outputfile, expInfo(1:5,:), 'Info');    
    
    % Quantify the data from each tray    
    fprintf('Processing File %d of %d: %s \n',f,numFiles,datafile{f});
    trayInfo = expInfo(6:17,1:end);
    for i = 1:12
        fprintf('Evaluating Tray %d... ',i);
        tray = trayInfo{i,1};
        scannerUsed = trayInfo{i,2};
        plateUsed = trayInfo{i,3};
        timeInit = trayInfo{i,4}/60;
        
        switch scannerUsed
            case {1,2}
                isModel = ismember(plates,plateUsed);
                
                if sum(isModel) == 0 && ~strcmp(plateUsed,'Empty')
                    error('Plate abreviation not valid. Please specify 96WSQ, 24WRO or 24GSQ.');
                elseif strcmp(plateUsed,'Empty')
                    fprintf('Tray empty. Quantification skipped. \n');
                    continue
                elseif isempty(calModels{scannerUsed,isModel})
                    warning('There is no calibration data available for the specified plate type. Conversion to OD equivalents skipped.');
                    continue
                else
                    fprintf('Quantifying data... ');
                    model = calModels{scannerUsed, isModel};
                    [~,~,rawdata] = xlsread(datafile{f},tray);                    
                    isBlank = ismember(rawdata(7,:),{'Empty','Blank','Error'});
                    rawdata(:,isBlank) = [];                    
                    time = rawdata(10:end,1); time = reshape([time{:}],size(time))./60;
                    gvalues = rawdata(10:end,2:end); gvalues = reshape([gvalues{:}],size(gvalues));
                    
                    gvalues(isnan(time),:) = [];
                    time(isnan(time)) = [];
                    
                    [time,ODcalc,numT] = quantify(time,gvalues,timeInit,model);                                       
                                        
                    str1 = rawdata(7,2:end)'; %Strains
                    str2 = rawdata(4,2:end)'; %Conditions
                    ltmp = rawdata(5,2:end)'; %Levels
                    
                    str3 = cell(length(ltmp),1);
                    for j = 1:length(ltmp)
                        str3{j,1} = num2str(ltmp{j},'%0.3f'); 
                    end
                    
                    %str4 = {[projectDate '.' projectStart]}; %Time stamps
                    str5 = {sprintf('T%02d',i)}; %Tray number
                    str6 = rawdata(3,2:end)'; %Wells
                    
                    S = strcat({'['},str1,{'].['},str2,{'].['},str3,{'].['},projectDate,{'].['},projectStart,{'].['},str5,{'].['},str6,{']'});
                                        
                    M1 = rawdata(1:8,:);                   
                    M1{1,1} = [tray ' - OD-values'];
                    M1(2,:) = [{'Tray:'},repmat(str5,1,size(rawdata,2)-1)];                    
                    M1(8,:) = [{'Experiment ID:'},S'];                   
                    M1(9,:) = cell(1,size(M1,2));
                    M2 = [{'Time'},repmat({'OD-value'},1,size(rawdata,2)-1)];                    
                    M3 = num2cell([time ODcalc]);
                    M = [M1;M2;M3];
                    
                    xlswrite(outputfile, M, tray);
                    fprintf('Done.\n')
                end
                
            otherwise
                error('The Growth Profiler only has two scanners. Please specify Scanner 1 or 2.');
        end
        trayInfo{i,4} = numT;
    end
    xlswrite(outputfile, trayInfo, 'Info','A6');
    fprintf('Processing complete. Results saved to: %s\n\n',outputfile)
end

%% AUXILIARY FUNCTIONS
% -------------------------------------------------------------------------
function [time,ODcalc,n] = quantify(time,gvalues,timeInit,model)

if timeInit ~= 0    
    t = find(time == timeInit);
    if isempty(t)
        error('Indicated time point not found. Program ended.')
    else
        time = (time(t:end) - timeInit);
        gvalues = gvalues(t:end,:);
    end
end
n = length(time);
ODcalc = slmeval(gvalues,model);

end

end

% % Set default values
% datafile = [];
% s = what('Growth Profiling Toolbox 1.0');
% calibrationfile = [s.path filesep 'calibration data' filesep 'GrowthProfiler_v1_CalibrationData.xlsx'];
% plotFlag = 'off';
% 
% varargin
% if ~isempty(varargin{1})
%     if mod(length(varargin),2)
%         error('Number of input variables is incorrect. Please provide at least one pair of inputs.');
%     else
%         for i = 1:2:length(varargin)
%             switch varargin{i}
%                 case 'datafile'
%                     datafile = varargin{i+1};
%                 case 'calibrationfile'
%                     calibrationfile = varargin{i+1};
%                 case 'plot'
%                     plotFlag = varargin{i+1};
%             end
%         end
%     end
% end
% 
% 
% if isempty(datafile)
%     [filename, pathname] = uigetfile( {'*.xls;*.xlsx', 'EXCEL Files (*.xls, *.xlsx)';'*.*','All Files (*.*)'},'Open Excel data file');
%     datafile = fullfile(pathname, filename);
% end
% 
% %% READ RAW DATA
% [~,~,expInfo] = xlsread(datafile,'Info');
% expDate = expInfo{2,2};
% trayInfo = expInfo(5:16,1:end);
% 
% outputfile = fullfile(pathname, [expDate ' Quantification results ' date '.xlsx']);
% xlswrite(outputfile, expInfo, 'Info');
% 
% %% GENERATE REQUIRED CALIBRATION MODELS
% 
% 
% %% QUANTIFY
% calModels = cell(6,2); calModels(:,1) = repmat({'E'},6,1);
% cc = 2;
% for i = 1:12
%     %trayData = trayInfo{i,1:4};
%     
%     if isempty(calModels)
%         % Generate calibration data
%         calModels{1,1} = [trayInfo{i,3} '-' trayInfo{i,2}];
%         model = genCalModel(calibrationfile,trayInfo(i,2:3),plotFlag);
%         calModels{1,2} = model;
%     else
%         % Check if calibration data already exist
%         [trayInfo{i,3} '-' num2str(trayInfo{i,2})]
%         id = strmatch([trayInfo{i,3} '-' num2str(trayInfo{i,2})], calModels(:,1), 'exact')
%         %id = find(ismember(calModels(:,1),[trayInfo{i,3} '-' trayInfo{i,2}]));
%         if isempty(id) &&
%             calModels{cc,1} = [trayInfo{i,3} '-' trayInfo{i,2}];
%             model = genCalModel(calibrationfile,trayInfo(i,2:3),plotFlag);
%             calModels{cc,2} = model;
%             cc = cc+1;
%         else
%             model = calModels{id,2};
%         end
%     end
%     
%     M = quantify(datafile,trayInfo(i,[1,4]),model);
%     xlswrite(outputfile, M, trayInfo{i,1});
% end
% 
% %% AUXILIARY FUNCTIONS
% % -------------------------------------------------------------------------
% function model = genCalModel(calfile,D,plotFlag)
% scannerUsed = D{1};
% plateUsed = D{2};
% 
% switch plateUsed
%     case {'96WSQ','24WRO','24GSQ'}
%         [~, ~, raw] = xlsread(calfile,plateUsed);
%         
%         switch scannerUsed
%             case 1, calData = raw(9:end,1:3);
%             case 2, calData = raw(9:end,5:7);
%             otherwise, error('The Growth Profiler only has two scanners. Please specify Scanner 1 or 2.');
%         end
%         
%         calData = reshape([calData{:}],size(calData));
%     otherwise, error('Calibration data for specified plate unavailable. Please specify 96WSQ, 24WRO or 24GSQ.');
% end
% 
% switch plateUsed
%     case '96WSQ', numKnots = 8;
%     case {'24WRO','24GSQ'}, numKnots = 10;
% end
% 
% model = slmengine(calData(:,2),calData(:,1),'plot',plotFlag,'concaveup','on','incr','on','minvalue',0,'knots',numKnots,'verbosity',1);
% 
% end
% % -------------------------------------------------------------------------
% function M = quantify(datafile,D,model)
% tray = D{1};
% [~,~,rawdata] = xlsread(datafile,tray);
% time = rawdata(10:end,1); time = reshape([time{:}],size(time))./60;
% gvalues = rawdata(10:end,2:end); gvalues = reshape([gvalues{:}],size(gvalues));
% 
% timeInit = D{2};
% if timeInit ~= 0
%     t = find(time == timeInit);
%     if isempty(t)
%         error('Indicated time point not found. Program ended.')
%     else
%         time = (time - timeInit);
%         gvalues = gvalues(t:end,:);
%     end
% end
% 
% ODcalc = slmeval(gvalues,model);
% 
% M1 = rawdata(1:9,:);
% M2 = num2cell([time ODcalc]);
% M = [M1;M2];
% end
% 
% 
% end
% 
