function [] = collectExperimentalData(condition,directory,varargin)
% collectExperimentalData - Collects the experimental data from a set of
% files
%
% collectExperimentalData(condition,directory)
% collectExperimentalData(condition,directory,levels,strains)
%
% The function searches through all files in the specified directory for
% experimental data collected from a specific condition. The function 
% collects data from specified concentrations/levels.
% The output is an Excel-file with all replicate experiments for one 
% concentration/level of the condition in separate sheets.
% 
% Required inputs
% condition         String indicating the condition to search for, e.g.
%                   'Acetic acid'
% directory         String indicating the directory containing files with
%                   experimental data
% 
% Optional input
% levels            Numeric vector with concentrations/levels, e.g. [1,1.5,3,6]
%                   Leave empty to extract all concentrations/levels.
% strains           Cell array of strings with names of strains

%% HANLDE INPUT
levels = 'all';
strains = 'all';

if nargin > 2 && ~isempty(varargin{1})
    levels = reshape(varargin{1},1,length(varargin{1}));
end
if nargin > 3 && ~isempty(varargin{2})
    strains = reshape(varargin{2},1,length(varargin{2}));
end
    
%% SEARCH FILES
dirContents = struct2cell(dir(directory));
isdir = cell2mat(dirContents(4,:));
dirContents(:,isdir) = [];
numFiles = size(dirContents,2);
fileData = cell(numFiles,1);
allLevels = [];
allStrains = {};
allIDs = {};

numXLSfiles = 0;
maxNumTPs = 0;


for i = 1:numFiles
    [p,f,x] = fileparts([directory filesep dirContents{1,i}]);   
    
    if ~strcmp(x,'.xlsx')
        continue
    else   
        fprintf('Searching file: %s%s  ...',f,x);
        
        % Import experimental information
        datafile = fullfile([p filesep f x]);
        [~,~,expInfo] = xlsread(datafile,'Info');        
        projectDate = expInfo{2,2};
        projectStart = datestr(expInfo{3,2},'HH:MM:SS');
        
        % Find which trays contain the condition
        C = expInfo(6:17,:);
        isCond = ismember(C(:,5:end),condition);
        isCond = any(isCond,2);
        
        if sum(isCond) == 0
            fprintf(' Done (file does not contain specified condition) \n');
            continue
        else
            numXLSfiles = numXLSfiles + 1;
            selectedTrays = C(isCond,1);
            
            fileData{numXLSfiles}.maxNumTPs = max(cell2mat(C(isCond,4)));            
            if fileData{numXLSfiles}.maxNumTPs > maxNumTPs, maxNumTPs = fileData{numXLSfiles}.maxNumTPs; end
            
            fileData{numXLSfiles}.trayData = cell(length(selectedTrays),1);            
            for j = 1:length(selectedTrays)
                [~,~,dataMatrix] = xlsread(datafile,selectedTrays{j});
                dataMatrix{3,1} = 'Date:';                
                dataMatrix(3,2:end) = repmat({[projectDate ' ' projectStart]},1,size(dataMatrix,2)-1);                
                
                trayTimeVec = dataMatrix(:,1);
                trayDataMatrix = dataMatrix(:,2:end);
                
                isCondInTray = ismember(trayDataMatrix(4,:),condition);
                fileData{numXLSfiles}.trayData{j} = [trayTimeVec trayDataMatrix(:,isCondInTray)];
                
                allLevels = [allLevels, cell2mat(trayDataMatrix(5,isCondInTray))]; %#ok<AGROW>
                allStrains = [allStrains,trayDataMatrix(7,isCondInTray)]; %#ok<AGROW> 
                allIDs = [allIDs; trayDataMatrix(8,isCondInTray)']; %#ok<AGROW> 
            end
        end
    end
    fprintf(' Done\n');
end

levelList = unique(allLevels);
strainList = unique(allStrains);
idList = unique(allIDs);

%% DETERMINE THE SIZE OF THE DATA
if strcmp(levels,'all')
    levels = levelList;
end

if strcmp(strains,'all')
    strains = strainList;
end

numLvls = length(levels);
numExps = zeros(numLvls,1);
expDates = cell(1,numXLSfiles);
expList = cell(size(idList));

for i = 1:numXLSfiles        
    for j = 1:length(fileData{i}.trayData)                   
        expDates(i) = fileData{i}.trayData{j}(3,2);
        trayDataMatrix = fileData{i}.trayData{j}(:,2:end);
        trayDataMatrix(5,:)
        trayDataMatrix(7,:)
        for k = 1:numLvls                        
            numExps(k) = numExps(k) + sum(cell2mat(trayDataMatrix(5,:)) == levels(k) & ismember(trayDataMatrix(7,:),strains));
        end
    end
end

%% EXTRACT THE DATA FOR EACH LEVEL
output = cell(1,numLvls); 
fprintf('\nProcessing data...');
expCnt = 0;
colCnt = num2cell(zeros(1,numLvls));

for i = 1:numLvls
    output{i} = cell(maxNumTPs+10,numExps(i));
end

for j = 1:length(strains)
    currStrain = strains{j};
    %colCnt = 0;    
    
    for i = 1:numLvls
        currLvl = levels(i);
        %output{i} = cell(maxNumTPs+10,numExps(i));           
        
        for k = 1:numXLSfiles            
            for p = 1:length(fileData{k}.trayData)
                trayTimeVec = fileData{k}.trayData{p}(10:end,1);
                trayDataMatrix = fileData{k}.trayData{p}(:,2:end);
                
                isStrain = ismember(trayDataMatrix(7,:),currStrain);
                isLevel = cell2mat(trayDataMatrix(5,:)) == currLvl;
                
                if sum(isStrain & isLevel) == 0
                    continue                
                else                    
                    strainDataVec = trayDataMatrix(:,isLevel & isStrain);                    
                    
                    [m,n] = size(strainDataVec);
                    for cc = 1:n
                        strainTimeVec = [strainDataVec(1:9,cc); trayTimeVec];
                        colCnt{i} = colCnt{i} + 1;                        
                        output{i}(1:length(strainTimeVec),colCnt{i}) = strainTimeVec;
                        colCnt{i} = colCnt{i} + 1;
                        output{i}(1:m,colCnt{i}) = strainDataVec(:,cc);
                        
                        expCnt = expCnt+1;
                        expList{expCnt} = strainDataVec{8,cc};
                    end                    
                                            
                end                
            end            
        end        
    end
    %output{i} = result;
end
fprintf(' Done\n');

%% SAVE DATA TO EXCEL
fprintf('Saving results... ')
[FileName,PathName] = uiputfile([datestr(date, 'yyyymmdd') ' ' condition ' Experiments.xlsx'], 'Save results to Excel');
%[FileName,PathName] = uiputfile([datestr(date, 29) ' ' condition ' Experiments.xlsx'], 'Save results to Excel');
if isequal(FileName,0) || isequal(PathName,0)
    error('No valid input given. Results not saved.')
else
    warning('off','MATLAB:xlswrite:AddSheet');
    outputfile = fullfile(PathName,FileName);
    
    info = cell(6,max([numXLSfiles,numLvls,length(strains)]));
    info{1,1} = 'EXPERIMENTS';
    info(2,1:2) = {'Date:' datestr(date, 29)};
    info(3,1:2) = {'Condition' condition};
    info(4,1:numXLSfiles+1) = [{'Experiment dates:'} expDates];
    info(5,1:numLvls+1) = [{'Concentrations/Levels:'} num2cell(levels)];    
    info(6,1:length(strains)+1) = [{'Strains:'} strains];
    
    xlswrite(outputfile, info, 'Info');
    
    if length(idList) > expCnt
        expList(expCnt+1:end) = [];
    end
    
    xlswrite(outputfile, expList,'Experiment list');
    
    for i = 1:numLvls
        sheetName = [condition ' ' num2str(levels(i))];
        output{i}{1,1} = sheetName;        
        xlswrite(outputfile, output{i}, sheetName);
    end    
end
fprintf('Results saved to: %s\n\n',FileName);


end
