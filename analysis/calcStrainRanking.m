function results = calcStrainRanking(datafile,property,varargin)
% calcStrainRanking - Calculates the scores and ranking of either
% Performance and Robustness or Conformance and Tolerance
%
% results = calcStrainRanking(datafile,property)
% results = calcStrainRanking(datafile,property,options)
%
% Required input
% datafile          There are three possible inputs:
%                   1) Empty vector will open a dialog window to choose an
%                      Excel-file generated by extractGrowthCurves or
%                      tuneGrowthParameters
%                   2) String with full path and file name of an Excel-file
%                      generated by analyzeGrowthCurves or tuneGrowthParameters
%                   3) Cell array generated by analyzeGrowthParameters
%                      or tuneGrowthKinetics
% property          String indicating which property to calculate, 'Robustness' or 'Tolerance'
%
% Optional input
% options           A structure contining any of the following fields:
%  levels           Numeric vector with concentrations/levels, e.g. [1 1.5 3 6]. 
%
%  strains          Cell array of strings with names of strains to include. Give an empty vector to include all strains.
%
%  parWeights       Numeric vector with the weight for each parameter:                   
%                   1) Lag phase (h)
%                   2) Growth duration (h)
%                   3) Average growth rate (1/h)
%                   4) Number of generations	
%                   5) Maximum specific growth rate (1/h)
%                   Default is to use equal weight for all paramteres, i.e. [1 1 1 1 1]
%                   Weights will automatically be adjusted to a sum of five.
%
%  refState         By default the function selects level = 0 as the reference state. 
%                   If this level does not exist the reference state is
%                   selected according to the following possibilities:
%                   1) 'best point' (default) - the best parameter value measured in
%                      any concentration/level is selected as the reference
%                   2) 'min' - minimum level 
%                   3) 'individual' - the best state is selected
%                      indvidually for each strain
%                   4) selected level specified as a scalar, e.g. [5.5]
%                   This setting is only applicable when property is set to 'Tolerance'
%                   
%  weightedTol      Calculate weighted tolerance scores, true or false (default)
%
%  runPIA           Perform Parameter Influence Analysis on the ranking, true (default) or false 
%
% Output
% results          Structure with the following fields:
%                   


%% HANLDE INPUT
results = [];
readExcelFlag = true;

if isempty(datafile)
    [filename, pathname] = uigetfile( {'*.xls;*.xlsx', 'EXCEL Files (*.xls, *.xlsx)';'*.*','All Files (*.*)'},'Open Excel data file');
    datafile = fullfile(pathname, filename);
    
    if isequal(filename,0) || isequal(pathname,0)
        disp('No valid file selected. Program ended.')
        return
    end
elseif iscell(datafile)
    readExcelFlag = false;    
end

switch property
    case {'Robustness','Tolerance'}
        %Do nothing
    otherwise
        disp('Property to be calculated is not specified correctly. Program ended.')
        return
end

% Set default settings
levelList = [];
strainList = {};
parWeights = [];
refState = 'best point';
weightedTol = false;
runPIA = true;

if nargin > 2
    if ~isstruct(varargin{1})
        warning('Options is not a structure variable. Default settings used.');
    else
        opt = varargin{1};        
        if isfield(opt,'levels'), levelList = unique(opt.levels); if size(levelList,2) == 1, levelList = reshape(levelList,1,length(levelList)); end; end
        if isfield(opt,'strains'), strainList = unique(opt.strains); if size(strainList,1) == 1, strainList = reshape(strainList,length(strainList),1); end; end
        
        if isfield(opt,'parWeights')
            parWeights = opt.parWeights;
            if length(parWeights) ~= 5
                disp('Parameter weights has to be of a length of five.');
                return;
            end
            
            if size(parWeights,2) == 1
                parWeights = reshape(parWeights,1,length(parWeights));
            end
            
            if sum(parWeights) ~= 5
                parWeights = parWeights./(sum(parWeights)./length(parWeights));
            end        
        end
        
        if isfield(opt,'refState')
            refState = opt.refState;
            if ~isscalar(refState)
                switch refState
                    case {'min','individual','best point'}
                        %Do nothing
                    otherwise
                        disp('Reference state is not specified correctly. Analysis ended.');
                        return
                end
            end
        end
        
        if isfield(opt,'weightedTol')
            weightedTol = opt.weightedTol; 
            if ~islogical(weightedTol)
                disp('The weighted tolerance option is not specified correctly. Analysis ended.');
                return
            end
        end     
        
        if isfield(opt,'runPIA')
            runPIA = opt.runPIA; 
            if ~islogical(runPIA)
                disp('The PIA option is not specified correctly. Analysis ended.');
                return
            end
        end   
        
    end
end

%options = struct('levelList',{levelList},'strainList',{strainList},'parWeights',{parWeights},'refState',{refState},'weightedTol',{weightedTol},'pia',{pia});

%% IMPORT DATA
if readExcelFlag
    [~,~,dataMatrix] = xlsread(datafile,1);    
end

strains = dataMatrix(2:end,1);
%condition = dataMatrix{2,2};
levels = dataMatrix(2:end,3); levels = reshape([levels{:}],size(levels));
%unit = dataMatrix{2,4};

allParameters = {'Lag phase (h)' 'Growth duration (h)'	'Average growth rate (1/h)'	'Number of generations'	'Maximum specific growth rate (1/h)'};
allParameterData = zeros(size(dataMatrix,1)-1,5);

for i = 1:length(allParameters)
    isPar = ismember(dataMatrix(1,:),allParameters{i});
    allParameterData(:,i) = cell2mat(dataMatrix(2:end,isPar));
end

%% DEFINE LISTS
if isempty(levelList)    
    levelList = unique(levels)';
end
%numLevels = length(levelList);
isLevels = ismember(levels,levelList);

if isempty(strainList)
    strainList = unique(strains);
end
numStrains = length(strainList);
isStrains = ismember(strains,strainList);

if isempty(parWeights)
    parWeights = ones(1,5);
end
%numPars = length(parameters);
%isPars = ismember(parameters,parList);

% Define new lists and data matrix
levels = levels(isLevels & isStrains);
strains = strains(isLevels & isStrains);
%parameters = parameters(isPars);
dataMatrix = allParameterData(isLevels & isStrains, :);

%% CALCULATE GROUP STATISTICS
% Calculate statistics for each strain and condition
[means,sds,gnames,numElements] = grpstats(dataMatrix,{strains,levels},{'mean','std','gname','numel'});
%levelsStats = str2double(gnames(:,2));

% Remove statistics based on less than 2 elements
isLow = numElements < 2;
means(isLow) = NaN;
sds(isLow) = NaN;
relSDs = sds./means*100;

%% GATHER RESULTS
results.strainList = unique(gnames(:,1),'stable');
results.parList = allParameters;
results.levelList = levelList;
results.options = struct('parWeights',{parWeights},'refState',{refState},'weightedTol',{weightedTol});

results.stats.means = means;
results.stats.sds = sds;
results.stats.relSDs = relSDs;
results.stats.gnames = gnames;
results.stats.numElements = numElements;

%% CALCULATE SCORES AND RANKING
numPars = length(parWeights);
switch property
    case 'Robustness'
        results = calcRobustnessParameterScores(results,numPars);
        results = calcRobustnessRanking(results,parWeights);
    case 'Tolerance'        
        results = calcToleranceParameterScores(results,numPars,refState);
        results = calcToleranceRanking(results,parWeights,weightedTol);
end

%% GENERATE OUTPUT
parScores = results.parameters.scores;
parErrors = results.parameters.errors;

switch property
    case 'Robustness'
        levelScores = results.performance.scores;
        levelErrors = results.performance.errors;
        
        propScore = results.robustness.scores;
        propError = results.robustness.errors;       
        
    case 'Tolerance'        
        levelScores = results.conformance.scores;
        levelErrors = results.conformance.errors;
        
        propScore = results.tolerance.scores;
        propError = results.tolerance.errors;
end

M1 = cell(numStrains,1);
for i = 1:numStrains
    M1{i,1} = sprintf('%.1f � %.1f', results.RVA.Ranks(i),results.RVA.RankBounds(i));
end

M2 = cell(numStrains,2);
for i = 1:numStrains
    M2{i,1} = sprintf('%0.1f � %0.1f', propScore(i),propError(i));
    M2{i,2} = sprintf('%0.1f � %0.1f', results.RVA.Scores(i),results.RVA.ScoreErrors(i));
end

parScores2 = [];
parErrors2 = [];
for p = 1:numPars
    parScores2 = [parScores2 parScores{p}]; %#ok<AGROW>
    parErrors2 = [parErrors2 parErrors{p}]; %#ok<AGROW>
end

pValues = results.RVA.pValues;
I = results.RVA.sortingIDX;

RVAresults1 = [num2cell((1:numStrains)') M1(I,:) results.strainList(I,:) num2cell(results.RVA.Ranks(I,:)) num2cell(results.RVA.RankBounds(I,:)) num2cell(pValues(I,:)) num2cell(levelScores(I,:)) M2(I,:) cell(numStrains,1) num2cell(levelErrors(I,:)) cell(numStrains,1) num2cell(parScores2(I,:)) cell(numStrains,1) num2cell(parErrors2(I,:)) ];
RVAresults2 = sortrows(RVAresults1,3);

results.RVAresults = [RVAresults1; cell(3,size(RVAresults1,2)); RVAresults2]; 

%% PERFORM PIA
if runPIA
    PIA = pia(results);    
    
    % Generate ranking order
    r = zeros(numStrains,1);
    for i = 1:numStrains
        r(i) = find(I == i);
    end
    
    % Generate ouput
    n = numStrains + 1;
    M1H = [{'Strain' 'Ranking order'} results.parList {'' 'Strain' 'Ranking order'} results.parList cell(1,n-1-2*(2+numPars))];
    M1D1 = [results.strainList num2cell(r) num2cell(PIA.InfluenceScores)];
    M1D2 = sortrows(M1D1,2);
    M1D = [M1D1 cell(numStrains,1) M1D2 cell(numStrains,n-1-2*(2+numPars))];
    M1 = [M1H; M1D; cell(1,n)];
    
    % Collect ranks for each parameter   
    M2 = [];
    for p = 1:numPars
        s = {sprintf('Ranks: %s',results.parList{p})};
        M2H1 = [s cell(1,n-1)];
        M2H2 = [{'Weight'} results.strainList' cell(1,n-1-numStrains)];
        M2D = [num2cell(PIA.weightVector)  num2cell(PIA.Ranks{p})];
        M2 = [M2; M2H1; M2H2; M2D; cell(1,n)];
    end
    
    % Collect scores for each parameter
    M3 = [];
    for p = 1:numPars
        s = {sprintf('Scores: %s',results.parList{p})};
        M3H1 = [s cell(1,n-1)];
        M3H2 = [{'Weight'} results.strainList' cell(1,n-1-numStrains)];
        M3D = [num2cell(PIA.weightVector)  num2cell(PIA.Ranks{p})];
        M3 = [M3; M3H1; M3H2; M3D; cell(1,n)];
    end
    
    results.PIA = PIA;
    results.PIAresults = [M1; M2; M3;];
end

%% AUXILIARY FUNCTIONS
%-------------------------------------------------------------------------
function results = calcRobustnessParameterScores(results,numPars)
%% DEFINE LISTS
levelList = results.levelList;
numLevels = length(levelList);

gnames = results.stats.gnames;
levelsStats = str2double(gnames(:,2));

%% CACLULATE PARAMETERS SCORES
means = results.stats.means;
relSDs = results.stats.relSDs;

parScores = cell(1,numPars);
parErrors = cell(1,numPars);

for p = 1:numPars
    for c = 1:numLevels
        levelData = means(levelsStats == levelList(c),p);
        data_relSDs = relSDs(levelsStats == levelList(c),p)./100;    
        
        % Calculate scores
       switch p
           case {1,2}
               levelData = 1./levelData;
       end
       maxValue = max(levelData);
       x = [0 maxValue];
       y = [0 20];
       
       b = polyfit(x,y,1);
       scores = polyval(b,levelData);
       scores(isnan(scores)) = 0;
       
       errors = scores.*data_relSDs;
       errors(isnan(errors)) = 0;
       
       parScores{p}(:,c) = scores;        
       parErrors{p}(:,c) = errors;
    end
end

%% GATHER RESULTS
results.parameters.scores = parScores;
results.parameters.errors = parErrors;

%-------------------------------------------------------------------------
function results = calcToleranceParameterScores(results,numPars,refState)
%% DEFINE LISTS
strainList = results.strainList;
levelList = results.levelList;
gnames = results.stats.gnames;

numStrains = length(strainList);
numLevels = length(levelList);

%% CALCULATE REFERENCE POINT AND RATIOS
means = results.stats.means;
relSDs = results.stats.relSDs;

refData = cell(numStrains,3);
parScores_tmp = cell(numStrains,1);
parErrors_tmp = cell(numStrains,1);

for s = 1:numStrains    
    isStrain = ismember(gnames(:,1),strainList{s});    
    dataS = means(isStrain,:);
    dataS_relSDs = relSDs(isStrain,:)./100;
    
    % Set the reference state for each strain
    if isscalar(refState)
        stateID = find(levelList == refState);
    elseif isequal(refState,'min')
        [~,stateID] = min(levelList);
    elseif isequal(refState,'individual')
        bestPoint = [min(dataS(:,1:2)) max(dataS(:,3:5))];
        normDataS = bsxfun(@rdivide, dataS, bestPoint);
        normDataS(:,[1,2]) = 1./normDataS(:,[1,2]);
        
        x = [0 1];
        y = [0 20];
        b = polyfit(x,y,1);
        scoresS = polyval(b,normDataS);
        scoresS(isnan(scoresS)) = 0;
        [~,stateID] = max(sum(scoresS,2));
    else
        stateID = [];
        [bestPoint1,bestPointIDs1] = min(dataS(:,1:2));    
        [bestPoint2,bestPointIDs2] = max(dataS(:,3:5));
        
        bestPoint = [bestPoint1 bestPoint2];
        bestPointIDs = [bestPointIDs1 bestPointIDs2];
    end
    
    % Calculate parameter ratios and maxValues for parameter scoring
    if ~isempty(stateID)
        refData{s,1} = strainList{s};
        refData{s,2} = levelList(stateID);
        refData{s,3} = dataS(stateID,:);        
        
        sdS = results.stats.sds(isStrain,:);
        bestValues = dataS(stateID,:) + [-1 -1 1 1 1].*sdS(stateID,:);
        maxValues = bestValues./dataS(stateID,:);
        maxValues(:,[1,2]) = 1./maxValues(:,[1,2]);
        
        refData{s,4} = maxValues;
    else
        refData{s,1} = strainList{s};
        refData{s,2} = levelList(bestPointIDs);
        refData{s,3} = bestPoint;        
        
        maxValues = ones(1,numPars);
        refData{s,4} = maxValues;
    end 
    
    % Calculate parameters relative to the reference state
    ratios = bsxfun(@rdivide,dataS,refData{s,3});
    %     ratios_errors = ratios.*dataS_relSDs;
    
    % Calculate ratio scores
    parScores_tmp{s} = zeros(size(dataS));
    for p = 1:numPars
        maxValue = maxValues(1,p);
        
        x = [0 maxValue];
        y = [0 20];
        b = polyfit(x,y,1);
        
        switch p
            case {1,2}
                parScores_tmp{s}(:,p) = polyval(b,1./ratios(:,p));
            otherwise
                parScores_tmp{s}(:,p) = polyval(b,ratios(:,p));
        end
    end
    
    parScores_tmp{s}(isnan(parScores_tmp{s})) = 0;
    if ~isempty(stateID)
        parScores_tmp{s}(stateID,:) = 20*ones(1,numPars);
    end
    
    parErrors_tmp{s} = parScores_tmp{s}.*dataS_relSDs;
    parErrors_tmp{s}(isnan(parErrors_tmp{s})) = 0;
end

%% GATHER RESULTS
parScores = cell(1,numPars);
parErrors = cell(1,numPars);

for p = 1:numPars
    parScores{p} = zeros(numStrains,numLevels);
    parErrors{p} = zeros(numStrains,numLevels);
    
    for i = 1:numStrains
        parScores{p}(i,:) = parScores_tmp{i}(:,p)';
        parErrors{p}(i,:) = parErrors_tmp{i}(:,p)';
    end
end

results.referenceData = refData;
results.parameters.scores = parScores;
results.parameters.errors = parErrors;
