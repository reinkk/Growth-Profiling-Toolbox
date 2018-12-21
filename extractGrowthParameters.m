function ResultsMatrix = extractGrowthParameters(input, varargin)
% extractGrowthParameters - Extracts the growth parameters from a structure
% variable created by analyzeGrowthCurves or tuneGrowthParameters and saves
% the results to Excel.
%
% ResultsMatrix = extractGrowthParameters(input)
% ResultsMatrix = extractGrowthParameters(input, excelFlag)
%
% Required input
% input             structure variable created by analyseGrowthCurves or tuneGrowthParameters
%
% Optional input
% excelFlag         'true' (default) or 'false'. Write the data to an
%                    Excel-file.
%
% Output
% ResultsMatrix     cell array with extracted data

%% HANDLE INPUT
excelFlag = true;
if nargin > 1 
    if islogical(varargin{1})
        excelFlag = varargin{1};
    else
        disp('Second input argument has to be logical true or false.');
        return
    end
end
    

%% DETERMINE THE SIZE OF THE DATA
numExp = length(input.experimentIDs);
maxNumPeaks = 0;
for i = 1:numExp
    numPeaks = size(input.kinetics{i}.Results.mu_max_lreg,1);
    if numPeaks > maxNumPeaks
        maxNumPeaks = numPeaks;
        I = i;
    end
end

%% DEFINE HEADINGS AND PREALLOCATE SPACE
M1_heading = {'Strain','Condition','Concentration/Level','Unit','Date','Tray','Well'};
sM1 = length(M1_heading);

M2_heading = {'Lag phase (h)', 'End of growth phase (h)', 'Growth duration (h)', 'Average growth rate (1/h)', 'Number of generations','Maximum OD','Time of maximum OD (h)','Maximum specific growth rate (1/h)','Maximum specific growth rate SD (1/h)','Peak selected'};
sM2 = length(M2_heading);

growthData = input.kinetics{I}.Results;

mu_max_lreg = growthData.mu_max_lreg;
sM3 = numel(mu_max_lreg);
[m,n] = size(mu_max_lreg);
M3_heading = cell(m,n);
for i = 1:m
    M3_heading{i,1} = sprintf('Maximum slope %d (1/h)',i);
    M3_heading{i,2} = sprintf('Maximum slope %d SD (1/h)',i);
    M3_heading{i,3} = sprintf('Regression interval %d start (h)',i);
    M3_heading{i,4} = sprintf('Regression interval %d end (h)',i);
    M3_heading{i,5} = sprintf('Regression interval %d length (h)',i);
    M3_heading{i,6} = sprintf('Number of generations in regression interval %d',i);
    M3_heading{i,7} = sprintf('Number of points used in regression %d',i);
end
M3_heading_cols = M3_heading';
M3_heading_new = M3_heading_cols(:)';

mu_max_der = growthData.mu_max_der;
sM4 = numel(mu_max_der);
[m,n] = size(mu_max_der);
M4_heading = cell(m,n);
for i = 1:m
    M4_heading{i,1} = sprintf('Maximum derivative %d (1/h)',i);    
    M4_heading{i,2} = sprintf('Time of maximum derivative %d start (h)',1);  
end
M4_heading_cols = M4_heading';
M4_heading_new = M4_heading_cols(:)';

M1 = cell(numExp,sM1);
M2 = NaN(numExp,sM2);
M3 = NaN(numExp,sM3);
M4 = NaN(numExp,sM4);


%% EXTRACT THE DATA
for i = 1:numExp
    M1(i,:) = input.kinetics{i}.ExperimentInfo;
    
    lagPhase = input.kinetics{i}.Results.lagPhase;
    growthPhase = input.kinetics{i}.Results.growthPhase;
    mu_mean = input.kinetics{i}.Results.mu_mean;
    numGen = input.kinetics{i}.Results.numGen;  
    maxOD = input.kinetics{i}.Results.maxOD;  
    mu_max_ID = input.kinetics{i}.Results.mu_max_ID;
    mu_max_lreg = input.kinetics{i}.Results.mu_max_lreg;
    M2(i,:) = [lagPhase  growthPhase  mu_mean  numGen maxOD  mu_max_lreg(mu_max_ID,[1,2])  mu_max_ID];
        
    mu_max_lreg_cols = input.kinetics{i}.Results.mu_max_lreg';
    mu_max_lreg = mu_max_lreg_cols(:)';
    for j = 1:length(mu_max_lreg)        
        M3(i,j) = mu_max_lreg(j);    
    end
    
    mu_max_der_cols = input.kinetics{i}.Results.mu_max_der';
    mu_max_der = mu_max_der_cols(:)';
    for j = 1:length(mu_max_der)        
        M4(i,j) = mu_max_der(j);
    end
end

headings = [M1_heading M2_heading M3_heading_new M4_heading_new];
M = [M1 num2cell([M2 M3 M4])];
ResultsMatrix = [headings; M];

%% SAVE RESULTS TO EXCEL
if excelFlag
    [filename, pathname] = uiputfile(['Parameters ' input.condition ' ' datestr(date,29) '.xlsx'], 'Save results to Excel');
    if isequal(filename,0) || isequal(pathname,0)
        disp('Results not saved to Excel')
    else
        file = fullfile(pathname, filename);
        xlswrite(file, ResultsMatrix, 1);
        disp(['Results saved to Excel in ' file])
    end
end

end