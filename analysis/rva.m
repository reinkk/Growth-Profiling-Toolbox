function RVA = rva(propScores,levelScores,levelErrors,strainList,numLevels)
% rva - Rank Variablity Analysis
%
% RVA = rva(propertyScore,levelScores,levelErrors,numStrains,numLevels);
%
% This file is used by calcRobustnessRanking and calcToleranceRanking.
%
% Required input
% propScores   Property scores, ie Robustness or Tolerance scores
% levelScores  Level scores, ie Performance or Conformance scores
% levelErrors  Level errors, ie SD values for the Performance or Conformance scores
% strainList   Cell array with names of strains analyzed
% numLevels    Number of inhibitory levels analyzed
%
% Output
% RVA           Structure with the following fields:
%   ranks       Matrix of size 1000-by-numStrains with samples of rankings
%   propScores  Matrix of size 1000-by-numStrains with samples of the property score
%   strainPDFs  Cell of size numStrains-by-2 with pdf for the normally
%               distributed property scores (column 1) and the custom pdf
%               for the rank values
%   Ranks       Vector of size numStrains-by-1 with rexpected rank value
%   RankBounds  Vector of size numStrains-by-1 with range of the rank value
%   Scores
%   ScoreErrors 
%   sortingIDX  Vector of size numStrains-by-1 with index numbers for sorting the original data
%   pValues     Matrix of size numStrains-by-2 with p-values for (1) the rank range, (2) being among Top10 strains and (3) among Bottom 10 strains
%   

%% DEFINE VARIABLES
numStrains = length(strainList);

%% RUN RVA
N = 1000;
ranks = zeros(N,numStrains);
scores = zeros(N,numStrains);

rng(1,'twister');
for i = 1:N    
    levelScores_tmp = levelErrors.*randn(numStrains,numLevels)+levelScores;
    propScores_tmp = mean(levelScores_tmp,2);
    [~,I] = sort(propScores_tmp,1,'descend');
    sortedResults = [num2cell(1:numStrains)' strainList(I,:)];
    s = sortrows(sortedResults,2);
    
    ranks(i,:) = cell2mat(s(:,1)');
    scores(i,:) = propScores_tmp';
end

RVA.ranks = ranks;
RVA.propScores = scores;

%% CALCULATE THE STATISTICS
% Determine the range of possible scores for each rank value
stats = (1:numStrains)';
for i = 1:numStrains
    incl = ranks == i;
    stats(i,2) = mean(scores(incl));
    stats(i,3) = std(scores(incl));
end
stats(:,4) = stats(:,2) + 1.96*stats(:,3); %Upper 95% CI
stats(:,5) = stats(:,2) - 1.96*stats(:,3); %Lower 95% CI

% Determine the range of rank values for each strain, the probablility for the range
% and p-value for being in Top10 and Bottom10
strainRankRanges = zeros(numStrains,2);
Scores = zeros(numStrains,1);
ScoreErrors = zeros(numStrains,1);
pValues = zeros(numStrains,3);
strainPDFs = cell(numStrains,2);

%L = stats(numStrains-10+1,5);
%U = stats(10,4);
limit = max([100 ceil(max(scores(:)))]);
xi = (0:0.5:limit)';
ri = (1:numStrains+1)';
for i = 1:numStrains
    [muhat,sigmahat] = normfit(scores(:,i));
    Scores(i,1) = muhat;
    ScoreErrors(i,1) = sigmahat;
    
    pdfN = normpdf(xi,muhat,sigmahat);
    strainPDFs{i,1} = [xi pdfN];
   
    % Calculate the probability for finding a score within the range associated with each rank value by integrating under the pdf
    pdfEst = zeros(numStrains+1,1);
    for j = 1:numStrains
        incl = xi >= stats(j,5) & xi <= stats(j,4);        
        pdfEst(j) = trapz(xi(incl),pdfN(incl));
    end
    %Normalize values to make sure the total integral of the P-values is 1
    pdfC = pdfEst./dint(ri,pdfEst);
    strainPDFs{i,2} = [ri pdfEst pdfC];    

    % Find rank range
    a = find(stats(:,5) <= muhat+sigmahat,1,'first');
    b = find(stats(:,4) >= muhat-sigmahat,1,'last');
    strainRankRanges(i,1) = a;
    strainRankRanges(i,2) = b;
    
    % Calculate p-values
%     if strainRankRanges(i,1) == strainRankRanges(i,2)
%         pValues(i,1) = pdfEst(incl)./sum(pdfEst);
%     else
%         incl = ri >= strainRankRanges(i,1) & ri <= strainRankRanges(i,2)+1;    
%         pValues(i,1) = trapz(pdfC(incl));                
%     end
    pValues(i,1) = dint(ri(a:b+1),pdfC(a:b+1));
    
    pValues(i,2) = dint(ri(ri <= 11), pdfC(ri <= 11)); %p-value for being in Top10
    pValues(i,3) = dint(ri(ri > numStrains-10), pdfC(ri > numStrains-10)); %p-value for being in Bottom10
end

RVA.Scores = Scores;
RVA.ScoreErrors = ScoreErrors;
RVA.strainPDFs = strainPDFs;
RVA.Ranks = mean(strainRankRanges,2);
RVA.RankBounds = RVA.Ranks - strainRankRanges(:,1);
RVA.pValues = pValues;

%% SORT RESULTS
[~,I] = sortrows([RVA.Ranks RVA.RankBounds pValues(:,1) propScores], [1 2 -3 -4]);
RVA.sortingIDX = I;