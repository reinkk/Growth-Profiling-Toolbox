function results = rankSensitivityAnalysis(results,property,options)



%% 
input = results;

numStrains = length(results.strainList);
numPars = length(results.parList);

% % Generate random seed of parameter weights
% s = 0:0.25:4.75;
% nseed = 10;
% seed = zeros(nseed*length(s),numPars);
% s2 = repmat(s,nseed,1);
% seed(:,1) = s2(:);
% n = size(seed,1);
% 
% for i = 2:4
%     remain = 5 - sum(seed,2);
%     seed(:,i) = zeros(n,1) + (remain-zeros(n,1)).*rand(n,1);
% end
% remain = 5 - sum(seed,2);    
% seed(:,5) = remain;
% 
% % Calculate the full matrix of weights
% N = 5*n;
% parWeights = zeros(N,5);
% cnt = 0;
% for i = 1:n
%     for j = 1:5
%         p = zeros(1,5);
%         while p(j) ~= 1;
%             p = randperm(5);
%         end
%         cnt = cnt+1;
%         parWeights(cnt,:) = seed(i,p);
%     end
% end    
% results.sensAnalysis.seed = seed;
% results.sensAnalysis.parWeights = parWeights;

N = 1000;
weightVector = 0 + (5-0).*rand(N,1);
%lower = [zeros(600,1);0.5*ones(500,1);ones(400,1);1.5*ones(400,1);2.5*ones(100,1)];
%upper = [0.5*ones(600,1);ones(500,1);1.5*ones(400,1);2.5*ones(400,1);5*ones(100,1)];
%parWeights(:,1) = lower + (upper-lower).*rand(2000,1);
weightMatrix = [weightVector zeros(N,4)];
for i = 2:4
    remain = 5 - sum(weightMatrix,2);
    weightMatrix(:,i) = zeros(N,1) + (remain-zeros(N,1)).*rand(N,1);
end
remain = 5 - sum(weightMatrix,2);
weightMatrix(:,5) = remain;

results.sensAnalysis.weightVector = weightVector;
results.sensAnalysis.weightMatrix = weightMatrix;


% for i = 1:N
%     parWeights(i,:) = parWeights(i,randperm(5));
% end    
% results.sensAnalysis.parWeights = parWeights;
% 
% 

% Run the Monte-Carlo simulation using the weights
%N = size(parWeights,1);
progress = (1:N)./N*100;
results.sensAnalysis.ranks = []; %zeros(N*numPars,numStrains);
results.sensAnalysis.parWeights = []; %zeros(N*numPars,numPars);

fprintf('Running sensitivity analysis...\n')
for j = 1:numPars
    fprintf('Parameter %d: ',j);
    p = zeros(1,numPars);
    while p(j) ~= 1
        p = randperm(numPars);
    end
    
    pw = weightMatrix(:,p);
    results.sensAnalysis.parWeights = [results.sensAnalysis.parWeights; pw];
    
    ranks = zeros(N,numStrains);
    for i = 1:N
        switch property
            case 'Robustness'
                input = calcRobustnessRanking(input,pw(i,:));
            case 'Tolerance'
                %parWeights(i,:)
                input = calcToleranceRanking(input,pw(i,:),options.refState,options.weightedTol);
        end
        sortedResults = input.sortedResults(:,1:2);
        s = sortrows(sortedResults,2);
        ranks(i,:) = cell2mat(s(:,1))';
        
        if i == find(progress >= 25,1,'first')
            fprintf('25%%... ');
        elseif i == find(progress >= 50,1,'first')
            fprintf('50%%... ');
        elseif i == find(progress >= 75,1,'first')
            fprintf('75%%... ');
        end
    end
    fprintf('100%%\n')
    results.sensAnalysis.ranks = [results.sensAnalysis.ranks; ranks];
end
fprintf('Done. \n')


results.sensAnalysis.medianRanks = zeros(numStrains,1);
results.sensAnalysis.medianRanksSDs = zeros(numStrains,1);
results.sensAnalysis.rankingStats = cell(numStrains,numPars);
results.sensAnalysis.parSensitivity = zeros(numStrains,numPars);

fprintf('Calculating results... ')
for s = 1:numStrains    
    % Calculate statistics and sensitivities
    max = TopExtreme(results.sensAnalysis.ranks(:,s));
    min = LowExtreme(results.sensAnalysis.ranks(:,s));    
    
    isIncl = results.sensAnalysis.ranks(:,s) >= min & results.sensAnalysis.ranks(:,s) <= max;
    
    results.sensAnalysis.medianRanks(s) = round(mean(results.sensAnalysis.ranks(isIncl,s)));
    results.sensAnalysis.medianRanksSDs(s) = std(results.sensAnalysis.ranks(isIncl,s));
    
    for p = 1:numPars        
        isWeightVector = ismember(results.sensAnalysis.parWeights(:,p),weightVector);
        
        ranks = results.sensAnalysis.ranks(isWeightVector,s);
        rankList = unique(ranks);
        counts = zeros(size(rankList));
        meanParWeights = zeros(size(rankList));
        %regrWeights = zeros(length(rankList),numPars);
        
        numPoints = length(rankList);
        
        for i = 1:numPoints
            isRank = ranks == rankList(i);
            counts(i) = sum(isRank);
            
            pw = weightVector(isRank);
            meanParWeights(i) = mean(pw);
        end
        
        denom = sum(counts)./size(counts,1);
        regrWeights = counts./denom;
        
        results.sensAnalysis.rankingStats{s,p} = [rankList counts meanParWeights regrWeights];
        
        if numPoints > 1
            mdl = LinearModel.fit(meanParWeights,rankList,'Weights',regrWeights);
            results.sensAnalysis.parSensitivity(s,p) = mdl.Coefficients{2,1};
            
        else
            results.sensAnalysis.parSensitivity(s,p) = 0;
        end
            
    end
end
fprintf('Done\n')

% -------------------------------------------------------------------------
function Y = TopExtreme(X)
%Maximum whisker length w. The default is a w of 1.5. Points are drawn as outliers if they are larger than q3 + w(q3 – q1) or smaller than q1 – w(q3 – q1), where q1 and q3 are the 25th and 75th percentiles, respectively. The default of 1.5 corresponds to approximately +/–2.7? and 99.3 coverage if the data are normally distributed. The plotted whisker extends to the adjacent value, which is the most extreme data value that is not an outlier.

P = prctile(X,[25,75]);
w = P(2) + 1.5*(P(2) - P(1));
if w > max(X);
    Y = max(X);
else
    Y = max(X(X <= w));
end

% -------------------------------------------------------------------------
function Y = LowExtreme(X)
%Maximum whisker length w. The default is a w of 1.5. Points are drawn as outliers if they are larger than q3 + w(q3 – q1) or smaller than q1 – w(q3 – q1), where q1 and q3 are the 25th and 75th percentiles, respectively. The default of 1.5 corresponds to approximately +/–2.7? and 99.3 coverage if the data are normally distributed. The plotted whisker extends to the adjacent value, which is the most extreme data value that is not an outlier.

P = prctile(X,[25,75]);
w =  P(1) - 1.5*(P(2) - P(1));
if w < min(X)
    Y = min(X);
else
    Y = min(X(X >= w));
end
% -------------------------------------------------------------------------
