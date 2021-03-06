function results = calcParameterInfluence(input,parWeightRefs)
% calcParameterInfluence - Calculates the influence of the growth
% parameters on the scores and ranking of the strains
%
% results = calcStrainRanking(input)
% results = calcStrainRanking(input,options)
%
% Required input
% input      Structure generated by calcStrainRanking
%
% Optional input
% options           Structure containing the following fields:
%  levels           Numeric vector with concentrations/levels, e.g. [1 1.5 3 6]. 
%
%  strains          Cell array of strings with names of strains to include.
%
%  parWeightRefs    Numeric vector with the weight for each parameter:                   
%                   1) Lag phase (h)
%                   2) Growth duration (h)
%                   3) Average growth rate (1/h)
%                   4) Number of generations	
%                   5) Maximum specific growth rate (1/h)
%                   Default is to use the weights from the rank analysis as
%                   reference. Weights will automatically be adjusted to a sum of five.
%                   
%  weightedTol      Calculate weighted tolerance scores, true or false (default)
%
% Output
% results       Structure with the following fields:
% 

%% HANLDE INPUT
% Set the property
if isfield(input,'robustness')
    property = 'Robustness';
elseif isfield(input,'tolerance');
    porperty = 'Tolerance';
else
    disp('Property to be calculated is not specified correctly. Program ended.')
    return
end

% Set default settings
levelList = input.leveList;
strainList = input.strainList;
parWeights = input.options.parWeights;
refState = input.options.refState;
weightedTol = input.options.false;   
    
    
    
    
    
    
    
    