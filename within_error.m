function [errorBars, normData] = within_error(data, alphaVal)
%WITHIN_ERROR: Compute standard error of the mean (SEM) or a confidence
%interval (CI) for within-subject designs. 
%
%Normal error bars also reflect the between-subject variability, which 
%plays no role in statistical tests incorporating repeated measures 
%(e.g. repeated measures ANOVAs). Normal error bars therefore reflect \
%sample variance that is too large, and not indicative of the statistical 
%test used.
%
%This function implements a simple correction to return error bars that
%only reflect the within-subject variability. It is based on the 
%recommendations in:
%
%Cousineau, D. (2005). 
%Confidence intervals in within-subject designs: A simpler solution to Loftus and Masson?s method. 
%Tutorial in Quantitative Methods for Psychology, 1(1), 42?45.
%
%Morey, R. (2008). 
%Confidence intervals from normalized data: A correction to Cousineau (2005). 
%Tutorial in Quantitative Methods for Psychology, 4(2), 61?64. 
%
%Usage:
%
% [errorBars, normData] = WITHIN_ERROR(data, alphaValVal)
%
%Inputs:
%
% -data: matrix where each row has data from one participant, and each
% column has data from one level in the design (i.e. a within-subject
% condition)
%
% -alphaVal: significance level to test at, to compute corresponding confidence
% interval. For example, will return a 95% confidence interval for an input
% of 0.05. If empty/unspecified, will return a standard error of the mean
% instead.
%
% Outputs:
%
% -errorBars: vector of within-subject CIs / SEMs, one for each condition.
%
% -normData: normalized data matrix (with between-subject variability
% removed).

subjMean = mean(data,2); % for each participant, mean over all within-subject conditions
grandMean = mean(mean(data)); % mean of entire data table
normData = bsxfun(@minus, data, subjMean) + grandMean; % remove between-subject variability, keeping condition (column) means intact (Cousineau et al., 2005)

SEM = std(normData,0,1) ./ sqrt(size(data,1)); %standard error of the mean
corrSEM = size(data,2)/(size(data,2)-1); % correction factor to remove bias in sample variance (Morey et al., 2008)

if ~exist('alphaVal', 'var') || isempty(alphaVal) % for SEMs
    errorBars = SEM * corrSEM;
else % for CIs
    tVal = tinv(1-alphaVal/2,size(data,1)-1); % t-value to obtain 95% confidence interval if alpha = 0.05
    errorBars = tVal * SEM * corrSEM;
end