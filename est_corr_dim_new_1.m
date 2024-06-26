% EST_CORR_DIM_NEW_1: estimate correlation dimension using four different
% methods. 
% 
% Uses the model
% c(s) ∝ s^(D-1),
% where D is correlation dimension and c(s) is correlation at distance s.
% 
% REQUIRED INPUTS:
% ss - Vector of network distances
% nn - Number of node pairs at each network distance
% 
% OPTIONAL INPUTS:
% DLims - Upper and lower bound for correlation dimension, entered as the
%         two-vector [DLower, DUpper] 
% 
% OUTPUTS:
% DVec - Vector of estimated correlation dimensions
% sMaxVec - Vector of estimated upper cutoffs sMax
% D2Mat - D2Mat(i, j) records network correlation dimension D estimated
%         using method i with upper cutoff sMax chosen using method j
% codeCell - Code for method of estimation
% descriptionCell - Description of method of estimation
% DMat - DMat(i, j) records estimate of network correlation dimension D
%        made using method i using upper cutoff number j
% objFuncMat - The function which is minimised (either locally or globally)
%              to choose the upper cutoff sMax
% bVec - Vector of estimated constants of proportionality b in
%        p(s) = b*s^(D-1), integer s ∈ [1, s_max],
%        where 1 = ∑_{s = 1}^{s_max} p(s) (i.e., p(s) is a probability)
%        Number of distinct nodes within distance of a typical node can be
%        estimated as
%        (N - 1)*c(s) = (sum(nn1)/N)*b*s^(D - 1), integer s ∈ [1, s_max]
% bMat - bMat(i, j) records estimate of constants of proportionality b
%        made using method i using upper cutoff number j
% 
% Example use:
% D = 2; N0 = 1000; k = 10; p = 0;
% A = small_world_manhattan(N0, k, D, p); N = size(A, 1);
% [ss, nn] = count_distances(A);
% DLims = [1, Inf]; [DVec, sMaxVec, D2Mat, codeCell, descriptionCell, DMat, objFuncMat, bVec, bMat] = est_corr_dim_new(ss, nn, DLims);
% 
%
% Moore et al. (2023), "Correlation dimension in empirical networks"
%
% Jack Moore, 2023
%
function [DVec, sMaxVec, D2Mat, codeCell, descriptionCell, DMat, objFuncMat, bVec, bMat] = est_corr_dim_new_1(ss, nn, varargin)
if (numel(varargin) == 0)
    DLims = [-Inf, Inf];
else
    DLims = varargin{1};
end
DLower = DLims(1); DUpper = DLims(2);

codeCell = {'CE', 'KS', 'RP', 'RC'};
descriptionCell = {...
    'Maximum coding efficiency/minimum coding inefficiency',...
    'Maximum likelihood corresponding to final local minimum in upper cut-off of Kolmogorov-Smirnoff distance',...
    'Linear regression to log-log plot of correlation',...
    'Linear regression to log-log plot of correlation sum',...
    };
numDTypes = numel(codeCell);

num_s = numel(ss);

DMat = NaN(numDTypes, num_s);
objFuncMat = NaN(numDTypes, num_s);
bMat = NaN(numDTypes, num_s);

objFuncMat(3, :) = -nn;
objFuncMat(4, :) = -nn;

options = optimoptions('fminunc', 'Display', 'off');
for ii_s = 2:num_s
    ss1 = ss(1:ii_s);
    nn1 = nn(1:ii_s);
    NN1 = cumsum(nn1);
    fun = @(g) -log_like(g, ss1, nn1);
    if sum(nn1) > 0
        [dML, ~] = fminunc(fun, 1, options);
    else
        dML = NaN;
    end
    DLimImposed = 0;
    if (dML + 1 < DLower); dML = DLower - 1; DLimImposed = 1; end
    if (dML + 1 > DUpper); dML = DUpper - 1; DLimImposed = 1; end
    negLogLike = fun(dML);
    scaledNegLogLike = negLogLike - log(1 + max(ss1) - min(ss1));
    if ~isequal(size(scaledNegLogLike), [1, 1]); scaledNegLogLike = NaN; end
    objFuncMat(1, ii_s) = scaledNegLogLike;%Objective function for CE
    unnormedPVec = ss1.^dML;
    bMat(1:2, ii_s) = 1/sum(unnormedPVec);
    
    DMat(1:2, ii_s) = dML + 1;%Dimension estimated from CE and KS
    if (~DLimImposed) && (ii_s <= 2)
        KS1 = 0;
    else
        M1 = sum(nn1);
        n1Fit = @(a) a.^dML*(sum(nn1)/sum(ss1.^dML));
        nn1Fit = n1Fit(ss1);
        nn1FitCum = cumsum(nn1Fit);
        nn1Cum = cumsum(nn1);
        KS1 = max(abs(nn1Cum - nn1FitCum), [], 'includenan')/M1;
    end
    objFuncMat(2, ii_s) = KS1;%Objective function for KS
    
    if (numel(ss1) >= 2)
        linFit = polyfit(log(ss1), log(nn1), 1);
        DRP = linFit(1) + 1;
        if (DRP < DLower); DRP = DLower; end
        if (DRP > DUpper); DRP = DUpper; end
        DMat(3, ii_s) = DRP;
        bMat(3, ii_s) = exp(linFit(2))/sum(nn1);
        
        linFit = polyfit(log(ss1), log(NN1), 1);
        DRC = linFit(1);
        if (DRC < DLower); DRC = DLower; end
        if (DRC > DUpper); DRC = DUpper; end
        DMat(4, ii_s) = DRC;
        bMat(4, ii_s) = exp(linFit(2))/sum(nn1);
    end
end

DVec = NaN(numDTypes, 1); sMaxVec = NaN(numDTypes, 1); bVec = NaN(numDTypes, 1);
D2Mat = NaN(numDTypes, numDTypes);%D2Mat(i, j) records D estimated using method i with sMax chosen using method j.
for iiDType = 1:numDTypes
    objFunc = objFuncMat(iiDType, :);
    try
    if (iiDType ~= 2)%Final global minimum:
        [~, ii_s] = min(fliplr(objFunc)); ii_s = numel(objFunc) - (ii_s - 1);
    else
        %Final final local minimum:
        [~, locMinInd] = find_local_minima(objFunc);
        ii_s = locMinInd(end);
    end
    DVec(iiDType) = DMat(iiDType, ii_s); sMaxVec(iiDType) = ss(ii_s); bVec(iiDType) = bMat(iiDType, ii_s);
    D2Mat(:, iiDType) = DMat(:, ii_s);
    catch
    end
end
end

% Log-likelihood per observation for model
% c(s) = s^d,
% where c(s) is correlation at distance s.
% 
% d is exponent (d + 1 is correlation dimension)
% ss is vector of distances from s = 1 to s = s_max
% nn is vector of counts of distances
%
function logLikePerObs = log_like(d, ss, nn)
unnormedPMat = ss.^d;
PMat = unnormedPMat./sum(unnormedPMat);
logLikePerObs = log(PMat)*(nn')/sum(nn);
end