function[I] = mL_moransI(featureMat,weightMat)
% Computes Moran's I of a spatial feature matrix [featureMat] using the weights
% [weightMat]. Moran's I is a measure of spatial autocorrelation. Put simply,
% for each index in matrix [featureMat], the difference between that value and
% each other value in the matrix is multiplied, and each of those products is
% weighted by [weightMat]. If you are only concerned with feature similarity
% between immediate neighbors, [weightMat] would contain 1s at all indices for
% adjacent indices in [featureMat] and 0s at all other indices. Note that
% [featureMat] is topographical; the spatial relationships of indices in
% [featureMat] reflect the spatial relationships of the features it represents.
% I = N./(SUMi(SUMj(WEIGHTij))).*...
% (SUMi(SUMj(WEIGHTij*(Xi-Xmean)*(Xj-Xmean))))/SUMi((Xi-Xmean)^2) , where N is
% the number of spatial units indexed by i and j, X is the variable of interest
% (ie feature value), Xmean is the mean of that variable (ie mean of
% featureMat), and WEIGHTij is the matrix of weights (ie weightMat).


% Set all nan indices in feature matrix to 0. Also set to 0 associated weights
% in weight matrix.
[nanInds] = find(isnan(featureMat)) ;

% Vectorize feature matrix
featureVector = featureMat(:) ;

% Number of non-nan indices (ie non-nan)
N = numel(featureMat) - length(nanInds);

% Center (subtract mean) from matrix
featureVectorCentered = featureVector - nanmean(featureVector) ;

% Set all nan indices in feature matrix to 0. Also set to 0 associated weights
% in weight matrix.
featureVectorCentered(nanInds) = 0 ;
weightMat(nanInds,:) = 0 ;
weightMat(:,nanInds) = 0 ;

% Calculate sum of weighted feature deviations. 
% We are rearranging
% SUMi(SUMj(WEIGHTij*(Xi-Xmean)*(Xj-Xmean)))to
% SUMij((WEIGHTij*(Xj-Xmean))*(Xi-Xmean)) and calculating
% (WEIGHTij*(Xj-Xmean)) first
weightedDeviationsI = (weightMat*featureVectorCentered) ;
% Now multiply by (Xi-Xmean) (ie centered feature vector)
weightedDeviationsIJ = weightedDeviationsI.*featureVectorCentered ;
% Now divide by SUMi((Xi-Xmean)^2) (ie sum of squared deviations)
weightedDeviations = nansum(weightedDeviationsIJ)/nansum(featureVectorCentered.^2);
% Calculate sample size and weight normalization term
nwWeight = (N/sum(nansum(weightMat))) ;
% Now apply sample size and weight normalization term
I = nwWeight*weightedDeviations ;