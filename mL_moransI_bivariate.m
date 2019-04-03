function[I] = mL_moransI_bivariate(featureMat1,featureMat2,weightMat)
% Computes Moran's I of a spatial feature matrix [featureMat] using the weights
% [weightMat]. Moran's I is a measure of spatial autocorrelation. See
% https://en.wikipedia.org/wiki/Moran%27s_I for a good explanation.
%
% Put simply, for each index in matrix [featureMat], the difference between
% that value and each other value in the matrix is multiplied, and each of
% those products is weighted by [weightMat]. If you are only concerned with
% feature similarity between immediate neighbors, [weightMat] would contain
% 1s at all indices for adjacent indices in [featureMat] and 0s at all
% other indices. Note that [featureMat] is topographical; the spatial
% relationships of indices in [featureMat] reflect the spatial
% relationships of the features it represents. I =
% N./(SUMi(SUMj(WEIGHTij))).*...
% (SUMi(SUMj(WEIGHTij*(Xi-Xmean)*(Xj-Xmean))))/SUMi((Xi-Xmean)^2) , where N
% is the number of spatial units indexed by i and j, X is the variable of
% interest (ie feature value), Xmean is the mean of that variable (ie mean
% of featureMat), and WEIGHTij is the matrix of weights (ie weightMat).

if size(featureMat1) ~= size(featureMat2)
    error('Feature matrices are not the same size') ;
end

% Find nan indices in feature vector
[nanInds] = find(isnan(featureMat1) | isnan(featureMat2)) ;

% Vectorize feature matrix
featureVector1 = featureMat1(:) ;
featureVector2 = featureMat2(:) ;

% Number of non-nan indices (ie non-nan)
N = numel(featureMat1) - length(nanInds);

% Center (subtract mean) from matrix
featureVectorCentered1 = featureVector1 - nanmean(featureVector1) ;
featureVectorCentered2 = featureVector2 - nanmean(featureVector2) ;

% Set all nan indices in feature matrix to 0. Also set to 0 associated weights
% in weight matrix.
featureVectorCentered1(nanInds) = 0 ;
featureVectorCentered2(nanInds) = 0 ;
weightMat(nanInds,:) = 0 ;
weightMat(:,nanInds) = 0 ;

% Calculate sum of weighted feature deviations.
% We are rearranging
% SUMi(SUMj(WEIGHTij*(Xi-Xmean)*(Xj-Xmean)))to
% SUMij((WEIGHTij*(Xj-Xmean))*(Xi-Xmean)) and calculating
% (WEIGHTij*(Xj-Xmean)) first
weightedDeviationsI1 = (weightMat*featureVectorCentered1) ;
weightedDeviationsI2 = (weightMat*featureVectorCentered2) ;
% Now multiply by (Xi-Xmean) (ie centered feature vector)
weightedDeviationsIJ1 = weightedDeviationsI1.*featureVectorCentered1 ;
weightedDeviationsIJ2 = weightedDeviationsI2.*featureVectorCentered2 ;
% Now divide by SUMi((Xi-Xmean)^2) (ie sum of squared deviations)
weightedDeviations = nansum(weightedDeviationsIJ1+weightedDeviationsIJ2)/...
    nansum(featureVectorCentered1.^2+featureVectorCentered2.^2);
% Calculate sample size and weight normalization term
nwWeight = (N/sum(nansum(weightMat))) ;
% Now apply sample size and weight normalization term
I = nwWeight*weightedDeviations ;