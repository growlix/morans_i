function[I pVal I_shuff] = mL_moransI_withShuffle...
    (arrayValues,weightMat,nReps)
% [arrayValues] is an m x n topographical matrix of some feature as it varies
% across a space, [weightMat] is a m*n x m*n marix of weights for each pairwise
% comparison

bivariate = false ;
% If arrayValues has a third dimension, we are doing a bivariate Moran's I
if size(arrayValues,3) > 1
    bivariate = true ;
    I = mL_moransI_bivariate(arrayValues(:,:,1),arrayValues(:,:,2),weightMat) ;
else
    I = mL_moransI(arrayValues,weightMat) ;
end

if ~exist('nReps','var')
    nReps = 1000 ;
end
I_shuff = nan(1,nReps) ;

if bivariate
    % Repeating myself is not a good idea, but i don't feel like making
    % this pretty right now
    arrayValuesVector1 = arrayValues(:,:,1) ;
    arrayValuesVector1 = arrayValuesVector1(:) ;
    arrayValuesVector2 = arrayValues(:,:,2) ;
    arrayValuesVector2 = arrayValuesVector2(:) ;
    nValues = length(arrayValuesVector1) ;
    arraySize = size(arrayValues(:,:,1)) ;
    for repN = 1:nReps
        shuffVals1 = reshape(randsample(arrayValuesVector1,nValues),arraySize) ;
        shuffVals2 = reshape(randsample(arrayValuesVector2,nValues),arraySize) ;
        I_shuff(repN) = mL_moransI_bivariate(shuffVals1,shuffVals2,weightMat) ;
    end
else
    arrayValuesVector = arrayValues(:) ;
    nValues = length(arrayValuesVector) ;
    arraySize = size(arrayValues) ;
    % Shuffle
    for repN = 1:nReps
        shuffVals = reshape(randsample(arrayValuesVector,nValues),arraySize) ;
        I_shuff(repN) = mL_moransI(shuffVals,weightMat) ;
    end
end

pVal = sum(rShuff>r)/nReps ;
if pVal > .5
    pVal = abs(pVal-1) ;
end
pVal = pVal .*2 ;