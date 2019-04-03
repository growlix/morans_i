function[output] = mL_moransI_master_acrossScales(valueMap,varargin)
% Computes Moran's I across spatial scales. Moran's I is a measure of
% spatial autocorrelation. See https://en.wikipedia.org/wiki/Moran%27s_I
% for a good explanation. Finds the set of all unique pairwise distances D
% between parcels (e.g. electrodes on an array). For each unique distance d
% in D, computes Moran's I for all parcels within distance ? d, and a
% distribution of shuffled values.
%
% OUTPUT:
%
% [output]: a structure of size (D) in which D is the number of unique
% distances between parcels in valueMap. Structure has fields:
% 'distances': a d x 1 vector of distances included in the the Moran's I
% calculation at the durrent distance.
% 'I': value of Moran's I.
% 'pVal': p-value
% 'I_shuff': distribution of shuffled Moran's I values.
%
% [valueMap] is matrix or cell array of two matrices in which each entry
% corresponds to a spatial location/parcel/area, and the value is some
% feature value of that spatial location (e.g. selectivity). [valueMap] is
% topographical; the spatial relationships of indices in [valueMap] reflect
% the spatial relationships of the features it represents. For a bivariate
% feature (e.g. spatial x- and y-position), the each entry in the cell
% array corresponds to one variable (e.g. valueMap{1} = x-position,
% valueMap{2} = y-position).
%
% Optional string/argument pairs:
%
% 'weightFun': a custom distance weighting function. Example = exponential
% weighting: @(d) exp(-1*(d-.4)) (exponential weighting). Default is no
% distance weighting (i.e. all included distances are weighted 1).
%
% 'percentile': a 2-tuple of values indicating the lower and upper
% percentiles of the confidence interval of shuffled values. Default = 95th
% percentile (i.e. [2.5 97.5]).
%
% 'electrodeSpacing': scalar. Space between electrodes on array. Default =
% 0.4.
%
% 'xMax': scalar. x-axis plot limit maximum. Default = 5.5
%
% 'plotMap': boolean. Determines whether valueMap is also plotted.
% Default = 0.
%
% 'nShuffles': scalar. Number of shuffle iterations to run. Default = 1000.

p = inputParser ;
p.addRequired('valueMap') ;
p.addParameter('weightFun',[]) ;
p.addParameter('percentile',[2.5 97.5]) ;
p.addParameter('electrodeSpacing',.4) ;
p.addParameter('xMax',5.5) ;
p.addParameter('plotMap',0) ;
p.addParameter('nShuffles',1000) ;

parse(p,valueMap,varargin{:}) ;

% If valueMap is a cell array, turn it into a regular array
if iscell(valueMap)
    valueMap = cat(3,valueMap{1},valueMap{2}) ;
end

% Determine if using bivariate data
bivariate = 0 ;
if size(valueMap,3) > 1
    bivariate = 1 ;
end

% Loop through each set of variables
% Compute distances
distanceMat = ...
    squareform(mL_distanceMat(valueMap(:,:,1))) ;

% Find unique distances
uniqueDistances = unique(distanceMat) ;
% Remove 0
uniqueDistances(uniqueDistances == 0) = [] ;
% Number of unique distances
nUniqueDists = length(uniqueDistances) ;

% Loop through, compute at each set of distances
for currDist = 1:nUniqueDists
    % Initialize current weight matrix
    currWeightMat = zeros(size(distanceMat)) ;
    % Included distances
    includedDists = uniqueDistances(1:currDist) ;
    % Indices of included distances
    currDistInds = ismember(distanceMat,includedDists) ;
    % Remove excluded distances
    currWeightMat(currDistInds) = distanceMat(currDistInds) ;
    
    % Distance weighting
    if ~isempty(p.Results.weightFun)
        currWeightMat = weightFun(currWeightMat) ;
    else
        currWeightMat = double(currWeightMat ~= 0) ;
    end
    
    % Weighting for zero distance should be 1, not zero
    currWeightMat(currWeightMat == inf) = 0 ;
    % Add distances to output structure
    output(currDist).distances = includedDists ;
    
    [output(currDist).I, output(currDist).pVal, output(currDist).I_shuff] = ...
        mL_wmData_moransI_withShuffle(valueMap,currWeightMat,p.Results.nShuffles) ;
    
    percentiles(:,:,currDist) = prctile(output(currDist).I_shuff,...
        p.Results.percentile,2) ;
end

iVals = [output.I] ;

% Plot it
figure ;
if p.Results.plotMap
    subplot(2,1,1)
end
hold on
lowerCI = squeeze(percentiles(1,1,:)) ;
upperCI = squeeze(percentiles(1,2,:)) ;
ciplot(lowerCI,upperCI,...
    p.Results.electrodeSpacing.*uniqueDistances,[.75 .75 .75]) ;
plot(.4.*uniqueDistances,iVals,'k','linewidth',2) ;
hold off
ylabel('Moran''s I') ;
xlabel('cluster radius (mm)') ;
xlim([0 p.Results.xMax]) ;
set(gca,'TickDir','out') ;
title(['Clustering of Preferred Locations Across Spatial Scales']) ;
% If also plotting valueMap
if p.Results.plotMap
    subplot(2,1,2)
    sanePColor(valueMap)
    %     set(gca,'XTick',[]) ;
    %     set(gca,'YTick',[]) ;
    axis square
    axis off
end