function[r, pVal, rShuff] = mL_mantel_test(valueMap,varargin)
% Mantel test computes a correlation, with the significance assesed by
% permutation testing, most often used to examine correlation between some
% feature and the spatial layout/position of those features (e.g. distance
% and selectivity). Distances are not independent of each other?changing
% the position of one object would change the distances between it and the
% other objects?thus significance is assesed using a permutation test.
% 
% See https://en.wikipedia.org/wiki/Mantel_test for more details.
%
% OUTPUT:
%
% r: Pearson's correlation
% pVal: p-value
% rShuff: shuffled correlation values
%
% INPUT:
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
% 'plot': boolean. Determines whether empirical and shuffled correlations
% are also plotted. Default = 1.
%
% 'nShuffles': scalar. Number of shuffle iterations to run. Default = 1000.

p = inputParser ;
p.addRequired('valueMap') ;
p.addParameter('plot',1) ;
p.addParameter('nShuffles',1000) ;

parse(p,valueMap,varargin{:}) ;

% If valueMap is a cell array, turn it into an n x d array, in which each
% row is a parcel and each column is a feature value.
if iscell(valueMap)
    parcelFeatureValues = [reshape(valueMap(:,:,1),[],1) ...
        reshape(valueMap(:,:,2),[],1)] ;
else
    parcelFeatureValues = valueMap(:) ;
end

% Compute physical distances between parcels on array
arrayDistances = mL_distanceMat(valueMap(:,:,1))' ;
% Compute feature distances
featureDistances = pdist(parcelFeatureValues)' ;
% Compute correlation
r = corrcoef(arrayDistances,featureDistances,'rows','complete') ;
r = r(1,2) ;

% Initialize vector of shuffled values
rShuff = nan(p.Results.nShuffles,1) ;
% Compute shuffled values
for i = 1:p.Results.nShuffles
    permutedDists = arrayDistances(randperm(length(arrayDistances))) ;
    ri = corrcoef(featureDistances,permutedDists,'rows','complete') ;
    rShuff(i) = ri(1,2) ;
end
% Compute p value
pVal = sum(rShuff>r)/p.Results.nShuffles ;
if pVal > .5
    pVal = abs(pVal-1) ;
end
pVal = pVal .*2 ;

% Plot results
if p.Results.plot
    figure('position',[1 787 200 550]);
    distributionPlot(rShuff,'showMM',3) ;
    hold on ;
    scatter(1,r,[],'filled','MarkerEdgeColor','none','MarkerFaceColor','r') ;
    ylabel('Pearson''s r') ;
    set(gca,'TickDir','out') ;
    set(gca,'XTickLabel',[])
    title({'distance-','feature','correlation'})
    legend({'null','median','empirical'}) ;
end