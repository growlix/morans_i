function[dists] = mL_distanceMat(mat,distance)
% Calculates the distance between each element of a matrix [mat], using
% euclidean (default) or other (see documentation for pDist) methods. [mat]
% is either a matrix, in which case size(mat) is used for the computations,
% a single number n, in which case the adjacency matrix is computed for a
% matrix of size n x n, or a two-element vector, in which case the
% adjacency matrix is computed for a matrix of size mat(1) x mat(2).
% [distance] is a string indicating what kinds of distance is used (see
% pdist documentation for details); defaults to 'euclidean'. The distances
% are arranged in the order (2,1), (3,1), ..., (m,1), (3,2), ..., (m,2),
% ..., (m,m?1)). D is commonly used as a dissimilarity matrix in clustering
% or multidimensional scaling. To save space and computation time, D is
% formatted as a vector. However, you can convert this vector into a square
% matrix using the squareform function so that element i, j in the matrix,
% where i < j, corresponds to the distance between objects i and j in the
% original data set.

% Parse input
[r,c] = size(mat);                        %# Get the matrix size
if r == 1 && c == 1
    r = mat ;
    c = mat ;
elseif (r == 1 && c == 2) || (r == 2 && c == 1)
    r = mat(1) ;
    c = mat(2) ;
end

[x, y] = meshgrid(1:c, 1:r) ;

if exist('distance','var')
   dists = pdist([x(:) y(:)],distance) ; 
else
    dists = pdist([x(:) y(:)]) ;
end