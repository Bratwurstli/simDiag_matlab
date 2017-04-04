function [basis values] = simDiag(matrices)
% SIMDIAG Calculate simultaneous eigenvalues and a joint basis to get upper right 
% triangular matrices.
%   [basis values] = simDiag([M1 M2]) calculates the simultaneous eigenvalues and 
%	returns them in a Matrix called "values" so that the first column of values 
%	contains the eigenvalues of M1 and the second column contains the eigenvalues 
%	of M2. On the other hand you will get upper right triangular matrices when 
%	you calculate the matrix product basis\M1*basis or basis\M2*basis. This also
%	works for k matrices: [basis values] = simDiag([M1 ... Mk]). 

% dimension is vector space dimension
dimension = length(matrices(:,1));
% kMatrices is number of matrices
kMatrices = length(matrices(1,:))/dimension;

% initialize return matrices with zeros
basis = zeros(dimension,dimension);
values = zeros(dimension,kMatrices);

% dimension 1 and you are done
if dimension==1
	basis = 1;
	values = matrices;
else
	% finds the first matrix for this algorithm see simDiag.pdf file
	startIndex = simDiagSelect(matrices);
	startMatrix = matrices(:,((startIndex-1)*dimension+1):((startIndex-1)*dimension+dimension));
	
	% uses Gram-Schmidt to get ONB from eigenvectors of the start matrix
	[V D] = simDiagEig(startMatrix);
	[Q R] = qr(V);
	
	% calculate dimensions of eigenspaces
	eigDimensions = unique(D);
	for iDimension = 1:length(eigDimensions)
		eigDimensions(iDimension) = length(find(D==eigDimensions(iDimension)));
	end
	
	% transform every matrix with ONB from above
	for iMatrix = 0:kMatrices-1
		matrices(:,(iMatrix*dimension+1):(iMatrix*dimension+dimension)) = Q'*matrices(:,(iMatrix*dimension+1):(iMatrix*dimension+dimension))*Q;
	end

	% put eigenvalues from start matrix which are fix from here into values matrix
	values(:,startIndex) = diag(matrices(:,(startIndex-1)*dimension+1:(startIndex-1)*dimension+dimension));

	% perform recursion with the other matrices 
	if kMatrices > 1
		index = setdiff(1:kMatrices,startIndex);
		for iDimension = 1:length(eigDimensions)
			if iDimension==1
				[basis(1:eigDimensions(1),1:eigDimensions(1)) values(1:eigDimensions(1),index)] = ...
					simDiag(matrices(1:eigDimensions(1), ...
						sort(repmat((index-1)*dimension,1,eigDimensions(1)))+ ...
							repmat(1:eigDimensions(1),1,kMatrices-1)));
			else
				[basis(sum(eigDimensions(1:iDimension-1))+1:sum(eigDimensions(1:iDimension)),sum(eigDimensions(1:iDimension-1))+1:sum(eigDimensions(1:iDimension))) values(sum(eigDimensions(1:iDimension-1))+1:sum(eigDimensions(1:iDimension)),index)] = ...
						simDiag(matrices(sum(eigDimensions(1:iDimension-1))+1:sum(eigDimensions(1:iDimension)), ...
							sort(repmat((index-1)*dimension,1,eigDimensions(iDimension)))+ ...
								repmat(1:eigDimensions(iDimension),1,kMatrices-1)+sum(eigDimensions(1:iDimension-1))));
			end
		end
		% TODO cant explain this at the moment because it is not documented in the simDiag.pdf 
		basis = Q*basis;
	else
		% TODO cant explain this at the moment because it is not documented in the simDiag.pdf
		basis = V;
	end
end

% nested functions

function matrixIndex = simDiagSelect(matrices)
% local function to get the best start matrix for the simultaneous diagonalization algorithm

% dimension is vector space dimension
dimension = length(matrices(:,1));
% kMatrices is number of matrices
kMatrices = length(matrices(1,:))/dimension;

% two eigenvalues will be seen equal when their distance is less than 10^-10
errorParam = 10;

% initialize eigenvalue frequency matrix
eigTableau = zeros(1,kMatrices);

% calculate values like it is described in simDiag.pdf
for iMatrix = 0:kMatrices-1
        matrix = matrices(:,(iMatrix*dimension+1):(iMatrix*dimension+dimension));
        eigValues = eig(matrix);
        eigValues = round(eigValues,errorParam);
        freqEig = mode(eigValues);
        freqInd = find(eigValues==freqEig);
        eigTableau(1,iMatrix+1) = length(freqInd);
end

% return index of matrix with lowest value
[Y matrixIndex] = min(eigTableau);

end

function [eigenVec eigenVal] = simDiagEig(matrix)
% local function with a small fix to the eig function for better performance with 
% the simultaneous diagonalization algorithm

% dimension is vector space dimension
dimension = length(matrix(:,1));

% two eigenvalues will be seen equal when their distance is less than 10^-10
errorParam = 10;

[eigenVec eigenVal] = eig(matrix,'vector');
eigenVal = round(eigenVal,errorParam);
eigenValSort = unique(sort(eigenVal));
permutation = zeros(1,dimension);
helper = 0;

for iValue = 1:length(eigenValSort)
        index = find(eigenVal==eigenValSort(iValue));
        indexLength = length(index);
        permutation(1,helper+1:helper+indexLength) = index;
        helper = helper + indexLength;
end

% returns eigenvalues and eigenvectors sorted by eigenvalue ocurrence instead of 
% eigenvalue norm 
eigenVal = sort(eigenVal);
eigenVec = eigenVec(:,permutation(1,:));

end

end
