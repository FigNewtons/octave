
% Returns a matrix of the distinct eigenvalues of n x n 
% matrix A and their corresponding algebraic multiplicity 

function evals = eigenvalues(A)
	
	% Check if A is square	
	if (size(A,1) != size(A,2))
		printf("Error: Not a square matrix.\n");
		return
	endif

	% Step 1: Compute characteristic polynomial f(t) of A
	coeff = charpoly(A);

	% Step 2: Get roots of f(t)
	eigs = roots(coeff);

	% There might be an issue with precision
	eigs = round(eigs .* 1000) ./ 1000;

	% Step 3: Format output / calculate multiplicities
	[u_eigs, ~, index] = unique(eigs);

	evals = [u_eigs, accumarray(index, 1)];

endfunction
