
% Given an n x n matrix A, computes the coefficients
% of A's characteristic polynomial using the Faddeev-Leverrier
% method

% The algorithm is as follows:
%
%	B_1 = A 						c_1 = trace(B_1)
%	B_2 = A(B_1 - c_1 * I)			c_2 = trace(B_2)
%	...
%	B_i = A(B_(i-1) - c_i * I)		c_i = trace(B_i)
%	...
%	B_n = A(B_(n-1) - c_(n-1) * I)  c_n = trace(B_n)
%
%	where I = n x n identity matrix 
%	and   c_i = coefficient at the ith component
%	
%   The resulting characteristic polynomial is: f(t) = x^n - c_1 x^(n-1) - ... - c_(n-1) x - c_n
%
%	Link: http://mathfaculty.fullerton.edu/mathews/n2003/FaddeevLeverrierMod.html

function coeff = charpoly(A)
	n = size(A, 1);
	coeff = zeros(1, size(A,1) + 1);

	% Initialize starting values
	I = eye(n,n); 
	B = 2 .* I;
	c = 1;
	coeff(1) = 1;

	for i = 1:n
		B = A * (B - c .* I);
		c = trace(B) / i; 		
		coeff(i + 1) = -1 * c;
	endfor
endfunction
