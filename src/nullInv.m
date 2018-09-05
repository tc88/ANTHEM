function AInv = nullInv(A)
% NULLINV creates the pseudo-inverse of a diagonal matrix (or vector) A.
% In other words, it inverts every entry on the diagonal (or of the vector)
% that does not equal zero.
%
% Input:
%   A      matrix (or vector) to be inverted
%
% Output:
%   AInv   pseudo-inverse of the matrix (or vector)
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

[indm, indn, values] = find(A);
AInv = sparse(indm, indn, 1./values, size(A,1), size(A,2));

end