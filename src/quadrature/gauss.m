% gauss.m  Gauss quadrature nodes and weights on [-1,1]
%          L. N. Trefethen 6/03

% This function computes nodes (column vector x) and
% weights (row vector w) for the N-point Gauss-Legendre
% quadrature formula:
%
% INT_{-1}^1 f(s) ds ~= w*f(x)
%
% x and w are obtained from a tridiagonal eigenvalue problem as
% proposed by Golub & Welsch (Math. Comp. 1969), an idea previously
% considered by by Goertzel 1954, Wilf 1960, and Gordon 1968; see
% Gautschi, Orthogonal Polynomials: Computation and Approximation, 
% Oxford 2004.  This code comes from Trefethen, Spectral Methods
% in MATLAB, 2000.  It is regrettably slow since MATLAB does not
% take advantage of tridiagonal structure for eigenvalue problems.

  function [x,w] = gauss(N)
  beta = .5./sqrt(1-(2*(1:N-1)).^(-2));   % 3-term recurrence coeffs.
  T = diag(beta,1) + diag(beta,-1);       % set up Jacobi matrix
  [V,D] = eig(T);                         % the eigenvalue problem
  x = diag(D); [x,i] = sort(x);           % nodes
  w = 2*V(1,i).^2;                        % weights
