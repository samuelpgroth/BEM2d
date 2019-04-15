% mesh_poly_degree.m
% S.P.Groth 24/09/12
% Function which defines the degree vector P defined as in notes (in notes
% it is written \mathbf{p}

function[P] = mesh_poly_degree(mu,p,n)
% 0<= mu <= 1 is a refinement parameter, mu=0 correspondes to a constant
% degree across the mesh.
% p is maximum degree, n is number of mesh elements.

P = zeros(1,n);
for i = 1:n-1
    P(i) = p - floor(mu*(n+1-i)*p/n);
end
P(n) = p;