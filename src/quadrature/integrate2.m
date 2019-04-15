% integrate2.m
% Samuel Groth 14/05/12

% Function to compute double integral in matrix/vector form 
% Double integral arises from inner product of integral operator (applied to
% basis function rho_y) and another basis function rho_x

% change from integrate.m: corrected so that integrates at each of gauss
% points, not just one as before (which was an error when using basis_func_2.m, ok for basis_func.m).


function[integral] = integrate2(ko,x_gauss,basis_x,basis_y,W)

% ko = kernel output, e.g. from B_kernel.m. Is vector of length n_gauss^2
n_gauss = length(x_gauss);
mat = reshape(ko,n_gauss,n_gauss); 
mat = mat.';
hx = basis_x.h; % length of element rho_x (rho_i)
hy = basis_y.h; % length of element rho_y (rho_j)

rho_nodes_x = basis_func_2(x_gauss,basis_x); % rho_x (rho_i) at Gauss nodes
rho_nodes_y = basis_func_2(x_gauss,basis_y);  % rho_y (rho_j) at Gauss nodes

integral = (0.5*hx*0.5*hy)*(W.*rho_nodes_x).'*mat*(W.*rho_nodes_y); % 0.5*hx*0.5*hy from Gauss quad formula
