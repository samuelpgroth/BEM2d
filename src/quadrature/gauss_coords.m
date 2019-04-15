function[X_gauss_coords] = gauss_coords(absc,P,tv,basis_mid,basis_h)
% [X_gauss_coords] = gauss_coords(absc,P,tv,basis_mid,basis_h)

% This function converts Gaussian quadrature abscissae on the interval
% [-1,1] into coordinates on support of basis function, i.e. interval
% [-h,h]
% Note: absc = abscissae

n_gauss = length(absc);                       % number of Gaussian quadrature points
x_mid = P + tv*basis_mid;                     % coordinates of centre of basis support
X_gauss_coords = repmat(x_mid,n_gauss,1)+...
    [tv(1)*absc tv(2)*absc]*basis_h*0.5;      % coordinates of Gauss quad points
