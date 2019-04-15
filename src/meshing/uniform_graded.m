function [mesh,P_VEC] = uniform_graded(L,hu,lambda_star,n,sigma,p_vec_count) 
% S.Groth 10/10/12
% Creates a mesh which is uniform to within lambda_star of ends, in these
% sections, the mesh is graded towards each end.

%           L = length of entire mesh
% lambda_star = length of graded sections at each end
%          hu = upper bound on uniform mesh step size
%           n = number of layers in graded sections
%       sigma = grading parameter

% M = ceil(L/hu);                        % number of elements, ceil is safer than round

% UNIFORM PORTION
L_uni = L-2*lambda_star;                 % length of uniform portion
M = round(L_uni/hu);                     % number of elements in uniform portion
if M==0 % then there is no uniform portion
    x_uni=[];
    h_uni=[];
    x_mid_uni=[];
else
    x_uni = L_uni/M*(0:1:M)' + lambda_star;  % element endpoints
    h_uni = diff(x_uni);                     % element length in uniform portion
    x_mid_uni = x_uni(1:end-1) + h_uni./2;   % midpoint of each element
end
p_vec_uni = ones(1,M)*max(p_vec_count);

% PORTION GRADED TOWARDS s=0
x_grad_0 = [0; sigma.^(n-1:-1:0)'*lambda_star];
h_0 = diff(x_grad_0);
x_mid_0 = x_grad_0(1:end-1) + h_0/2;
p_vec_left = p_vec_count;

% PORTION GRADED TOWARDS s=L
x_grad_L = L - flipud(x_grad_0);
h_L = diff(x_grad_L);
x_mid_L = x_grad_L(1:end-1) + h_L/2;
p_vec_right = fliplr(p_vec_count);

% CONCATENATE PORTIONS
if lambda_star==0
    x = x_uni';
    h = h_uni';
    x_mid = x_mid_uni';
    P_VEC = p_vec_uni;
else
    x = [x_grad_0(1:end-1);x_uni;x_grad_L(2:end)]';
    h = [h_0;h_uni;h_L]';
    x_mid = [x_mid_0;x_mid_uni;x_mid_L]';
    P_VEC = [p_vec_left p_vec_uni p_vec_right];
end

m = M + 2*n;  % total number of elements in mesh
mesh = struct('x',x,'h',h,'mid',x_mid,'m',m);
% keyboard