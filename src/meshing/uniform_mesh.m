% uniform_mesh.m
% Constructs a uniform mesh on a side with length L. hu is the upper bound
% on the mesh element length.

function[mesh] = uniform_mesh(L,hu) 

% M = ceil(L/hu);                 % number of elements
% x = L/M*(0:1:M)';               % element endpoints
% h = L/M*ones(M,1);              % element length
% x_mid=(x(2:end)+x(1:end-1))/2;  % midpoint of each element
% 
% mesh = struct('x',x,'h',h,'mid',x_mid,'m',M);

M = round(L/hu);                % number of elements
x = L/M*(0:1:M)';               % element endpoints
h = L/M*ones(M,1);              % element length
x_mid=(x(2:end)+x(1:end-1))/2;  % midpoint of each element

mesh = struct('x',x,'h',h,'mid',x_mid,'m',M);