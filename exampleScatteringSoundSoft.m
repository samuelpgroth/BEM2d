%-------------------------------------------------------------------------%
%
%        Time-harmonic wave scattering by a sound-soft polygon      
%
%        Solve the single-layer integral equation -  quick and easy
%
%        S.P.Groth 15th April 2019
%
%-------------------------------------------------------------------------%

clear
close all
addpath(genpath('src'))
%======================== PROBLEM PARAMETERS =============================%

xInc = [0 1000];  % location of point source
A = 1;                             % incident wave amplitude
d = [0 1];                        % incident wave vector
k0 = 10;                           % reference wavenumber
shape = 8;                         % specify shape. 0=square, 1=equilateral triangle, 2=hexagon
res = 2;%RES(i_p);                         % resolution, i.e. elements per wavelength
n_gauss = 5;                       % number of abscissae in Gauss quadrature
n_grading = 3;
sigma = 0.15;
c_star = pi/2;                     % controls distance from corner over which the mesh is graded. Set to 2*pi to get 1 wavelength from corner.
MU = 1;                            % MU in [0,1], MU=0 corresponds to a constant poly degree across mesh, MU>0 makes poly degree vary.
p_max = n_grading-1;               % maximum polynomial degree on mesh.
%=========================================================================%
d = d/norm(d);                     % normalise d
[n,tv,N,P,L] = geometry(d,shape);  % n unit norms, s tangents, N no. sides, P vertices
lambda=2*pi/k0;                    % compute wavelength in the exterior
lambda_star = c_star/k0;            % distance from corner over which mesh is graded.

hu=lambda/res;                     % upper bound on uniform mesh step size

if shape == 7
    N=1;
end

[p_vec] = mesh_poly_degree(MU,p_max,n_grading);% polynomial degree vector \bf{p} in notes.
p_vec_count = p_vec + 1;                       % tells us how many polynomials of different degree are on each element. Useful for loops.
%=================== Construct mesh/meshes and basis =====================%
tic
mesh = struct('x',cell(1),'h',cell(1),'mid',cell(1),'m',zeros(1,1));
basis = struct('side',zeros(1,1),'const',zeros(1,1),'deg',zeros(1,1),'phase',zeros(1,1),'where',zeros(1,1),'h',zeros(1,1),'which',zeros(1,1),'corner',zeros(1,1),...
    'mesh_no',zeros(1,1));
grid = zeros(N,1);
for i = 1:N
    [mesh(i),P_VEC{i}] = uniform_graded(L(i),hu,lambda_star,n_grading,sigma,p_vec_count);
    grid(i) = mesh(i).m;
    for j = 1:grid(i)
        for p = 1:P_VEC{i}(j)
            num_temp=0;
            for ii_temp=1:i-1
                num_temp=num_temp+sum(P_VEC{ii_temp});
            end
            index = num_temp+sum(P_VEC{i}(1:j-1))+p;
            scaling = sqrt(2)/sqrt(mesh(i).h(j));
            basis(index) = struct('side',i,'const',scaling,...
                'deg',p-1,'phase',0,'where',mesh(i).mid(j),'h',mesh(i).h(j),'which',0,'corner',0,'mesh_no',j);
        end
    end
end
%======================= Quadrature parameters ===========================%
[x_gauss,W]=gauleg(n_gauss);


%% Assemble RHS
 Ui_vec = zeros(length(basis),1);

 for i = 1:length(basis)
     side = basis(i).side;            % which side does basis live on
     [X_gauss_coords] = gauss_coords(x_gauss,P(side,:),tv(side,:),basis(i).where,basis(i).h); % coords of Gauss quad points     
     Ui = ui(X_gauss_coords,k0,d,A);  % incident field at Gauss quad points
     Ui_vec(i) = 0.5*basis(i).h*(W.*basis_func_2(x_gauss,basis(i)))'*Ui;           % Gaussian quadrature   0.5*basis.h from Gauss quad formula
 end
 b=-Ui_vec;
 
 %% Assemble matrix
[AA] = single_layer_operator(n_gauss, tv, P, basis, k0, N, L);

%% Solve system
dudn=AA\b;       
toc

%% Plot solution on boundary and in domain
h_plot = lambda/20; 

for i = 1:N
    [mesh_plot(i)] = uniform_mesh(L(i),h_plot);
    grid(i) = mesh_plot(i).m;
end

% Evaluate solution
clear Acols
for i = 1:length(basis)
    x = mesh_plot(basis(i).side).mid;
    temp = zeros(1,sum(grid));
    a = basis(i).where-basis(i).h/2;
    b = basis(i).where+basis(i).h/2;
    where = find(x>a&x<=b);
    x_eval = (x(where)-a)*2/basis(i).h-1;
    values = basis_func_eval(x_eval,basis(i),P,tv);
    temp(where+sum(grid(1:basis(i).side-1))) = values;
    Acols(:,i) = temp;
end
dudnB=Acols*dudn;

% Boundary solution
plot(real(dudnB))

nGauss = 10;
sigma = 0.15;
n_grading = nGauss;
[xG,wG]=gauss(nGauss);

Nsides=N;

Xquad=[];
w_quad=[];
n_count=0;
for i=1:Nsides
    [x,w] = hp_1d_quad_composite(0,L(i)/2,10,xG,wG',n_grading);
    x_eval = [x;L(i)-flipud(x)];
    w_eval = [w;flipud(w)];
    x_quad{i} = x_eval;
    w_quad = [w_quad;w_eval];
    Xquad = [Xquad;P(i,:)+tv(i,:).*x_eval]; 
    n_count = n_count+length(x_eval);
    N_COUNT(i) = length(x_eval);
end

% scatter(Xquad(:,1),Xquad(:,2),'x');  % plot quad points on boundary to check

%% Evaluation points in the domain
Ndom=201;
x_yo = linspace(-7,7,Ndom);
y_yo = x_yo;

[X,Y] = meshgrid(x_yo,y_yo);

X_v(:,1) = reshape(X,Ndom^2,1);
X_v(:,2) = reshape(Y,Ndom^2,1);

% quad coords replicated Ndom^2 times
X_new(:,1)=reshape(repmat(Xquad(:,1)',Ndom^2,1),length(Xquad)*Ndom^2,1);
X_new(:,2)=reshape(repmat(Xquad(:,2)',Ndom^2,1),length(Xquad)*Ndom^2,1);

% domain points replicated length(Xquad) times
X_rv = repmat(X_v,length(Xquad),1);

diff = X_new-X_rv;
dist = sqrt(diff(:,1).^2+diff(:,2).^2);

temp = 1i/4*besselh(0,1,k0*dist);

hey=reshape(temp,Ndom^2,length(w_quad));

% keyboard
clear Acols
% Evaluate solution
for i = 1:length(basis)
    x = x_quad{basis(i).side};
    temp = zeros(1,n_count);
    a = basis(i).where-basis(i).h/2;
    b = basis(i).where+basis(i).h/2;
    where = find(x>a&x<=b);
    x_eval = (x(where)-a)*2/basis(i).h-1;
    values = basis_func_eval(x_eval,basis(i),P,tv);
%     temp(where+(basis(i).side-1)*length(x)) = values;
    temp(where+sum(N_COUNT(1:basis(i).side-1))) = values;
    Acols(:,i) = temp;
end
dudnref=Acols*dudn;

int = hey*(w_quad.*dudnref);

field = reshape(int,Ndom,Ndom);

uinc = ui(X_v,k0,d,1);
uincField = reshape(uinc,Ndom,Ndom);
uTot = uincField+field;

figure
pcolor(x_yo,y_yo,real(uTot))
shading interp
colormap(bluewhitered)
% caxis([-1 1])
colorbar
axis equal
axis tight
axis off

hold on
for i=1:Nsides
    plot([P(i,1) P(i+1,1)],[P(i,2) P(i+1,2)],'-k','LineWidth',2), hold on
end

set(gcf,'color','w')
% export_fig irregular_scatterer.pdf