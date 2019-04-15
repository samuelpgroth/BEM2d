function[S] = single_layer_operator(n_gauss, tv, P, basis, k0, N, L)

[x_gauss,W]=gauleg(n_gauss);

vec=[];
 vec_sing=[];
 for i = 1:n_gauss
     vec = [vec i*ones(1,n_gauss)]; %#ok<AGROW>
 end

nGaussGen = 40;
[x_sing_shift, y_sing_shift, r_temp, w_temp] = DuffyLog(-1, 1, nGaussGen, 2);  % generate points between -1 and 1 for basis_func_2 input
[x_s, y_s, r_s, w_s] = DuffyLog(0, 1, nGaussGen, 1);

progressbar;
 
 S = zeros(length(basis));

for i = 1:length(basis)            % rho_i
    side_i = basis(i).side;        % which side is rho_i on
    ai = basis(i).where-basis(i).h/2;  % first end of element (going anti-clockwise)
    bi = basis(i).where+basis(i).h/2;  % other end of element
    [i_gauss_coords] = gauss_coords(x_gauss,P(side_i,:),tv(side_i,:),basis(i).where,basis(i).h); % coords of Gauss quad points for rho_i
    X = i_gauss_coords(vec,:);     % replicates each coord in i_gauss_coords n_gauss times
    for j = 1:length(basis)        % rho_j is basis which kernels are applied to
        side_j = basis(j).side;    % which side is basis rho_j on
        if basis(j).deg==0
            aj = basis(j).where-basis(j).h/2; % first end of element j (going anti-clockwise)
            bj = basis(j).where+basis(j).h/2; % other end of element j
            [j_gauss_coords] = gauss_coords(x_gauss,P(side_j,:),tv(side_j,:),basis(j).where,basis(j).h); % coords of Gauss quad points for rho_j
            Y = repmat(j_gauss_coords,n_gauss,1);
            if side_i==side_j
                if basis(i).where==basis(j).where % C_kernel is singular
                    S_ker_temp = SLpotSing(basis(i).h*r_s',k0);
                    S(i,j) = basis(i).h*basis(j).h*sum(S_ker_temp.'.*w_s.*basis_func_2(x_sing_shift,basis(i)).*basis_func_2(y_sing_shift,basis(j)));
                elseif abs(bi-aj)<5e-15    % elements touch at one end (xb touches ya)
                    xleft = basis(i).where-basis(i).h/2;
                    hx=basis(i).h;
                    hy = basis(j).h;
                    yleft = xleft+hx;
                    [x_yo,y_yo,r_yoBIAJ,w_BIAJ] = RquadLog( ai, bi, bj,  nGaussGen);
                    x_BIAJ = -1+2*(x_yo-xleft)/hx; y_BIAJ = -1+2*(y_yo-yleft)/hy;
                    S_ker_temp = SLpotSing(r_yoBIAJ',k0);
                    S(i,j) = sum(S_ker_temp.'.*w_BIAJ.*basis_func_2(x_BIAJ,basis(i)).*basis_func_2(y_BIAJ,basis(j)));
                elseif abs(ai-bj)<5e-15   % elements touch at one end (xa touches yb)
                    xleft = basis(i).where-basis(i).h/2;
                    hx=basis(i).h;
                    hy = basis(j).h;
                    yleft = xleft-hy;
                    [y_yo,x_yo,r_AIBJ,w_AIBJ] = RquadLog( aj, bj, bi,  nGaussGen);
                    x_AIBJ = -1+2*(x_yo-xleft)/hx; y_AIBJ = -1+2*(y_yo-yleft)/hy;
                    S_ker_temp = SLpotSing(r_AIBJ',k0);
                    S(i,j) = sum(S_ker_temp.'.*w_AIBJ.*basis_func_2(x_AIBJ,basis(i)).*basis_func_2(y_AIBJ,basis(j)));
                elseif abs(ai-bj)>5e-15 && abs(ai-bj)<0.15
                    [x_ho,y_ho,r_ho,w_ho] = near_sing_2d(ai,bi,aj,bj,x_gauss,W,16);
                    x_ho = -1+2*(x_ho-ai)./basis(i).h; y_ho = -1+2*(y_ho-aj)./basis(j).h;
                    S_ker_temp = SLpotSing(r_ho',k0);
                    S(i,j) = sum(S_ker_temp.'.*w_ho.*basis_func_2(x_ho,basis(i)).*basis_func_2(y_ho,basis(j)));
                    
                elseif abs(aj-bi)>5e-15 && abs(bi-aj)<0.15
                    [x_ho,y_ho,r_ho,w_ho] = near_sing_2d(ai,bi,aj,bj,x_gauss,W,16);
                    x_ho = -1+2*(x_ho-ai)./basis(i).h; y_ho = -1+2*(y_ho-aj)./basis(j).h;
                    S_ker_temp = SLpotSing(r_ho',k0);
                    S(i,j) = sum(S_ker_temp.'.*w_ho.*basis_func_2(x_ho,basis(i)).*basis_func_2(y_ho,basis(j)));  
                else  % not touching nor near
                    S_ker_temp = 1i/4*besselh(0,1,k0*sqrt((X(:,1)-Y(:,1)).^2 + (X(:,2)-Y(:,2)).^2));
                    S(i,j) = integrate2(S_ker_temp,x_gauss,basis(i),basis(j),W);
                end
            else % i.e on different sides
                if side_i==modulo(side_j-1,N) && abs(basis(i).where+basis(i).h/2-L(side_i))<5e-15 && abs(basis(j).where-basis(j).h/2)<5e-15
                    % elements touch at corner (xb touches ya)
                    
                    hx=basis(i).h;
                    xleft = -hx;
                    hy = basis(j).h;
                    yleft = 0;
                    [x_yo,y_yo,r_XBYA,w_XBYA] = RquadLog( -hx, 0, hy,  nGaussGen);
                    [ymx_coords_XBYA] = y_minus_x_calc(x_yo,y_yo,tv,side_i,side_j,N);
                    x_XBYA = -1+2*(x_yo-xleft)/hx; y_XBYA = -1+2*(y_yo-yleft)/hy;
                    S_ker_temp = SLpotSing(r_XBYA',k0);
                    S(i,j) = sum(S_ker_temp.'.*w_XBYA.*basis_func_2(x_XBYA,basis(i)).*basis_func_2(y_XBYA,basis(j)));
                elseif side_j==modulo(side_i-1,N) && abs(basis(j).where+basis(j).h/2-L(side_j))<5e-15 && abs(basis(i).where-basis(i).h/2)<5e-15
                    % elements touch at corner (xa touches yb)
                    xleft = 0;
                    hx=basis(i).h;
                    hy = basis(j).h;
                    yleft = xleft-hy;
                    [y_yo,x_yo,r_XAYB,w_XAYB] = RquadLog( -hy, 0, hx,  nGaussGen);
                    [ymx_coords_XAYB] = y_minus_x_calc(x_yo,y_yo,tv,side_i,side_j,N);
                    x_XAYB = -1+2*(x_yo-xleft)/hx; y_XAYB = -1+2*(y_yo-yleft)/hy;
                    S_ker_temp = SLpotSing(r_XAYB',k0);
                    S(i,j) = sum(S_ker_temp.'.*w_XAYB.*basis_func_2(x_XAYB,basis(i)).*basis_func_2(y_XAYB,basis(j)));
                elseif side_i==modulo(side_j-1,N) && (L(side_i)-bi+aj)<0.15 && (L(side_i)-bi+aj)>5e-15
                    [x_ho,y_ho,r_XBYA,w_XBYA] = near_sing_2d(ai-L(side_i),bi-L(side_i),aj,bj,x_gauss,W,16);
                    [ymx_coords_XBYA] = y_minus_x_calc(x_ho,y_ho,tv,side_i,side_j,N);
                    x_XBYA = -1+2*(x_ho-(ai-L(side_i)))./basis(i).h; y_XBYA = -1+2*(y_ho-aj)./basis(j).h;
                    r_XBYA = sqrt(sum((ymx_coords_XBYA.').^2));
                    
                    S_ker_temp = SLpotSing(r_XBYA,k0);
                    S(i,j) = sum(S_ker_temp.'.*w_XBYA.*basis_func_2(x_XBYA,basis(i)).*basis_func_2(y_XBYA,basis(j)));
                elseif side_j==modulo(side_i-1,N) && (L(side_j)-bj+ai)<0.15 && (L(side_j)-bj+ai)>5e-15
                    [x_ho,y_ho,r_XAYBnear,w_XAYBnear] = near_sing_2d(ai,bi,aj-L(side_j),bj-L(side_j),x_gauss,W,16);
                    [ymx_coords_XAYBnear] = y_minus_x_calc(x_ho,y_ho,tv,side_i,side_j,N);
                    x_XAYBnear = -1+2*(x_ho-(ai))./basis(i).h; y_XAYBnear = -1+2*(y_ho-aj+L(side_j))./basis(j).h;
                    r_XAYBnear = sqrt(sum((ymx_coords_XAYBnear.').^2));
                    S_ker_temp = SLpotSing(r_XAYBnear,k0);
                    S(i,j) = sum(S_ker_temp.'.*w_XAYBnear.*basis_func_2(x_XAYBnear,basis(i)).*basis_func_2(y_XAYBnear,basis(j))); 
                else
                    S_ker_temp = 1i/4*besselh(0,1,k0*sqrt((X(:,1)-Y(:,1)).^2 + (X(:,2)-Y(:,2)).^2));
                    S(i,j) = integrate2(S_ker_temp,x_gauss,basis(i),basis(j),W);
                end
            end
        else  % degree is >0 and can re-use kernels
            if side_i==side_j
                if basis(i).where==basis(j).where % C_kernel is singular
                    S(i,j) = basis(i).h*basis(j).h*sum(S_ker_temp.'.*w_s.*basis_func_2(x_sing_shift,basis(i)).*basis_func_2(y_sing_shift,basis(j)));
                    
                elseif abs(bi-aj)<5e-15    % elements touch at one end (xb touches ya)
                    S(i,j) = sum(S_ker_temp.'.*w_BIAJ.*basis_func_2(x_BIAJ,basis(i)).*basis_func_2(y_BIAJ,basis(j)));
                elseif abs(ai-bj)<5e-15   % elements touch at one end (xa touches yb)
                    S(i,j) = sum(S_ker_temp.'.*w_AIBJ.*basis_func_2(x_AIBJ,basis(i)).*basis_func_2(y_AIBJ,basis(j)));
                elseif abs(ai-bj)>5e-15 && abs(ai-bj)<0.15
                    S(i,j) = sum(S_ker_temp.'.*w_ho.*basis_func_2(x_ho,basis(i)).*basis_func_2(y_ho,basis(j)));
                elseif abs(aj-bi)>5e-15 && abs(bi-aj)<0.15
                    S(i,j) = sum(S_ker_temp.'.*w_ho.*basis_func_2(x_ho,basis(i)).*basis_func_2(y_ho,basis(j)));
                else
                    S(i,j) = integrate2(S_ker_temp,x_gauss,basis(i),basis(j),W);
                end
            else % i.e on different sides
                if side_i==modulo(side_j-1,N) && abs(basis(i).where+basis(i).h/2-L(side_i))<5e-15 && abs(basis(j).where-basis(j).h/2)<5e-15
                    % elements touch at corner (xb touches ya)
                    S(i,j) = sum(S_ker_temp.'.*w_XBYA.*basis_func_2(x_XBYA,basis(i)).*basis_func_2(y_XBYA,basis(j)));
                elseif side_j==modulo(side_i-1,N) && abs(basis(j).where+basis(j).h/2-L(side_j))<1e-15 && abs(basis(i).where-basis(i).h/2)<1e-15
                    % elements touch at corner (xa touches yb)
                    S(i,j) = sum(S_ker_temp.'.*w_XAYB.*basis_func_2(x_XAYB,basis(i)).*basis_func_2(y_XAYB,basis(j)));
                elseif side_i==modulo(side_j-1,N) && (L(side_i)-bi+aj)<0.15 && (L(side_i)-bi+aj)>5e-15
                    S(i,j) = sum(S_ker_temp.'.*w_XBYA.*basis_func_2(x_XBYA,basis(i)).*basis_func_2(y_XBYA,basis(j)));
                    
                elseif side_j==modulo(side_i-1,N) && (L(side_j)-bj+ai)<0.15 && (L(side_j)-bj+ai)>5e-15
                    S(i,j) = sum(S_ker_temp.'.*w_XAYBnear.*basis_func_2(x_XAYBnear,basis(i)).*basis_func_2(y_XAYBnear,basis(j)));
                else
                    S(i,j) = integrate2(S_ker_temp,x_gauss,basis(i),basis(j),W);
                end
            end
        end
    end
    progressbar(i/length(basis)) % Update figure
end