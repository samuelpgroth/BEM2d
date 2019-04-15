%  geometry.m   A function to define the geometry of the scatterer.                                                           %
%               Samuel Groth   16/12/11

% function[n,s,N,V,tv] = geometry(incident,geom)
function[n,tv,N,V,L] = geometry(incident,geom)
a=0;
theta=0;
if geom == 0           % square
    a = sqrt(2)*pi;
    theta = (pi/4:pi/2:2*pi);
elseif geom == 0.1
    a = sqrt(2);%sqrt(2)*pi;
    theta = (pi/4:pi/2:2*pi);
elseif geom == 0.2
    a = sqrt(2)*2.5;
    theta = (pi/4:pi/2:2*pi);
elseif geom == 1       % triangle
    a=2*pi*sind(30)/sind(120);
%         theta=(pi/12:2*pi/3:2*pi);
    theta=(-pi/6:2*pi/3:7*pi/6);
elseif geom == 1.1
    a=2*sind(30)/sind(120);
%     theta=(-pi/6:2*pi/3:7*pi/6);
%     theta=[pi/2,7*pi/6,11*pi/6];
    theta=[pi/6,5*pi/6,9*pi/6];
%     P=2*sind(30)/sind(120)*(cos(theta)+1i*sin(theta));
 elseif geom==8  % pentagon
    a=2*pi*sin(3*pi/10)/sin(2*pi/5);
%     theta=[pi/2:2*pi/5:pi/2+8*pi/5];  
    theta=[pi/2-2*pi/10:2*pi/5:2*pi];
        
elseif geom ==2       % hexagon
    a=2*pi;
    theta=(pi/3:pi/3:2*pi);
elseif geom == 2.1
    a = 2;
    theta=(pi/3:pi/3:2*pi);
    
end
V = a*[cos(theta);sin(theta)]; V = V';  % define vertices
if geom == 4
    clear V
    %     V = [5 0;5 1;-5 1;-5 0];
    V = [3 0;3 4;-3 4;-3 0];
elseif geom == 5
    V = [5.5 2;7 2.5;7.5 4.5;5 5;3.5 3.5];
elseif geom == 6
    V = [3.141592653589793  -1.813799364234217;0   3.627598728468435;
        -3.141592653589793 3.627598728468435;-3.141592653589793  -1.813799364234216];
elseif geom == 7   % screen
    V= [1 0;-1 0];
    
elseif geom == 11
    V = [1 sqrt(39)/10;-1 sqrt(39)/10;0 -sqrt(39)/10]; 
elseif geom == 12
    V = [1 sqrt(11)/10;-1 sqrt(11)/10;0 -sqrt(11)/10];
    
elseif geom==10
    for i=1:10
        theta=2*pi/10*(i-1);
        if mod(i,2)==1
            V(i,1)=2*cos(theta);
            V(i,2)=2*sin(theta);
        else
            V(i,1)=cos(theta);
            V(i,2)=sin(theta);
        end
    end
%     keyboard
end
N = length(V);                          % number of vertices
% keyboard
% shift = repmat(V(3,:),N,1);
% V = V - shift;      % shifting to make P1 at origin

V(N+1,:) = V(1,:);

s = zeros(N+1,2);     % side vectors (pointing anticlockwise)
n = zeros(N+1,2);     % normals to sides
tv = zeros(N+1,2);    % tangents to sides
L = zeros(N+1,1);     % side lengths
% quiver(-4.5,4.5,incident(1),incident(2),'g','Linewidth',2),hold on;       % Plot incident light vector
figure(1)
for j=1:N
%     plot([V(j,1) V(j+1,1)],[V(j,2) V(j+1,2)],'k','LineWidth',2),axis equal, axis off, hold on; % Plot polygon
    s(j,:) = V(j+1,:)-V(j,:);                     % Side vectors
    L(j) = norm(s(j,:));                          % Side lengths
%     if j>1
%         L(j)=L(j-1);    % slight hack to avoid rounding errors
%     end
    tv(j,:) = s(j,:)/L(j);                        % Normalised side vectors (tangents)
    n(j,:) = [s(j,2),-s(j,1)];                    % Outward normals of sides
    n(j,:) = n(j,:)/norm(n(j,:));                 % Normalise
end
s(N+1,:)=s(1,:);
n(N+1,:)=n(1,:);
L(N+1)=L(1);
tv(N+1,:) = tv(1,:);