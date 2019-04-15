% y_minus_x_calc.m

% This function converts Gaussian quadrature abscissae on the interval
% [0,h] into coordinates of Y-X
% Note: absc = abscissae


function[y_minus_x_coords] = y_minus_x_calc(x,y,tv,side_x,side_y,N)

% if side_y==modulo(side_x-1,N)   % elements touch at corner (xa touches yb)
    y_minus_x_coords = y*tv(side_y,:) - x*tv(side_x,:);
%     keyboard
% else
%     y_minus_x_coords = y*tv(side_y,:) + x*tv(side_x,:);
%     keyboard
% end
   
end