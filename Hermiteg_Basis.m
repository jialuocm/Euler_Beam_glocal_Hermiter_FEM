% =========================================================================
% This is the shape function routine for one-dimensional finite element code
%
% In this function, we give the Lagrange type basis function over the
% computational domain [-1,1]. Interior nodes are uniformly distributed in
% the reference domain.
%
% degree: the degree of our interpolation basis function space. In our code
%         we give the basis function up to degree six.
% i     : the number of the basis function. i takes value from 1 to degree+1.
% der   : if der == 0, return the value of the basis function;
%         if der == 1, return the 1st derivative of the basis function;
%         if der == 2, return the 2nd derivative of the basis funciton.
% x     : the point we which to perform the evaluation.
%
% Output: the value of the basis function or the 1st and 2nd derivatives of
%         the basis function.
% -------------------------------------------------------------------------
% By Ju Liu, CAM student, 2009 Dec. 24th.
% =========================================================================

function poly = Hermiteg_Basis(i, der, x, xa, xb)

h = xb - xa;

if i == 1
    if der == 0
        poly = -(x-xb)^2*(-h + 2*(xa-x))./h^3;
    elseif der == 1
        poly = 0;
    elseif der == 2
        poly = (-2*(-h+2*(xa-x))+8*(x-xb))./h^3;
    end
end
if i == 2
    if der == 0
        poly = (x - xa)*(x - xb)^2./h^2;
    elseif der == 1
        poly = 0;
    elseif der == 2
        poly = 2*(3*x - 2*xb - xa)./h^2;
    end
end
if i == 3
    if der == 0
        poly = (x-xa)^2*(h + 2*(xb-x))./h^3;
    elseif der == 1
        poly = 0;
    elseif der == 2
        poly = (2*(h+2*(xb-x))-8*(x-xa))./h^3;
    end
end
if i == 4
    if der == 0
        poly = (x - xa )^2 * (x - xb)./h^2;
    elseif der == 1
        poly = 0;
    elseif der == 2
        poly =  2*(3*x - 2*xa - xb)./h^2;
    end
end

% EOF