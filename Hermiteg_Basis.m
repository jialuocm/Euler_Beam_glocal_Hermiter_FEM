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

function poly = Hermiteg_Basis(ii, der, xx, xa, xb)

h = xb - xa;

if ii == 1
    if der == 0
        poly = -(xx-xb)^2*(-h + 2*(xa-xx))./h^3;
    elseif der == 1
        poly = 0;
    elseif der == 2
        poly = (-2*(-h+2*(xa-xx))+8*(xx-xb))./h^3;
    end
end
if ii == 2
    if der == 0
        poly = (xx - xa)*(xx - xb)^2./h^2;
    elseif der == 1
        poly = 0;
    elseif der == 2
        poly = 2*(3*xx - 2*xb - xa)./h^2;
    end
end
if ii == 3
    if der == 0
        poly = (xx-xa)^2*(h + 2*(xb-xx))./h^3;
    elseif der == 1
        poly = 0;
    elseif der == 2
        poly = (2*(h+2*(xb-xx))-8*(xx-xa))./h^3;
    end
end
if ii == 4
    if der == 0
        poly = (xx - xa )^2 * (xx - xb)./h^2;
    elseif der == 1
        poly = 0;
    elseif der == 2
        poly =  2*(3*xx - 2*xa - xb)./h^2;
    end
end

% EOF