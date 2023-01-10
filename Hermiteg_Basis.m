% =========================================================================
% This is the shape function routine for one-dimensional finite element code
%
% In this function, adopte the Hermite type basis function.
% Interior nodes are uniformly distributed in the reference domain.
%
% degree: the interpolation basis function is a cubic polynomial, degree three.
% ii    : the number of the basis function. i takes value from 1 to degree+1.
% der   : if der == 0, return the value of the basis function;
%         if der == 1, return the 1st derivative of the basis function;
%         if der == 2, return the 2nd derivative of the basis funciton.
% xa,xb : the point we which to perform the evaluation.
%
% Output: the value of the basis function or the 1st and 2nd derivatives of
%         the basis function.
% -------------------------------------------------------------------------
% By Jia Luo, 2023 Jan. 5th.
% =========================================================================

function poly = Hermiteg_Basis(ii, der, xx, xa, xb)

hh = xb - xa;

switch ii
    case 1
        if der == 0
            poly = -(xx - xb)^2*(-hh + 2*(xa - xx))./hh^3;
        elseif der == 1
            poly = (2*(xb - xx)^2)/hh^3 - ((2*xb - 2*xx)*(hh - 2*xa + 2*xx))/hh^3;
        elseif der == 2
            poly = (-2*(-hh + 2*(xa - xx))+8*(xx - xb))./hh^3;
        end
        
    case 2
        if der == 0
            poly = (xx - xa)*(xx - xb)^2./hh^2;
        elseif der == 1
            poly =(xb - xx)^2/hh^2 + ((2*xb - 2*xx)*(xa - xx))/hh^2;
        elseif der == 2
            poly = 2*(3*xx - 2*xb - xa)./hh^2;
        end

    case 3
        if der == 0
            poly = (xx - xa)^2*(hh + 2*(xb - xx))./hh^3;
        elseif der == 1
            poly = - (2*(xa - xx)^2)/hh^3 - ((2*xa - 2*xx)*(hh + 2*xb - 2*xx))/hh^3;
        elseif der == 2
            poly = (2*(hh + 2*(xb - xx))-8*(xx - xa))./hh^3;
        end

    case 4
        if der == 0
            poly = (xx - xa )^2 * (xx - xb)./hh^2;
        elseif der == 1
            poly = (xa - xx)^2/hh^2 + ((2*xa - 2*xx)*(xb - xx))/hh^2;
        elseif der == 2
            poly =  2*(3*xx - 2*xa - xb)./hh^2;
        end
        
    otherwise
        disp('Hermite polynomial is considered out of the normal range ');

end

% EOF