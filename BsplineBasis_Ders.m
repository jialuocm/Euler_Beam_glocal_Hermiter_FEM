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
function ders = BsplineBasis_Ders(i, xi, p, U, n)
% Preallocation.
left = zeros(p+1, 1);
right = zeros(p+1, 1);
% Computation of the inverse triangular table starting from degree 0;
% Array structure ndu to store basis function and knot difference;
ndu(1,1) = 1;
% Computate degree 1,2,...,p;
for j = 1:p
    % j meaning degree of basis function;
    left(j)  = xi-U(i+1-j); % Notation of left knote difference;
    right(j) = U(i+j)-xi;   % Notation of right knote difference;
    saved = 0.0;
    for r = 0 : j-1
        % Compute and store Lower triangle
        ndu(j+1,r+1) = right(r+1)+left(j-r);
        temp = ndu(r+1,j) / ndu(j+1,r+1);
        % Compute and store upper triangle
        ndu(r+1,j+1) = saved + right(r+1).*temp;
        saved = left(j-r).*temp;
    end
    ndu(j+1,j+1) = saved;
end
%--------------------------------------------------------------------------
% load the basis function
ders = zeros(1,p+1);
for j = 0:p
    ders(1,j+1) = ndu(j+1,p+1);
end
%Computes the derivatives based on EQ2.10 of Piegl,Tiller book;
for r = 0 : p % Loop over function index;
    s1 = 0; s2 = 1;
    a(1,1) = 1.0;
    % Loop to compute kth derivative of basis;
    for k = 1 : n
        d = 0.0;
        rk = r - k; pk = p - k;
        if r >= k
            a(s2+1,1) = a(s1+1,1) / ndu(pk+1+1,rk+1);
            d = a(s2+1,1) * ndu(rk+1,pk+1);
        end
        if rk >= -1
            j1 = 1;
        else
            j1 = -rk;
        end
        if r-1 <= pk
            j2 = k-1;
        else
            j2 = p -r;
        end
        for j = j1 : j2
            a(s2+1,j+1) = (a(s1+1,j+1) - a(s1+1,j)) / ndu(pk+1+1,rk+j+1);
            d = d + a(s2+1,j+1) * ndu(rk+j+1,pk+1);
        end
        if r <= pk
            a(s2+1,k+1) = -a(s1+1,k) / ndu(pk+1+1,r+1);
            d = d + a(s2+1,k+1) * ndu(r+1,pk+1);
        end
        ders(k+1,r+1) = d;
        j = s1 ;
        s1 = s2;
        s2 = j;  %swich row;
    end
end
%Multiply through by the correct factors
r = p;
for k = 1 : n
    for j = 0 : p
        ders(k+1,j+1) = ders(k+1,j+1) * r;
    end
    r = r*(p-k);
end

end