% Compute the nonvanishing basis functions.
% Input:
%   i: index of the basis function to compute (this value cannot be smaller
%       than p);
%   xi: value in which the basis function is being evaluated;
%   p: degree of the basis function;
%   Xi: knot vector over which the basis function is being built.
% Output:
%   N: vector containing the value of all the nonvanishing basis functions:
%   N(1)=N_{i-p},...,N(p+1)=N_{i}.

%U  = [0,0,0,0,2,4,4,6,6,6,8,8,8,8];
U  = [0,0,0,1,2,3,4,4,5,5,5];
p = 2;
%xi = 2.5;
xi_vector = linspace(1,2,150);
nm = length(xi_vector);
M_data = zeros(p+1,nm);

for nn = 1:nm
    xi = xi_vector(nn);
    i  = FindSpan(p,U,xi);
    Ni_p = Bspline_Basis(i,xi,p,U);
    M_data(:,nn) = Ni_p;
end
plot(xi_vector,M_data(1,:),'--rO','LineWidth',2);
hold on;
plot(xi_vector,M_data(2,:),'--b*','LineWidth',2);
plot(xi_vector,M_data(3,:),'--k*','LineWidth',2);
hold off;

% =========================================================================
function N = Bspline_Basis(ii, xi, pp, U)
    % Preallocation.
    N = zeros(pp+1, 1);
    left = zeros(pp+1, 1);
    right = zeros(pp+1, 1);
    % Computation of the inverse triangular table starting from degree 0;
    N(1) = 1; % Ni,0, degree equal 0 basis function;
    % computate degree 1,2,...,p;
    for j = 1:pp
        % j meaning degree of basis function;
        left(j)  = xi-U(ii+1-j); % notation of left knote difference;
        right(j) = U(ii+j)-xi;   % notation of right knote difference; 
        saved = 0.0;
        for r = 0 : j-1
            temp = N(r+1)./(right(r+1)+left(j-r));
            N(r+1) = saved + right(r+1).*temp;
            saved = left(j-r).*temp;
        end
        N(j+1) = saved;
    end
end

% =========================================================================
function mid = FindSpan(pp,U,xi)
% returns the knot span index
n = length(U)-pp;
if xi==U(end)
    % special case
    mid = n-1;
    return
end
low = pp;
high = n;
mid = floor((low+high)/2);
while xi<U(mid) || xi>= U(mid+1)
    if xi<U(mid)
        high = mid;
    else
        low = mid;
    end
    mid = floor((low+high)/2);
end
end
% =========================================================================
%END






