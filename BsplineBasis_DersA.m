% Compute the nonvanishing basis functions.
% Input:
%   i: index of the basis function to compute (this value cannot be smaller
%       than p+1);
%   xi: value in which the basis function is being evaluated;
%   pp: degree of the basis function;
%   U: knot vector over which the basis function is being built.
%   k: defined accordingly to the knot vector;
% Output:
%   N: vector containing the value of all the nonvanishing basis functions:
%   N(1)=N_{i-p},...,N(p+1)=N_{i}.
%   ders: bidimensional matrix where the element in position
%         (k,j) is the kth derivative of the function
%         N_{i-p+j,p} with 0<=k<=n and 0<=j<=p.
% =========================================================================
%U  = [0,0,0,0,1,3,4,6,6,6,8,8,8,8];
 U  =  [0,0,0,0,1/4,1/2,3/4,1,1,1,1];
%U  =  [0,0,1/2,1,1];
%U  =  [0,0,0,1/2,1,1,1];
pp = 3;
k  = 2;
% xi = 2.5;
xi_vector = linspace(0,1,500);
nm = length(xi_vector);
M_data = zeros(pp+1,nm);
Ders_data = zeros(pp+1,nm);

for nn = 1:nm
    xi = xi_vector(nn);
    i  = FindSpan(pp,U,xi);
    % Ni_p = Bspline_Basis(i,xi,pp,U);
    Ni_p_k = BsplineBasis_Ders(i,xi,pp,U,k);
    % M_data(:,nn) = Ni_p;
    M_data(:,nn) = Ni_p_k(1,:);
    Ders_data(:,nn) = Ni_p_k(2,:);
end
% =========================================================================
%Postprocessing:plotting to check validation of basis and deriv of functions;
close all;

plot(xi_vector,M_data(1,:),'-r','LineWidth',2);
hold on;
plot(xi_vector,M_data(2,:),'-b','LineWidth',2);
plot(xi_vector,M_data(3,:),'-k','LineWidth',2);
plot(xi_vector,M_data(4,:),'-g','LineWidth',2);
hold off;

figure
hold on;
plot(xi_vector,Ders_data(1,:),'-r','LineWidth',2);
plot(xi_vector,Ders_data(2,:),'-b','LineWidth',2);
plot(xi_vector,Ders_data(3,:),'-k','LineWidth',2);
plot(xi_vector,Ders_data(4,:),'-g','LineWidth',2);
hold off;
% =========================================================================

function ders = BsplineBasis_Ders(i, xi, p, U,n)
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

% =========================================================================
function mid = FindSpan(p,U,xi)
% returns the knot span index
n = length(U)-p;
if xi==U(end)
    % special case
    mid = n-1;
    return
end
low = p;
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