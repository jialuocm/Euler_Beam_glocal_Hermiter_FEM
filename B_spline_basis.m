% =========================================================================
% This is the B-spline function routine for finite element code
% Input: knot vector, order to construct a B-spline
% degree: the order of basis function.
% knot vector: open knot vector.
% der   : if der == 0, return the value of the basis function;
%         if der == 1, return the 1st derivative of the basis function;
%         if der == 2, return the 2nd derivative of the basis funciton.
% xa,xb : the point we which to perform the evaluation.
%
% Output: 
% Nip B : B-spline basis functions
% dNip1 : 1st derivative of B-spline;
% dNip2 : 2st derivative of B-spline;
% -------------------------------------------------------------------------
% By Jia Luo, 2023 Jan. 16th.
% =========================================================================

function [Nip,dNip1,dNip2] = B_spline_basis(knot,pp,xx)
%-------------------------------------------------------------

n = length(knot)-pp-1; % n: number of basis functions

%calculate B-spline basis functions at a point==

Nip = zeros(length(knot)-1,pp+1);

for i = 1:length(knot)-1 % 0 oder calculation
    if xx >= knot(i) && xx < knot(i+1)
        Nip(i,1) = 1;
    end
end

% calculaate B-spline basis by Cox-de Boor recursion formula 
for k = 2:pp+1 % Higher-oder calculation, k=2, 1st-oder; k=3 2nd-oder;
    for i = 1:(length(knot)-k)
        pod = k - 1; % polynomial oder
        c1 = knot(i+pod) - knot(i);
        if c1 ~= 0
            c1 = (xx - knot(i))/c1;
        end
        c2 = knot(i+pod+1) - knot(i+1);
        if c2 ~= 0
            c2 = (knot(i+pod+1)-xx)/c2;
        end  
        Nip(i,k) = c1 * Nip(i,k-1) + c2 * Nip(i+1,k-1);
    end
end
Nip = Nip(i:n,pp+1);
%== compute 1st derivatives==

dNip1 = zeros(n,1);
for i = 1:n
    c1 = knot(i+pod) - knot(i);
    if c1 ~= 0
        c1 = pp/c1;
    end
    c2 = knot(i+pod+1) - knot(i+1);
    if c2 ~= 0
        c2 = pp/c2;
    end

    dNip1(i) = c1 * Nip(i,pp-1) - c2 * Nip(i+1,pp-1);
end

%== compute 2st derivatives==

dNip2 = zeros(n,1);
for i = 1:n
    c1 = (knot(i+pp) - knot(i))*(knot(i+pp-1) - knot(i));
    if c1 ~= 0
        c1 = pp*(pp-1)/c1;
    end
    c2a = (knot(i+pp) - knot(i))*(knot(i+pp) - knot(i));
    if c2a ~= 0
        c2a = pp*(pp-1)/c2a;
    end

    c2b = (knot(i+pp+1) - knot(i+1))*(knot(i+pp) - knot(i+1));
    if c2b ~= 0
        c2b = pp*(pp-1)/c2b;
    end 

    c2 = c2a + c2b;

    c3 = (knot(i+pp+1) - knot(i+1)) * (knot(i+pp+1) - knot(i+2));
    if c3 ~= 0
        c3 = pp*(pp-1)/c3;
    end 
    dNip2(i) = c1 * Nip(i,pp-2) - c2 * Nip(i+1,pp-2) + c3 * Nip(i+2,pp-2);
end

% EOF