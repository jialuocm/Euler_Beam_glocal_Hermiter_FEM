% simulation the beam with Hermite polynomial
clear all; clc;clear;close all;

EI = 2 ; % bend module stiffness

omega_L = 0; % phsical length of beam
omega_R = 1; % position length of beam
L = omega_R - omega_L; %  phsical length of beam

%q = @(x) x./((omega_R - omega_L))*(1-x/(omega_R - omega_L));
q = @(x) -48;    %constant distributive load
g_f = @(x) 0;  %essential BC transverse displacement
dg_f = @(x) 0; %essential BC derivative of disp slope;
%h_f = @(x) 1;% Natual BC prescribed Moment
%dh_f = @(x) -2;% Natual BC prescribed shear force

% exact solution

U_fx = @(x) (q(x)/24.*x.^4 - 1/6.*q(x).*L.^3.*x + 1/8.*q(x).*L.^4)./EI;

fai_fx = @(x) (1/6.*q(x).*x.^3 -  1/6.*q(x).*L.^3)./EI;

U_f2x = @(x) (1/2.*q(x).*x.^2)./EI;

M_fx = @(x) 1/2.*q(x).*x.^2;

%-------------------------------------------------------------------

nElem_x = 3: 6: 36;

%nElem_x = 12;

cn = length(nElem_x);


hh_x = zeros(cn,1);
error_l2_x = zeros(cn,1);
error_h2_x = zeros(cn,1);

hh_lg = zeros(cn,1);
slop_l2 = zeros(cn,1);
slop_h2 = zeros(cn,1);

for num = 1 : cn 

%Hermite-cubic-interpolation degree
pp = 3;

%number of elements
nElem = nElem_x(num);

n_en = 2; % node number of element
n_np = nElem + 1;
d_node = 2; % node degree of freedom: displacements and slopes
d_en = n_en*2; % element DOF 
nLocBas = d_en;% local Basis equal to DOF of element
nFunc = nLocBas+2*(nElem-1);

%assemble
IEN  = zeros(d_en, nElem);
for ee = 1 : nElem
    for aa = 1 : d_en
        IEN (aa, ee) = (ee -1 ) * d_node + aa;
    end
end

%mesh

hh_x(num) = (omega_R - omega_L) / nElem;

x_coor = omega_L : hh_x(num) : omega_R;

%setup ID array based on BC

ID = 1 : nFunc;
ID(end-1:end) = -1; %  the right of beam is fix, transverse dis and slope equal 0
n_eq = nFunc - 2; % Dirichlet nodes is kown
T_IEN = table(IEN,'variableNames',{'IEN'});
writetable(T_IEN);
%allocate an empty stiffness matrix and load vector

K = sparse(nFunc, nFunc);
F = zeros(nFunc,1);
matrix_Ksize = size(K);
matrix_Fsize = size(F);
%quadrature rule
nqp = pp + 1;
%assembly the stiffness matrix and load vector
for ee = 1: nElem
    %allocate zero element stiffness matrix and element load vector
    % k_ele = zeros(d_en, d_en);
    % f_ele = zeros(d_en,1);
    x_ele = zeros(d_en,1);
    for aa = 1 : d_en

        if mod(IEN(aa,ee),2)==0 % check if node number is even
            x_ele(aa) = x_coor(IEN(aa,ee)/2);
        else
            x_ele(aa) = x_coor((IEN(aa,ee)+1)/2);
        end

    end
    %disp(ee);

    [qp, wq] = Gauss(nqp, x_ele(1), x_ele(end));

    %disp(wq);

    for qua = 1 : nqp
    
        % x_qua = 0;
        dx_dxi = 0;
        d2x_dxi = 0;
        e_l = x_ele(end)- x_ele(1);

        for aa = 1 : nLocBas
            Na = Hermiteg_Basis(pp,aa,0,qp(qua),e_l,x_ele(1),x_ele(end));
            Na_2x = Hermiteg_Basis(pp,aa,2,qp(qua),e_l,x_ele(1),x_ele(end));
            %calculate the global K matrix and force vector
            AA = ID(IEN(aa,ee));
            if( AA > 0 )
                F(AA) = F(AA) + wq (qua) * Na *q(qp(aa));
                for bb = 1: nLocBas
                    BB = ID(IEN(bb,ee));
                    if (BB > 0)
                        %  Nb_xi = HPolyBasis(pp, bb, 1, qp(qua),e_l);
                        Nb_2x = Hermiteg_Basis(pp,bb,2,qp(qua),e_l,x_ele(1),x_ele(end));
                        K(AA,BB) = K(AA,BB) + wq(qua) * EI* Na_2x * Nb_2x;
                    end

                end
           % disp(K);
           

            else
                if mod(IEN(aa,ee)+1,2)==0 % check if node number is even
                    x_p = x_coor((IEN(aa,ee)+1)/2);
                    g_ubc = g_f(x_p);   % esstential BC : transverse displace
                    % g_uxbc = dg_f(x_p); % esstential BC : slope
                    % AA = IEN(aa,ee);
                    AA =  IEN(aa,ee);
                    %disp(AA);
                    K( AA, AA ) = 1.0;
                    F(AA) = F(AA) + g_ubc;

                else
                    x_p = x_coor(IEN(aa,ee)/2);
                    % g_ubc = g_f(x_p);   % esstential BC : transverse displace
                    g_uxbc = dg_f(x_p);   % esstential BC : slope
                    %AA = IEN(aa,ee);
                    AA =  IEN(aa,ee);
                    %disp(AA);
                    K( IEN(aa,ee), IEN(aa,ee) ) = 1.0;
                    F(AA) = F(AA) + g_uxbc;
                end
            end

        end
        %End up quadrature loop

    end

    % modified the load vector by natural BC
    %  if ee == 1
    %  F(ID(IEN(1,ee))) =  F(ID(IEN(1,ee))) + dh_f(x_coor(IEN(1,ee)));% prescribed shear force
    %  F(ID(IEN(2,ee))) =  F(ID(IEN(2,ee))) + h_f(x_coor(IEN(2,ee)));% prescribed Moment force
    %end
end

%disp(K);disp(F);
%solve the stiffness matrix problem

 Uh = K \ F;

 uh_od = Uh(1:2:nFunc);
 fai_even = Uh(2:2:nFunc);

 figure
 uhp = plot(x_coor,uh_od,'--r*','linewidth',2);
 hold on
 u_exact = U_fx(x_coor);
 up = plot(x_coor, u_exact,'b--','LineWidth',2);

 legend([uhp,up],'uh FEM','u Exact');
 xlabel('x coor');
 ylabel('u');
 hold off
 exportgraphics(gca,['file uh' '.jpg']);


 figure
 faih_p = plot(x_coor,fai_even,'--r*','linewidth',2);

 hold on
 p_M_exact = plot(x_coor,fai_fx(x_coor),'b--','LineWidth',2);
 legend([faih_p,p_M_exact],'\phi_h FEM','\phi_{Exact}');
 xlabel('x coor');
 ylabel('\phi ');
 hold off
 exportgraphics(gca,['file_fai' '.jpg']);


% error convergence analysis
error_l2  = 0.0;
error_H2  = 0.0;

j=1;
M_h = zeros(nElem*nqp,1);
M_exact = zeros(nElem*nqp,1);
x_m = zeros(nElem*nqp,1);
x_b = zeros(nElem+1,1);
U_h_p = zeros(nElem*nqp,1);

 for ee = 1: nElem
     u_ele = zeros(nLocBas,1);
     x_ele = zeros(nLocBas,1);

     for aa = 1 : nLocBas

         u_ele(aa) = Uh(IEN(aa,ee));

         if mod(IEN(aa,ee)+1,2)==0 % check if node number is even
            
             x_ele(aa) = x_coor((IEN(aa,ee)+1)/2);
         else
            
             x_ele(aa) = x_coor(IEN(aa,ee)/2);
         end

     end
 
     e_l = x_ele(end)- x_ele(1);

     [qp, wq] = Gauss(nqp, x_ele(1), x_ele(end));
           %disp(qp);
           qp = qp(end:-1:1);
           %disp(qp);

     for qua = 1 : nqp

         U_h    = 0.0;
         u_dx2  = 0.0 ;
         m_dx2  = 0.0 ;

         for aa = 1 : nLocBas

             U_h    = U_h   + u_ele(aa) * Hermiteg_Basis(pp,aa,0,qp(qua),e_l,x_ele(1),x_ele(end));
             u_dx2  = u_dx2 + u_ele(aa) * Hermiteg_Basis(pp,aa,2,qp(qua),e_l,x_ele(1),x_ele(end));
             m_dx2  = m_dx2 + u_ele(aa) * Hermiteg_Basis(pp,aa,2,x_ele(1),e_l,x_ele(1),x_ele(end));

         end

         x_m(j)   =  qp(qua);
         U_h_p(j) =  U_h;

         u_exact    = U_fx(qp(qua));
         u_dx2_exact = U_f2x(qp(qua));

         x_b(j)   =  x_ele(1);
         M_h(j)   =  EI.* m_dx2;
         M_exact(j) = M_fx(x_ele(1));

         j = j+1 ;

         errorl2 = error_l2 + wq(qua) * (U_h - u_exact) * (U_h - u_exact);
         errorH2 = error_H2 + wq(qua) * (u_dx2 - u_dx2_exact) * (u_dx2 - u_dx2_exact);

     end
 end

  error_l2_x(num,1)  = log(sqrt(errorl2));
  error_h2_x(num,1)  = log(sqrt(errorH2));
 hh_lg(num) = log(hh_x(num))./L;
 % hh_lg(num) = log(hh_x(num));


     dx_h = zeros( 2,1);
     dy_e = zeros( 2,1);
     if (num >= 2)
         dx_h = [hh_lg(num);hh_lg(num-1)];

         dy_l2e = [error_l2_x(num);error_l2_x(num-1)];
         dy_h2e = [error_h2_x(num);error_h2_x(num-1)];

         slop_l2(num) = (dy_l2e(2) - dy_l2e(1)) / (dx_h(2)-dx_h(1));
         slop_h2(num) = (dy_h2e(2) - dy_h2e(1)) / (dx_h(2)-dx_h(1));

     else

         slop_l2(num) = 0;
         slop_h2(num) = 0;
     end

 figure

 p_M_h = plot(x_b,M_h,'--rO','LineWidth',2);
 hold on
 p_M_exact = plot(x_b,M_exact,'--b*','LineWidth',2);
 legend([p_M_h ,p_M_exact],'M_h FEM','M_{Exact}');
 xlabel('x coor');
 ylabel('M )');
 hold off
 exportgraphics(gca,['file_M' num2str(num) '.jpg']);

end


 figure
 error_h = plot(hh_lg,error_l2_x,'--rO','LineWidth',2);

 legend(error_h,'L_2 error convergence analysis');
 xlabel('log ||hh|| ');
 ylabel('log ||e||_l2 ');
 exportgraphics(gca,['error_u_l2' '.jpg']);

 figure
 error_h = plot(hh_lg,error_h2_x,'--rO','LineWidth',2);

 legend(error_h,'H_2error convergence analysis');
 xlabel('log ||hh|| ');
 ylabel('log ||e||_H2 ');
 exportgraphics(gca,['error_u_h2' '.jpg']);

T = table(hh_x,error_l2_x,error_h2_x,slop_l2,slop_h2,'variableNames',{'hh_mesh','error_l2','error_H2','L2 convergence rate','H2 convergence rate'});
writetable(T);
T

% END 