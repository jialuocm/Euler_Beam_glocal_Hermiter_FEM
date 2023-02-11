% =========================================================================
% This is the FEM Matlab main code for one-dimensional Bernulli-Euler beam.
% In this main, adopte the B-spine basis function.
% Input:
%   i: index of the basis function to compute (this value cannot be smaller
%       than p+1);
%   xi: value in which the basis function is being evaluated;
%   pp: degree of the basis function;
%   U: knot vector over which the basis function is being built.
%   k: defined accordingly to the knot vector;
% -------------------------------------------------------------------------
% By Jia Luo, 2023 Jan. 11th.
% =========================================================================
clear all; clc; close all;

pi  = atan(1) * 4;

EI = 1 ;                         % Bend module stiffness
%M  = 3 ;                         % Natual BC prescribed Moment
M  = 0 ;                         % Natual BC prescribed Moment
%Q  = -12 ;                       % Natual BC prescribed shear force
Q  = 0 ;                       % Natual BC prescribed shear force

omega_L = 0;                     % Phsical location of beam
omega_R = pi;                    % Phsical location of beam
LL = omega_R - omega_L;          % Phsical length of beam

f_q = @(x) sin(x);               % Distribution body force

u_g = @(x) 0;                    % Essential BC transverse displacement
u_x_g = @(x) 0;                  % Essential BC derivative of disp slope

M_f = @(x) M / EI;               % Natual BC prescribed Moment
Q_f = @(x) Q / EI;               % Natual BC prescribed shear force
%-------------------------------------------------------------------
% Exact solution
c1 = Q + 1.0;
c2 = M;
c3 = 1.0 - 0.5 * (Q + 1.0) * pi^2 - M * pi ;
c4 = (Q + 1.0) * pi^3 / 3 - pi + 0.5 * M * pi^2;

u_d0x = @(x) (sin(x) + c1.*x.^3 / 6 + 0.5 * c2.*x.^2 + c3.*x + c4)./EI; % Displacements
u_d1x = @(x) (cos(x) + 0.5 * c1 * x.^2 + c2 * x + c3)./EI;            % First derivative
u_d2x = @(x) (-sin(x) + c1 .* x + c2) ./ EI;                          % Second derivative
%-------------------------------------------------------------------

nElem_x = 3 : 1: 9;
%nElem_x = 3;

cn         = numel(nElem_x);
hh_x       = zeros(cn,1);
error_l2_x = zeros(cn,1);
error_h1_x = zeros(cn,1);
error_h2_x = zeros(cn,1);

hh_lg   = zeros(cn,1);
slop_l2 = zeros(cn,1);
slop_h1 = zeros(cn,1);
slop_h2 = zeros(cn,1);

for num = 1 : cn
    nElem   = nElem_x(num);   % Number of elements
    disp(nElem);
    hh_x(num) = (omega_R - omega_L) / nElem;
    %U  = [0,0,0,0,2,4,4,6,6,6,8,8,8,8].*pi/8; % representing mesh of problem;
    %U  = [0,0,0,0,1/3,2/3,1,1,1,1].*pi;
    pp = 3;
    k  = 2;
    nqp = pp + 1;  % Quadrature rule
    U = zeros(1,2 * pp + nElem + 1);
    U(1,end-pp:end) = omega_R;

    for ii = 1 : nElem - 1
        U(pp+1+ii) = ii * hh_x(num);
    end
    disp(U);
    
    x_coor = zeros(1,nElem + 1);
    I_u = length(U) ;
    N_i = 1; 
    x_coor(N_i) = U(pp+1);
    for ii = pp+2 : I_u
        x_a = x_coor(N_i); x_b = U(ii);
        if x_b > x_a
            N_i = N_i+1;
            x_coor(N_i) = x_b;
        end
    end

    disp(x_coor);

    n_np    = length(x_coor);     % Number of global node
    nElem   = length(x_coor) - 1; % Number of elements

    n_en    = 2;   % Number of local node of element
    n_dof   = 2;   % Degree of freedom number per node, disp and slopes

    nLocBas = n_en * n_dof; % Local equation number
    nFunc   = length(U)-pp-1; % n_bs: number of global basis functions;
    disp(nFunc);

    % Assemble IEN
    IEN  = zeros(n_en, nElem);
    for ee = 1 : nElem
        for aa = 1 : n_en
            IEN (aa, ee) = ee - 1 + aa;
        end
    end
    disp(IEN);
    % Mesh of geometric domain
    % Setup ID array and modify the ID array based on BC
    ID = zeros(n_dof,n_np);
    for g_np = 1: n_np
        for ii = 1 : n_dof
            ID(ii,g_np) = (g_np-1) * n_dof + ii;
        end
    end
    % Assign the ID for the Dirichlet node to be -1
    ID(:,n_np) = -1;
    disp(ID);
    n_eq = n_np * n_dof - 2;
    % LM matrix
    LM = zeros(nLocBas,nElem);
    for ee = 1 : nElem
        for aa = 1 : n_en
            A = IEN(aa,ee);
            for ii = 1 : n_dof
                n_eq_i = n_dof * (aa - 1) + ii;
                LM(n_eq_i,ee) = ID(ii,A) ;
            end
        end
    end
    disp(LM);
    % Allocate an empty stiffness matrix and load vector
    K = sparse(n_eq, n_eq);
    F = zeros(n_eq, 1);

    %-------------------------------------------------------------------
    % Assembly the stiffness matrix and load vector
    for ee = 1 : nElem
        % Allocate zero element stiffness matrix and element load vector
        x_ele  = zeros(n_en,1);
        u_ebc  = zeros(nLocBas,1);
        % New data structure: mapping local element information to global
        for aa = 1 : n_en
            x_ele(aa) = x_coor( IEN(aa,ee) ) ; % Physical geometrical position
        end

        [qp, wq] = Gauss(nqp, x_ele(1), x_ele(end));
        % Find index of Pysical mesh in index space;
        i_id  = FindSpan(pp,U,x_ele(1));

        for qua = 1 : nqp

            Ni_p_k = BsplineBasis_Ders(i_id, qp(qua), pp, U, k);
            for aa = 1 : nLocBas
                Na    = Ni_p_k(1,aa);
              % Na_x  = Ni_p_k(2,aa);
                Na_2x = Ni_p_k(3,aa);
                % Calculate the global K matrix and force vector
                AA = LM(aa,ee);
                if( AA > 0 )
                    F(AA) = F(AA) + wq (qua) * Na * f_q(qp(qua));
                    for bb = 1 : nLocBas
                        BB = LM(bb,ee);
                        if ( BB > 0 )
                            Nb_2x = Ni_p_k(3,bb);
                            K(AA,BB) = K(AA,BB) + wq(qua) * EI * Na_2x * Nb_2x;
                        else

                            Nb_2x = Ni_p_k(3,bb);
                            % Obtain the Dirichlet node's physical coordinates
                            if mod(bb,2)==0 % calculate Local dof i;
                                i = 2;
                            else
                                i = 1;
                            end
                            a_i = (bb - i) / n_dof + 1;
                            x_g_u = x_coor(IEN(a_i,ee)) ;
                            % Obtain the boundary data at this point
                            % Esstential BC: Node displacement; Node slope
                            u_ebc = [0; 0; u_g(x_g_u); u_x_g(x_g_u) ];

                            F( AA ) = F( AA ) -  wq(qua) * EI * Na_2x * Nb_2x * u_ebc(bb);

                        end
                    end
                end
            end

        end % End of quadrature loop

        % Modified the load vector by natural BC
        if ee == 1
            x_coor_q = x_coor(IEN(1,ee));
            x_coor_m = x_coor(IEN(1,ee));
            Ni_p_k = BsplineBasis_Ders(i_id,x_coor_q, pp, U, k);
            Na    = Ni_p_k(1,aa);
            Na_x  = Ni_p_k(2,aa);
            F(ID(1,IEN(1,ee))) =  F(ID(1,IEN(1,ee))) + Na * Q_f(x_coor_q);  % prescribed shear force
            F(ID(2,IEN(1,ee))) =  F(ID(2,IEN(1,ee))) - Na_x * M_f(x_coor_m); % prescribed Moment force
        end

    end

    % Solve the stiffness matrix problem

    d_Ah = K \ F;
    % Append the displacement vector by the Dirichlet data
    d_Ah = [ d_Ah; u_g(x_g_u); u_x_g(x_g_u) ];
    disp(d_Ah);
    % =========================================================================
    % Displacement distribution
   
    %     uh_od    = Uh(1:2:nFunc - 1);
    %     fai_even = Uh(2:2:nFunc);

    ID = zeros(n_dof,n_np);
    for g_np = 1: n_np
        for ii = 1 : n_dof
            ID(ii,g_np) = (g_np-1) * n_dof + ii;
        end
    end
    disp(ID);

    for ee = 1:nElem
        x_ele = zeros(n_en,1);
        d_ele = zeros(nLocBas,1);
        %d_u_ele   = zeros(n_en,1);
       % d_fai_ele = zeros(n_en,1);

        for aa = 1 : n_en
            x_ele(aa)     = x_coor(IEN(aa,ee)) ;
            for ii = 1 : n_dof
                n_locBas_i = n_dof * (aa - 1) + ii;
                d_ele(n_locBas_i) = d_Ah(ID(ii,IEN(aa,ee)));
            end
%           d_u_ele(aa)   = d_Ah(ID(1,IEN(aa,ee)));
%           d_fai_ele(aa) = d_Ah(ID(2,IEN(aa,ee)));
        end
        % disp(d_ele);
        npl = 1;
        %x_coor_p = linspace( x_ele(1), x_ele(end), npl);
         x_coor_p = 0.5 * (x_ele(1) + x_ele(end));
        %disp(x_coor_p);

        for x_pl = 1 : npl
            U_h     = 0.0;
            U_hdx1  = 0.0;
            Ni_p_k  = BsplineBasis_Ders(i_id, x_coor_p(x_pl), pp, U, k);

            for aa = 1 : n_en
                for ii = 1: n_dof
                    nLocBas_j = n_dof * (aa - 1) + ii;
                    Na     = Ni_p_k(1, nLocBas_j);
                    Na_1x  = Ni_p_k(2, nLocBas_j);
                    Na_2x  = Ni_p_k(3, nLocBas_j);
                    U_h    = U_h + d_ele(aa) * Na;
                    U_hdx1 = U_hdx1 + d_ele(aa) * Na_1x;
                end

            end
            j = (ee - 1 ) * npl + x_pl;
            x_sample(num, j) = x_coor_p(x_pl);
            U_hp(num, j)     = U_h;
            U_hdx1p(num, j)  = U_hdx1;
            u_exact(num, j)  = u_d0x(x_coor_p(x_pl));
            u_dx1exact(num, j)  = u_d1x(x_coor_p(x_pl));
        end
    end
    disp(x_sample); disp(U_hp); disp(u_exact)

    % =========================================================================
        figure;
        uh_p = plot(x_sample, U_hp(num,:),'--r*','linewidth',2);
        hold on
        u_p = plot(x_sample, u_exact(num,:),'--bO','LineWidth',2);
    
        %legend([uh_p,u_p],'uh FEM','u Exact');
        xlabel('x coor');
        ylabel('u');
        hold off
        exportgraphics(gca,['file uh' '.jpg']);

    %------------------------------------------------------------------ 
   % slope distribution
        figure
        faih_p = plot(x_sample,U_hdx1p,'--r*','linewidth',2);
    
        hold on
        p_M_exact = plot(x_sample, u_dx1exact,'--bO','LineWidth',2);
  %     legend([faih_p,p_M_exact],'\phi_h FEM','\phi_{Exact}');
        xlabel('x coor');
        ylabel('\phi ');
        hold off
        exportgraphics(gca,['file_fai' '.jpg']);

    % =========================================================================
    % Now to postprocessing, Error convergence analysis
    nqp = 5 ;
    error_l2  = 0.0;
    error_H1  = 0.0;
    error_H2  = 0.0;

    N_i=1;
    M_h     = zeros(nElem*nqp,1);
    M_exact = zeros(nElem*nqp,1);
    x_m     = zeros(nElem*nqp,1);
    x_b     = zeros(nElem+1,1);
    U_h_p   = zeros(nElem*nqp,1);


    for ee = 1: nElem
        x_ele = zeros(n_en,1);
        d_ele = zeros(nLocBas,1);
        d_u_ele = zeros(n_en,1);
        d_fai_ele = zeros(n_en,1);

        for aa = 1 : n_en
            x_ele(aa) = x_coor(IEN(aa,ee)) ;
            for ii = 1 : n_dof
                n_locBas_i = n_dof * (aa - 1) + ii;
                d_ele(n_locBas_i) = d_Ah(ID(ii,IEN(aa,ee)));
            end                 
        end
        disp(x_ele);
        [qp, wq] = Gauss(nqp, x_ele(1), x_ele(end));
        % disp(qp);
        % Find index of Pysical mesh in index space;
        i_id  = FindSpan(pp, U, x_ele(1));

        for qua = 1 : nqp
            U_h     = 0.0;
            U_hdx1  = 0.0;
            U_hdx2  = 0.0;

            Ni_p_k = BsplineBasis_Ders(i_id, qp(qua), pp, U, k);

            for aa = 1 : n_en
                for ii = 1 : n_dof
                    nLocBas_i = n_dof * (aa - 1) + ii;
                    Na     = Ni_p_k(1, nLocBas_i);
                    Na_1x  = Ni_p_k(2, nLocBas_i);
                    Na_2x  = Ni_p_k(3, nLocBas_i);

                    U_h    = U_h + d_ele(aa) * Na;
                    U_hdx1 = U_hdx1 + d_ele(aa) * Na_1x;
                    U_hdx2 = U_hdx2 + d_ele(aa) * Na_2x;
                   
                end

            end

            u_exact     = u_d0x(qp(qua));
            u_dx1_exact = u_d1x(qp(qua));
            u_dx2_exact = u_d2x(qp(qua));

            error_l2 = error_l2 + wq(qua) * (U_h - u_exact) * (U_h - u_exact);
            error_H1 = error_H1 + wq(qua) * (U_hdx1 - u_dx1_exact) * (U_hdx1 - u_dx1_exact);
            error_H2 = error_H2 + wq(qua) * (U_hdx2 - u_dx2_exact) * (U_hdx2 - u_dx2_exact);

        end
    end

    error_l2 = sqrt(error_l2);
    error_H1 = sqrt(error_H1);
    error_H2 = sqrt(error_H2);

    error_l2_x(num,1)  = log(error_l2);
    error_h1_x(num,1)  = log(error_H1);
    error_h2_x(num,1)  = log(error_H2);
    hh_lg(num) = log(hh_x(num) / LL); % Normlized the mesh length size

    dx_hh  = zeros( 2,1);
    dy_l2e = zeros( 2,1);
    dy_h1e = zeros( 2,1);
    dy_h2e = zeros( 2,1);
    if (num >= 2)
        dx_hh  = [hh_lg(num); hh_lg(num-1)];

        dy_l2e = [error_l2_x(num); error_l2_x(num-1)];
        dy_h1e = [error_h1_x(num); error_h1_x(num-1)];
        dy_h2e = [error_h2_x(num); error_h2_x(num-1)];

        slop_l2(num) = (dy_l2e(2) - dy_l2e(1)) / (dx_hh(2)-dx_hh(1));
        slop_h1(num) = (dy_h1e(2) - dy_h1e(1)) / (dx_hh(2)-dx_hh(1));
        slop_h2(num) = (dy_h2e(2) - dy_h2e(1)) / (dx_hh(2)-dx_hh(1));

    else
        slop_l2(num) = 0;
        slop_h1(num) = 0;
        slop_h2(num) = 0;
    end

end
% =========================================================================
% Postprocessing
figure
yyaxis left
error_h_l2 = plot(hh_lg,error_l2_x,'--rO','LineWidth',2);
%ylim([-18 -2]);

legend(error_h_l2,'L_2 error convergence analysis');
xlabel('log ||hh||/L ');
ylabel('log ||e||_{l2} ');
% exportgraphics(gca,['error_u_l2' '.jpg']);
hold on;

yyaxis right
error_h_H1 = plot(hh_lg,error_h1_x,'--g*','LineWidth',2);
%ylim([-18 -2]);

legend('L_2 error convergence rate','H_2error convergence rate');
xlabel('log ||hh||/L ');
ylabel('log ||e||_{H1} ');



% figure
yyaxis right
error_h_H2 = plot(hh_lg,error_h2_x,'--b*','LineWidth',2);
%ylim([-18 -2]);

legend('L_2 error convergence rate','H_1error convergence rate',...
       'H_2error convergence rate');
xlabel('log ||hh||/L ');
ylabel('log ||e||_{H2} ');
exportgraphics(gca,['error_u_l2_h2' '.jpg']);

T = table(hh_x,error_l2_x,error_h1_x,error_h2_x,slop_l2,slop_h1,slop_h2,...
    'variableNames',{'hh_mesh','error_l2','error_H1','error_H2',...
    'l2 convergence rate','H1 convergence rate','H2 convergence rate'});
T
writetable(T);


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
% END