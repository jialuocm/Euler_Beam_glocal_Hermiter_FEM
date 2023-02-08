% =========================================================================
% This is the FEM Matlab main code for one-dimensional Bernulli-Euler beam.
%
% In this main, adopte the Hermite type basis function.
%
% -------------------------------------------------------------------------
% By Jia Luo, 2023 Jan. 11th.
% =========================================================================
clear all; clc; close all;

pi  = atan(1) * 4;
nqp = 10;                        % Quadrature rule

EI = 0.5 ;                         % Bend module stiffness
M  = 0 ;                         % Natual BC prescribed Moment
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
U_fx = @(x) (sin(x) + c1.*x.^3 / 6 + 0.5 * c2.*x.^2 + c3.*x + c4)./EI; % displacements

U_d1x = @(x) (cos(x) + 0.5 * c1 * x.^2 + c2 * x + c3)./EI;            %first derivative
U_d2x = @(x) (-sin(x) + c1 .* x + c2) ./ EI;                          %second derivative
%-------------------------------------------------------------------


nElem_x = 3 : 1 : 9;
%nElem_x = 2;

cn         = length(nElem_x);
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
    n_en    = 2;              % Node number of local element
    n_np    = nElem + 1;      % Number nodes of elements
    d_node  = 2;              % Node degree of freedom: displacements and slopes
    d_en    = n_en*2;         % Element DOF
    nLocBas = d_en;           % Local Basis equal to DOF of element
    nFunc   = nLocBas + d_node * (nElem-1); % 

    % Assemble IEN
    IEN  = zeros(d_en, nElem);
    for ee = 1 : nElem
        for aa = 1 : d_en
            IEN (aa, ee) = (ee - 1) * d_node + aa;
        end
    end

    % Mesh of geometric domain
    hh_x(num) = (omega_R - omega_L) / nElem;
    x_coor    = omega_L : hh_x(num) : omega_R;

    % Setup ID array based on BC
    ID = 1 : nFunc;
    % Assign the ID for the Dirichlet node to be -1
    ID(end-1:end) = -1;
    n_eq = nFunc - 2;  

    % Allocate an empty stiffness matrix and load vector
    K = sparse(n_eq, n_eq);
    F = zeros(n_eq, 1);

    %-------------------------------------------------------------------
    % Assembly the stiffness matrix and load vector
    for ee = 1 : nElem
        % Allocate zero element stiffness matrix and element load vector
        x_ele = zeros(n_en,1);
        u_ebc  = zeros(nLocBas,1);

        % New data structure: mapping local element information to global

        for aa = 1 : n_en
            x_ele(aa) = x_coor(ee+aa-1) ; % Physical geometrical position
        end

        [qp, wq] = Gauss(nqp, x_ele(1), x_ele(end));

        for qua = 1 : nqp

            for aa = 1 : nLocBas
                Na    = Hermiteg_Basis(aa, 0, qp(qua), x_ele(1), x_ele(end));
                Na_x  = Hermiteg_Basis(aa, 1, qp(qua), x_ele(1), x_ele(end));
                Na_2x = Hermiteg_Basis(aa, 2, qp(qua), x_ele(1), x_ele(end));
                % Calculate the global K matrix and force vector

                % LM,Location matrix, given a particular degree of freedom number and an element number,
                % return the corresponding global number equation number.
                AA = ID( IEN(aa,ee) );

                if( AA > 0 )
                    F(AA) = F(AA) + wq (qua) * Na * f_q(qp(qua));
                    for bb = 1: nLocBas
                        BB = ID( IEN(bb,ee) );
                        if (BB > 0)
                            Nb_2x = Hermiteg_Basis(bb, 2, qp(qua), x_ele(1), x_ele(end));
                            K(AA,BB) = K(AA,BB) + wq(qua) * EI * Na_2x * Nb_2x;
                        else
                 
                            Nb_2x = Hermiteg_Basis(bb, 2, qp(qua), x_ele(1), x_ele(end));
                            % Obtain the Dirichlet node's physical coordinates
                            x_g_u = x_coor(ee+1) ;

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
            F(ID(IEN(1,ee))) =  F(ID(IEN(1,ee))) + Q_f(x_coor(ee));  % prescribed shear force
            F(ID(IEN(2,ee))) =  F(ID(IEN(2,ee))) - M_f(x_coor(ee));  % prescribed Moment force
        end

    end

    % Solve the stiffness matrix problem
  
    Uh = K \ F;
    % Append the displacement vector by the Dirichlet data
    Uh = [ Uh; u_g(x_g_u); u_x_g(x_g_u) ];
    %-------------------------------------------------------------------
    % Displacement distribution
    uh_od    = Uh(1:2:nFunc - 1);
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

    %-------------------------------------------------------------------
    % slope distribution
    figure
    faih_p = plot(x_coor,fai_even,'--r*','linewidth',2);

    hold on
    p_M_exact = plot(x_coor,U_d1x(x_coor),'b--','LineWidth',2);
    legend([faih_p,p_M_exact],'\phi_h FEM','\phi_{Exact}');
    xlabel('x coor');
    ylabel('\phi ');
    hold off
%     exportgraphics(gca,['file_fai' '.jpg']);

    %-------------------------------------------------------------------
    % Error convergence analysis
    error_l2  = 0.0;
    error_H1  = 0.0;
    error_H2  = 0.0;

    jj=1;
    M_h     = zeros(nElem*nqp,1);
    M_exact = zeros(nElem*nqp,1);
    x_m     = zeros(nElem*nqp,1);
    x_b     = zeros(nElem+1,1);
    U_h_p   = zeros(nElem*nqp,1);

    for ee = 1: nElem
        u_ele = zeros(nLocBas,1);
        x_ele = zeros(n_en,1);

        for aa = 1 : nLocBas
            u_ele(aa) = Uh(IEN(aa,ee));
        end

        for aa = 1 : n_en
            x_ele(aa) = x_coor(ee+aa-1) ;
        end

        [qp, wq] = Gauss(nqp, x_ele(1), x_ele(end));

        for qua = 1 : nqp

            U_h    = 0.0;
            u_hdx1  = 0.0 ;
            u_hdx2  = 0.0 ;

            for aa = 1 : nLocBas

                U_h    = U_h   + u_ele(aa) * Hermiteg_Basis(aa,0,qp(qua),x_ele(1),x_ele(end));
                u_hdx1  = u_hdx1 + u_ele(aa) * Hermiteg_Basis(aa,1,qp(qua),x_ele(1),x_ele(end));
                u_hdx2  = u_hdx2 + u_ele(aa) * Hermiteg_Basis(aa,2,qp(qua),x_ele(1),x_ele(end));

            end

            u_exact     = U_fx(qp(qua));
            u_dx1_exact = U_d1x(qp(qua));
            u_dx2_exact = U_d2x(qp(qua));

            error_l2 = error_l2 + wq(qua) * (U_h - u_exact) * (U_h - u_exact);
            error_H1 = error_H1 + wq(qua) * (u_hdx1 - u_dx1_exact) * (u_hdx1 - u_dx1_exact);

            error_H2 = error_H2 + wq(qua) * (u_hdx2 - u_dx2_exact) * (u_hdx2 - u_dx2_exact);

        end
    end

    error_l2 = sqrt(error_l2);
    error_H1 = sqrt(error_H1);
    error_H2 = sqrt(error_H2);

    error_l2_x(num,1)  = log(error_l2);
    error_h1_x(num,1)  = log(error_H1);
    error_h2_x(num,1)  = log(error_H2);
    hh_lg(num) = log(hh_x(num) / LL); % Normlized the mesh length size

    dx_hh   = zeros( 2,1);
    e_dy_l2 = zeros( 2,1);
    e_dy_H1 = zeros( 2,1);
    e_dy_H2 = zeros( 2,1);
    if (num >= 2)
        dx_hh = [hh_lg(num);hh_lg(num-1)];

        e_dy_l2 = [error_l2_x(num);error_l2_x(num-1)];
        e_dy_H1= [error_h1_x(num);error_h1_x(num-1)];
        e_dy_H2 = [error_h2_x(num);error_h2_x(num-1)];

        slop_l2(num) = (e_dy_l2(2) - e_dy_l2(1)) / (dx_hh(2)-dx_hh(1));
        slop_h1(num) = (e_dy_H1(2) - e_dy_H1(1)) / (dx_hh(2)-dx_hh(1));
        slop_h2(num) = (e_dy_H2(2) - e_dy_H2(1)) / (dx_hh(2)-dx_hh(1));

    else

        slop_l2(num) = 0;
        slop_h1(num) = 0;
        slop_h2(num) = 0;
    end

end
%-------------------------------------------------------------------
% Postprocessing
figure
yyaxis left
error_h_l2 = plot(hh_lg,error_l2_x,'--rO','LineWidth',2);
ylim([-18 -2]);

legend(error_h_l2,'L_2 error convergence analysis');
xlabel('log ||hh||/L ');
ylabel('log ||e||_{l2} ');
% exportgraphics(gca,['error_u_l2' '.jpg']);

hold on;
% figure
yyaxis right
error_h_H2 = plot(hh_lg,error_h2_x,'--b*','LineWidth',2);
ylim([-18 -2]);

legend('L_2 error convergence rate','H_2error convergence rate');
xlabel('log ||hh||/L ');
ylabel('log ||e||_{H2} ');
exportgraphics(gca,['error_u_l2_h2' '.jpg']);

T = table(hh_x,error_l2_x,error_h1_x,error_h2_x,slop_l2,slop_h1,slop_h2,...
    'variableNames',{'hh_mesh','error_l2','error_H1','error_H2',...
    'l2 convergence rate','H1 convergence rate','H2 convergence rate'});
writetable(T);
disp(T);

% END