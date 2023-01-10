% simulation the beam with Hermite polynomial
clear all; clc; close all;

pi  = atan(1) * 4;
nqp = 10; %quadrature rule

EI = 1 ; % bend module stiffness
M  = 3 ; % Natual BC prescribed Moment
Q  = -12 ; % Natual BC prescribed shear force

omega_L = 0; % phsical length of beam
omega_R = pi; % position length of beam
LL = omega_R - omega_L; %  phsical length of beam

f_q = @(x) sin(x); %distribution body force

u_g = @(x) 0;   %essential BC transverse displacement
u_x_g = @(x) 0; %essential BC derivative of disp slope;

M_f = @(x) M / EI;% Natual BC prescribed Moment
Q_f = @(x) Q / EI;% Natual BC prescribed shear force

% exact solution

c1 = Q + 1.0;
c2 = M;
c3 = 1.0 - 0.5 * (Q + 1.0) * pi^2 - M * pi ;
c4 = (Q + 1.0) * pi^3 / 3 - pi + 0.5 * M * pi^2;
U_fx = @(x) (sin(x) + c1.*x.^3 / 6 + 0.5 * c2.*x.^2 + c3.*x + c4)./EI; % displacements

U_d1x = @(x) (cos(x) + 0.5 * c1 * x.^2 + c2 * x + c3)./EI;%first derivative
U_d2x = @(x) (-sin(x) + c1 .* x + c2) ./ EI;%second derivative
U_d3x = @(x) (-cos(x) + c1)./EI;%third derivative

%-------------------------------------------------------------------

nElem_x = 2 : 4 : 32;
%nElem_x = 2;

cn         = length(nElem_x);

hh_x       = zeros(cn,1);
error_l2_x = zeros(cn,1);
error_h2_x = zeros(cn,1);

hh_lg   = zeros(cn,1);
slop_l2 = zeros(cn,1);
slop_h2 = zeros(cn,1);

for num = 1 : cn

    nElem   = nElem_x(num); %number of elements
    n_en    = 2;      % node number of element
    n_np    = nElem + 1;
    d_node  = 2;      % node degree of freedom: displacements and slopes
    d_en    = n_en*2; % element DOF
    nLocBas = d_en;   % local Basis equal to DOF of element
    nFunc   = nLocBas + d_node * (nElem-1);

    %assemble IEN
    IEN  = zeros(d_en, nElem);
    for ee = 1 : nElem
        for aa = 1 : d_en
            IEN (aa, ee) = (ee - 1) * d_node + aa;
        end
    end

    %mesh
    hh_x(num) = (omega_R - omega_L) / nElem;

    x_coor = omega_L : hh_x(num) : omega_R;

    %setup ID array based on BC

    ID = 1 : nFunc;
    % assign the ID for the Dirichlet node to be -1.
    % the right of beam is fixed, transverse dis and slope equal = 0.
    ID(end-1:end) = -1;
   
    n_eq = nFunc - 2;   % Dirichlet nodes is known
 
    %allocate an empty stiffness matrix and load vector
    %     K = sparse(nFunc, nFunc);
    %     F = zeros( nFunc, 1);
    K = sparse(n_eq, n_eq);
    F = zeros(n_eq, 1);

    %-------------------------------------------------------------------
    %Assembly the stiffness matrix and load vector
    for ee = 1 : nElem
        %allocate zero element stiffness matrix and element load vector
        x_ele = zeros(n_en,1);
        u_ebc  = zeros(nLocBas,1);

        % New data structure: mapping local element information to global.

        for aa = 1 : n_en
            x_ele(aa) = x_coor(ee+aa-1) ; % physical geometrical position;
        end

        [qp, wq] = Gauss(nqp, x_ele(1), x_ele(end));

        for qua = 1 : nqp

            for aa = 1 : nLocBas
                Na    = Hermiteg_Basis(aa, 0, qp(qua), x_ele(1), x_ele(end));
                Na_x  = Hermiteg_Basis(aa, 1, qp(qua), x_ele(1), x_ele(end));
                Na_2x = Hermiteg_Basis(aa, 2, qp(qua), x_ele(1), x_ele(end));
                %calculate the global K matrix and force vector

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
                            % TO BE FILLED
                            Nb_2x = Hermiteg_Basis(bb, 2, qp(qua), x_ele(1), x_ele(end));
                            %obtain the Dirichlet node's physical coordinates
                            x_g_u = x_coor(ee+1) ;
                            % Obtain the boundary data at this point

                            % esstential BC: Node displacement; Node slope
                            u_ebc = [0; 0; u_g(x_g_u); u_x_g(x_g_u) ];

                            F( AA ) = F( AA ) -  wq(qua) * EI * Na_2x * Nb_2x * u_ebc(bb); 

                        end
                    end
                end
            end

        end %End of quadrature loop

        % modified the load vector by natural BC
        if ee == 1
            F(ID(IEN(1,ee))) =  F(ID(IEN(1,ee))) + Q_f(x_coor(ee));  % prescribed shear force
            F(ID(IEN(2,ee))) =  F(ID(IEN(2,ee))) - M_f(x_coor(ee));  % prescribed Moment force
        end

    end

    %disp(K);disp(F);
    %solve the stiffness matrix problem
  
    Uh = K \ F;
    % Append the displacement vector by the Dirichlet data
    Uh = [ Uh; u_g(x_g_u); u_x_g(x_g_u) ];
    %-------------------------------------------------------------------
    % displacement distribution
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
    exportgraphics(gca,['file_fai' '.jpg']);

    %-------------------------------------------------------------------
    % error convergence analysis
    error_l2  = 0.0;
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
            u_dx2  = 0.0 ;

            for aa = 1 : nLocBas

                U_h    = U_h   + u_ele(aa) * Hermiteg_Basis(aa,0,qp(qua),x_ele(1),x_ele(end));
                u_dx2  = u_dx2 + u_ele(aa) * Hermiteg_Basis(aa,2,qp(qua),x_ele(1),x_ele(end));

            end

            u_exact     = U_fx(qp(qua));
            u_dx2_exact = U_d2x(qp(qua));

            error_l2 = error_l2 + wq(qua) * (U_h - u_exact) * (U_h - u_exact);
            error_H2 = error_H2 + wq(qua) * (u_dx2 - u_dx2_exact) * (u_dx2 - u_dx2_exact);

        end
    end

    error_l2 = sqrt(error_l2);
    error_H2 = sqrt(error_H2);

    error_l2_x(num,1)  = log(error_l2);
    error_h2_x(num,1)  = log(error_H2);
    hh_lg(num) = log(hh_x(num) / LL); % normlized the mesh length size


    dx_hh  = zeros( 2,1);
    dy_l2e = zeros( 2,1);
    dy_h2e = zeros( 2,1);
    if (num >= 2)
        dx_hh = [hh_lg(num);hh_lg(num-1)];

        dy_l2e = [error_l2_x(num);error_l2_x(num-1)];
        dy_h2e = [error_h2_x(num);error_h2_x(num-1)];

        slop_l2(num) = (dy_l2e(2) - dy_l2e(1)) / (dx_hh(2)-dx_hh(1));
        slop_h2(num) = (dy_h2e(2) - dy_h2e(1)) / (dx_hh(2)-dx_hh(1));

    else

        slop_l2(num) = 0;
        slop_h2(num) = 0;
    end

end


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

T = table(hh_x,error_l2_x,error_h2_x,slop_l2,slop_h2,...
    'variableNames',{'hh_mesh','error_l2','error_H2',...
    'l2 convergence rate','H2 convergence rate'});
writetable(T);
T

% END