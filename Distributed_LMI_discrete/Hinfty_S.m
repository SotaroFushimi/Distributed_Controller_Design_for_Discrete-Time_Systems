function [gamma_opt,K_opt,P_opt] = Hinfty_S(params,flag)

    n = params.n;
    G = params.G;

    A = params.A;
    B = params.B;
    C = params.C;
    D = params.D;

    Bw = params.Bw;
    Dw = params.Dw;

    E = generate_Ematrix(n,G);
    EE = generate_Ematrix_cell(n,G);
    M = eye(size(E,1)) - E*inv(E'*E)*E';

    %% dilation
    A_e = dilation_by_E(A,E);
    B_e = dilation_by_E(B,E);
    C_e = C * inv(E'*E) * E';
    D_e = D * inv(E'*E) * E';
    Bw_e = E * Bw;
    % Dw_e = Dw;
    
    cliques = maximalCliques(adjacency(G)); 

    %%  -------- solve LMI -------- 
    yalmip('clear')
    ops = sdpsettings('solver',params.solver_chosens);
    ops.verbose = flag;

    Z_e = [];
    Q_e = []; 
    gamma= sdpvar(1,1);
    eta = sdpvar(1,1);
    rho1 = sdpvar(1,1);
    eps = 1e-5;

    for l = 1:length(cliques)
        tmp_q = sdpvar(size(EE{l},1),size(EE{l},1),'symmetric');
        tmp_z = sdpvar(size(EE{l},1),size(EE{l},1),'full');
        Q_e = blkdiag(Q_e, tmp_q);
        Z_e = blkdiag(Z_e, tmp_z);
    end
    
    dim = size(Q_e,1);
    n   = size(C, 1);
    % Objective and constraints

    LMI = [];

    PHI = A_e*Q_e+B_e*Z_e;

    LMI_1 = [-Q_e, PHI, Bw_e, zeros(dim,n);
             PHI', -Q_e, zeros(dim,n), (C_e*Q_e+D_e*Z_e)';
            Bw_e', zeros(n,dim), -gamma * eye(n), Dw';
            zeros(n,dim), C_e*Q_e + D_e*Z_e, Dw, -gamma * eye(n)];

    Theta = blkdiag(M,M,zeros(n),zeros(n));

    % LMI = [LMI_1-rho1*Theta <= -eps*eye(2*dim+2*n), Q_e >= eps*eye(dim), gamma>=0, M*Q_e+Q_e*M>=eta*M, eta >= eps];
    LMI = [LMI_1-rho1*Theta <= 0, Q_e >= 0, gamma>=0, M*Q_e+Q_e*M>=eta*M, eta >= 0];
    result = optimize(LMI,gamma,ops)
    %%  --------  end LMI  --------
    gamma_opt = value(gamma);
    if result.problem == 1
        gamma_opt = 0;
    end
    Z_opt = value(Z_e);
    Q_opt = value(Q_e);
    K_opt = inv(E'*E)*E'*(Z_opt*inv(Q_opt))*E;
    P_opt = E' * inv(Q_opt) * E;
    % eig(P_opt)
    % eig((A+B*K_opt)'*P_opt*(A+B*K_opt)-P_opt)
    

    fprintf('-------------------------------------------\n');
    fprintf('------ Clique-wise decomposition K_S ------\n')
    fprintf(' gamma_opt                      : %8.3e \n', gamma_opt);
    fprintf(' min of Ps eigval               : %8.2e \n', min(eig(P_opt)));
    fprintf(' condition number of P          : %8.2e \n', max(eig(P_opt))/min(eig(P_opt)));
    fprintf(' Norm of K                      : %8.2e \n', norm(K_opt));
    fprintf(' max of A+BKs eigval (abs)      : %8.2e \n', max( abs(eig( A + B*K_opt )) ));
    fprintf('-------------------------------------------\n');

end