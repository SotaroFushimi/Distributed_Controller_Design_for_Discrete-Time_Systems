function [K_opt,P_opt,eig_max] = Stbl_S(params,flag)

    n = params.n;
    G = params.G;

    A = params.A;
    B = params.B;
    C = params.C;
    D = params.D;

    Bw = params.Bw;
    Dw = params.Dw;
    eps = params.eps;

    E = generate_Ematrix(n,G);
    EE = generate_Ematrix_cell(n,G);
    M = eye(size(E,1)) - E*inv(E'*E)*E';

    %% dilation
    A_e = dilation_by_E(A,E);
    B_e = dilation_by_E(B,E);
    % Dw_e = Dw;
    
    cliques = maximalCliques(adjacency(G)); 

    %%  -------- solve LMI -------- 
    ops = sdpsettings('solver',params.solver_chosens);
    ops.verbose = flag;

    Z_e = [];
    Q_e = []; 
    rho = sdpvar(1,1);
    eta = sdpvar(1,1);

    for l = 1:length(cliques)
        tmp_q = sdpvar(size(EE{l},1),size(EE{l},1),'symmetric');
        tmp_z = sdpvar(size(EE{l},1),size(EE{l},1),'full');
        Q_e = blkdiag(Q_e, tmp_q);
        Z_e = blkdiag(Z_e, tmp_z);
    end
    
    % Objective and constraints

    LMI = [];

    PHI = A_e*Q_e+B_e*Z_e;

    PHI = [Q_e, PHI';
           PHI, Q_e];
    Fins = rho*blkdiag(M,M);

    % LMI = [PHI+Fins>= 0, Q_e >= 0, M*Q_e+Q_e*M>=eta*M, eta >= 0];
    dim = size(Q_e, 1);
    LMI = [PHI+Fins>= eps*eye(2*dim), Q_e >= eps*eye(dim), M*Q_e+Q_e'*M-eta*M>=eps*eye(dim), eta >= eps];
    
    optimize(LMI,0,ops)
    %%  --------  end LMI  --------


    %% results 
    Z_opt = value(Z_e);
    Q_opt = value(Q_e);
    K_opt = inv(E'*E)*E'*(Z_opt*inv(Q_opt))*E;
    P_opt = E' * inv(Q_opt) * E;
    eig_max = max( abs(eig( A+B*K_opt)) );
    if min(eig(P_opt))<0 % (if the optimization fails)
        eig_max = 10;
    end

    fprintf('-------------------------------------------\n');
    fprintf('------------- Stabilizer K_S --------------\n')
    fprintf(' min of Ps eigval               : %8.2e \n', min(eig(P_opt)));
    fprintf(' condition number of P          : %8.2e \n', max(eig(P_opt))/min(eig(P_opt)));
    fprintf(' Norm of K                      : %8.2e \n', norm(K_opt));
    fprintf(' max of A+BKs eigval (abs)      : %8.2e \n', eig_max);
    fprintf('-------------------------------------------\n');

end