function [K_opt,P_opt,eig_max,X_opt] = Stbl_proposed(params,flag)

    n = params.n;
    G = params.G;

    A = params.A;
    B = params.B;
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
    X_e = []; 
    rho = sdpvar(1,1);
    eta = sdpvar(1,1);

    for l = 1:length(cliques)
        tmp_x = sdpvar(size(EE{l},1),size(EE{l},1),'full');
        tmp_z = sdpvar(size(EE{l},1),size(EE{l},1),'full');
        X_e = blkdiag(X_e, tmp_x);
        Z_e = blkdiag(Z_e, tmp_z);
    end

    dim = size(Z_e,1);
    Q_e = sdpvar(dim);
    
    % Objective and constraints

    LMI = [];

    PHI = A_e*X_e+B_e*Z_e;

    PHI = [X_e+X_e'-Q_e, PHI';
           PHI, Q_e];
    Fins = rho*blkdiag(M,M);

    % LMI = [PHI+Fins>= 0, Q_e >= 0, M*X_e+X_e'*M>=eta*M, eta >= 0];
    LMI = [PHI+Fins>= eps*eye(2*dim), Q_e >= eps*eye(dim), M*X_e+X_e'*M-eta*M>=eps*eye(dim), eta >= eps];
    
    optimize(LMI,0,ops)
    %%  --------  end LMI  --------


    %% results 
    Z_opt = value(Z_e);
    Q_opt = value(Q_e);
    X_opt = value(X_e);
    K_opt = inv(E'*E)*E'*(Z_opt*inv(X_opt))*E;
    P_opt = (E'*inv(X_opt)*E)* inv(E'*inv(X_opt')*Q_opt*inv(X_opt)*E)* (E'*inv(X_opt)*E)';
    eig_max = max( abs(eig( A+B*K_opt)) );

    if min(eig(P_opt))<0 % (if the optimization fails)
        eig_max = 10;
    end

    fprintf('-------------------------------------------\n');
    fprintf('------- Stabilizer Proposed method --------\n')
    fprintf(' min of Ps eigval               : %8.2e \n', min(eig(P_opt)));
    fprintf(' min of Xs eigval               : %8.2e \n', min(eig(X_opt)));
    fprintf(' condition number of P          : %8.2e \n', max(eig(P_opt))/min(eig(P_opt)));
    fprintf(' Norm of K                      : %8.2e \n', norm(K_opt));
    fprintf(' max of A+BKs eigval (abs): %8.2e \n', eig_max);
    fprintf('-------------------------------------------\n');

end