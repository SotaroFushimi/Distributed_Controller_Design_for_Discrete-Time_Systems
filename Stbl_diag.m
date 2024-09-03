function [K_opt,P_opt,eig_max] = Stbl_diag(params,flag)

    n = params.n;
    G = params.G;

    A = params.A;
    B = params.B;
    C = params.C;
    D = params.D;

    Bw = params.Bw;
    Dw = params.Dw;
    eps = params.eps;

    % E = generate_Ematrix(n,G);
    % EE = generate_Ematrix_cell(n,G);
    % M = eye(size(E,1)) - E*inv(E'*E)*E';

    % %% dilation
    % A_e = dilation_by_E(A,E);
    % B_e = dilation_by_E(B,E);
    % C_e = C * inv(E'*E) * E';
    % D_e = D * inv(E'*E) * E';
    % Bw_e = E * Bw;
    % % Dw_e = Dw;
    
    % cliques = maximalCliques(adjacency(G)); 

    %%  -------- solve LMI -------- 
    ops = sdpsettings('solver',params.solver_chosens);
    ops.verbose = flag;

    Z = [];
    Q = diag(sdpvar(n,1));
    
    L = laplacian(G);
    for i = 1:n
        tmp = [];
        for j = 1:n
            if L(i,j) ~= 0
                tmp = [tmp sdpvar(1,1)];
            else
                tmp = [tmp 0];
            end
        end
        Z = [Z;tmp];
    end
    
    % Objective and constraints
    
    LMI = [];

    LMI_1 = [Q,(A*Q+B*Z)';
             A*Q+B*Z,Q];
    
    % LMI = [LMI_1 >= 0, Q >= 0];
    LMI = [LMI_1 >= eps*eye(2*n), Q >= eps*eye(n)]; 


    optimize(LMI,[],ops)
    %%  --------  end LMI  --------


    %% results 
    Z_opt = value(Z);
    Q_opt = value(Q);
    K_opt = Z_opt*inv(Q_opt);
    P_opt = inv(Q_opt);
    eig_max=max( abs(eig( A + B*K_opt )));

    if min(eig(P_opt))<0 % (if the optimization fails)
        eig_max = 10;
    end

    fprintf('-------------------------------------------\n');
    fprintf('----Stabilizer Block-diagonal relaxation---\n')
    fprintf(' min of Ps eigval               : %8.2e \n', min(eig(P_opt)));
    fprintf(' condition number of P          : %8.2e \n', max(eig(P_opt))/min(eig(P_opt)));
    fprintf(' Norm of K                      : %8.2e \n', norm(K_opt)); 
    fprintf(' max of A+BKs eigval (abs)      : %8.2e \n', eig_max);
    fprintf('-------------------------------------------\n');

end