function [gamma_opt,K_opt,P_opt] = Hinfty_diag(params,flag)

    n = params.n;
    G = params.G;

    A = params.A;
    B = params.B;
    C = params.C;
    D = params.D;

    Bw = params.Bw;
    Dw = params.Dw;

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
    gamma= sdpvar(1,1);
    eps = 10e-4;
    % gamma = 1.0981;
    
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

    LMI_1 = [-Q, A*Q+B*Z, Bw, zeros(n);
             (A*Q+B*Z)', -Q, zeros(n), (C*Q+D*Z)';
             Bw', zeros(n), -gamma * eye(n), Dw';
             zeros(n), C*Q + D*Z, Dw, -gamma * eye(n)];
    
    LMI = [LMI_1 <= 0, Q >= 0, gamma >=0];


    result = optimize(LMI,gamma,ops)
    %%  --------  end LMI  --------
    gamma_opt = value(gamma);
    if result.problem == 1
        gamma_opt = 0;
    end

    Z_opt = value(Z);
    Q_opt = value(Q);
    K_opt = Z_opt*inv(Q_opt);
    P_opt = inv(Q_opt);

    fprintf('-------------------------------------------\n');
    fprintf('-------- Block-diagonal relaxation --------\n')
    fprintf(' gamma_opt                      : %8.3e \n', gamma_opt);
    fprintf(' min of Ps eigval               : %8.2e \n', min(eig(P_opt)));
    fprintf(' condition number of P          : %8.2e \n', max(eig(P_opt))/min(eig(P_opt)));
    fprintf(' Norm of K                      : %8.2e \n', norm(K_opt)); 
    fprintf(' max of A+BKs eigval (abs)      : %8.2e \n', max( abs(eig( A + B*K_opt )) ));
    fprintf('-------------------------------------------\n');

end