function [gamma_opt,K_opt,P_opt,X_opt] = Hinfty_ext(params,flag)

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
    Q = sdpvar(n);
    gamma= sdpvar(1,1);
    X = diag(sdpvar(n,1));
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
    
    LMI_1 = [-Q, A*X+B*Z, Bw, zeros(n);
             (A*X+B*Z)', Q-X-X', zeros(n), (C*X+D*Z)';
             Bw', zeros(n), -gamma*eye(n), Dw';
             zeros(n), C*X + D*Z, Dw, -gamma * eye(n)];
    
    LMI = [LMI_1<= 0, Q >= 0, gamma>=0];

    result = optimize(LMI,gamma,ops)
    %%  --------  end LMI  --------
    gamma_opt = value(gamma);
    if result.problem == 1
        gamma_opt = 0;
    end
    Z_opt = value(Z);
    X_opt = value(X);
    K_opt = Z_opt*inv(X_opt);
    P_opt = inv(value(Q));
    % eig(P_opt)
    % eig((A+B*K_opt)'*P_opt*(A+B*K_opt)-P_opt)



    fprintf('-------------------------------------------\n');
    fprintf('-------------- Extended LMI ---------------\n')
    fprintf(' gamma_opt                      : %8.3e \n', gamma_opt);
    fprintf(' min of Ps eigval               : %8.2e \n', min(eig(P_opt)));
    fprintf(' condition number of P          : %8.2e \n', max(eig(P_opt))/min(eig(P_opt)));
    fprintf(' Norm of K                      : %8.2e \n', norm(K_opt)); 
    fprintf(' max of A+BKs eigval (abs)      : %8.2e \n', max( abs(eig( A + B*K_opt )) ));
    fprintf('-------------------------------------------\n');

end