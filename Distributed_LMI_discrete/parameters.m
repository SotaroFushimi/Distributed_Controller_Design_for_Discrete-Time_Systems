%% define parameters
%%%%%%%%% --------------- start ---------------
params.n = 8; %% the number of nodes 偶数とする
params.d = 1; %% dimension
params.p = 0.8;
params.G = generate_graph(params.n,params.p);
params.eps = 0;
% plot(params.G);
% system's parametersx

% define A's eigen values
eig_A_theta     = 2*pi* rand(params.n/2,1);
eig_A_theta     = [eig_A_theta;-eig_A_theta];
eig_A_r         = 1+rand(params.n/2,1);
eig_A_r         = [eig_A_r;eig_A_r];
eig_A           = complex(eig_A_r.*cos(eig_A_theta), eig_A_r.*sin(eig_A_theta));

% tmp = rand(params.n,params.n) .* laplacian(params.G);
tmp = rand(params.n,params.n);
params.A = tmp - place(tmp,eye(params.n), eig_A);
params.A = tmp \ params.A * tmp;

dd = 1;
%params.B = diag([ones(params.n-dd,1);zeros(dd,1)]);
params.B = [rand([params.n-dd,params.n]);zeros(dd,params.n)];
% params.B = eye(params.n);
% params.B = rand(params.n,params.n);

params.C = eye(params.n);
params.D = eye(params.n);
%% exogenous noises
% params.Bw = rand(params.n,params.n);
% params.Dw = rand(params.n,params.n);
params.Bw = eye(params.n);
params.Dw = eye(params.n);

%%%%%%%%% --------------  end  ---------------
fprintf('Rank deficiency of the contrability matrix - n:%8.2e \n', rank(ctrb(params.A,params.B))-params.n);


% params.solver_chosens = 'sedumi';
params.solver_chosens = 'sdpt3';