%% Stabilization
 
% block-diagonal relaxation:
[K_opt_diag, P_opt_diag, eig_max_diag] = Stbl_diag(params,0);

% extended LMI:
[K_opt_ext, P_opt_ext, eig_max_ext, X_opt_ext] = Stbl_ext(params,0);

% K_S:
[K_opt_S, P_opt_S, eig_max_S] = Stbl_S(params,0);

% proposed method:
[K_opt_proposed, P_opt_proposed, eig_max_proposed, X_opt_proposed] = Stbl_proposed(params,0);

% centralized controller:
[K_opt_cen, P_opt_cen, eig_max_cen] = Stbl_centralized(params,0);