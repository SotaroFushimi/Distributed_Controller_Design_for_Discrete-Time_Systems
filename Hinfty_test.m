gamma_result = [];
stab_result = [];

%% Hinfty
 
% block-diagonal relaxation:
[gamma_opt_diag,K_opt_diag,P_opt_diag] = Hinfty_diag(params,0);

% extended LMI
[gamma_opt_ext,K_opt_ext,P_opt_ext,X_opt_ext] = Hinfty_ext(params,0);

% K_S:
[gamma_opt_S,K_opt_S,P_opt_S] = Hinfty_S(params,0);

% K_proposed:
[gamma_opt_proposed,K_opt_proposed,P_opt_proposed,X_opt_proposed] = Hinfty_proposed(params,0);

% centralized controller:
[gamma_opt_cen,K_opt_cen,P_opt_cen] = Hinfty_centralized(params,0);

gamma_result = [gamma_opt_diag, gamma_opt_ext, gamma_opt_S, gamma_opt_proposed, gamma_opt_cen];
stab_result = [check_stab(params,K_opt_diag)
               check_stab(params,K_opt_ext)
               check_stab(params,K_opt_S)
               check_stab(params,K_opt_proposed)
               check_stab(params,K_opt_cen)];
% 
function output = check_stab(params,K)
    output =  max(abs(eig(params.A+params.B*K)));
end