clear all;

num_sim = 100;


eig_result_hist = zeros(num_sim,5);
sim_results(num_sim) = struct();
for iii = 1:num_sim
    iii
    parameters;
    Stbl_test;
    eig_result_hist(iii,:) = [eig_max_diag, eig_max_ext, eig_max_S, eig_max_proposed, eig_max_cen]; 
    sim_results(iii).A = params.A;
    sim_results(iii).G = params.G;
    sim_results(iii).P_diag = P_opt_diag;
    sim_results(iii).K_diag = K_opt_diag;
    sim_results(iii).P_ext = P_opt_ext;
    sim_results(iii).K_ext = K_opt_ext;
    sim_results(iii).X_ext = X_opt_ext;
    sim_results(iii).P_S = P_opt_S;
    sim_results(iii).K_S = K_opt_S;
    sim_results(iii).P_proposed = P_opt_proposed;
    sim_results(iii).K_proposed = K_opt_proposed;
    sim_results(iii).X_proposed = X_opt_proposed;
    sim_results(iii).P_cen = P_opt_cen;
    sim_results(iii).K_cen = K_opt_cen;
end
%% 

num_of_success = [length(find(eig_result_hist(:,1)<1))
                  length(find(eig_result_hist(:,2)<1))
                  length(find(eig_result_hist(:,3)<1))
                  length(find(eig_result_hist(:,4)<1))
                  length(find(eig_result_hist(:,5)<1))];

fprintf('-------------------------------------------\n');
fprintf('------ Stabilizer Feasibility Result ------\n');
fprintf('  Num of Agents                  : %d \n', params.n);
fprintf('  Edge Generation Probability    : %.2f \n', params.p);
fprintf('  Num of Simulations             : %d \n', num_sim);
fprintf('  Num of Success   - K_diag      : %d \n', num_of_success(1));
fprintf('                   - K_ext       : %d \n', num_of_success(2)); 
fprintf('                   - K_S         : %d \n', num_of_success(3));
fprintf('                   - K_proposed  : %d \n', num_of_success(4));
fprintf('                   - K_cen       : %d \n', num_of_success(5));
fprintf('-------------------------------------------\n');

% save('data/stbl/stabilization_result_n10_p04_unstbl.mat', 'sim_results', 'num_of_success');
