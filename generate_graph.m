function G = generate_graph(n,prob)

    % arguments
    %     prob = 0.3 % probability to generate edges
    % end
    
    kkk = 20;
    ll = 0;
    
    Adj = zeros(n,n); % adjacency matrix
    pp = prob; % probability of generating an edge 
    
    for iii = 1:kkk
        for i = 1:n
            for  j = 1:n
                tmp = rand;
                if i > j && tmp >= 1-pp 
                    Adj(i,j) = 1;
                    Adj(j,i) = 1;
                    
                    ll = ll + 1;
                end
            end
        end
        eig_lap = eig(laplacian(graph(Adj)));
        if eig_lap(2) >0 %% connectivity
            break
        end
    end
    
    G = graph(Adj)
    if ll == kkk
        fprintf('error -- unconnected');
    end
end