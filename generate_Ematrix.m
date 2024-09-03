function E = generate_Ematrix(n,G)
    E = []; %% CD matrix for maximal cliques

    Adj = adjacency(G);
    cliques = maximalCliques(Adj);
    for l = 1:length(cliques)
        tmp = zeros(length(cliques{l}),n);
    
        for j = 1:length(cliques{l})
            tmp(j, cliques{l}(j) ) = 1; 
        end
        E = [E; tmp];
        % EE{l} = tmp;
    end
end