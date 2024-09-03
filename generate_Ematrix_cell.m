function EE = generate_Ematrix(n,G)
    EE = {}; %% {E_1, E_2, ....}    
    
    Adj = adjacency(G);
    cliques = maximalCliques(Adj);
    for l = 1:length(cliques)
        tmp = zeros(length(cliques{l}),n);
    
        for j = 1:length(cliques{l})
            tmp(j, cliques{l}(j) ) = 1; 
        end
        
        EE{l} = tmp;
    end
end