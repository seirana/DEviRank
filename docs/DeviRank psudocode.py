Input:
    G = (V, E)           // PPI network with edge weights PPI_e
    T ⊆ V               // Disease target proteins
    S ⊆ V               // Drug target proteins
    DGI(s)              // Drug–gene interaction score for each s ∈ S

Output:
    w_d(d)              // Final drug score

// Step 1: Disease protein weights
for each t ∈ T do
    w_t(t) ← sum of PPI weights of edges incident to t
end for

// Step 2: Drug target protein weights
for each s ∈ S do
    w_s(s) ← 0
    for each t ∈ T do
        for each path p ∈ P[s → t] do
            path_score ← product of PPI weights along p
            w_s(s) ← w_s(s) + path_score × (1 + w_t(t))
        end for
    end for
end for

// Step 3: Final drug score
w_d(d) ← 0
for each s ∈ S do
    w_d(d) ← w_d(d) + DGI(s) × w_s(s)
end for

return w_d(d)
