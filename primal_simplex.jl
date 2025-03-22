function primal_admissible_pivots(BT::BasisTableau)
#    @assert !primal_unbounded_certificate(BT) # too strong a check for learning purposes
    @assert !is_optimal(BT)
    (m, n) = size(BT.A);  btab = BT.btab; s_N = -btab[1, (1+m) .+ (1:n-m)]
    #display("s_N is $s_N")
    Q = [𝓝(BT,j) for j in 1:(n-m) if s_N[j] < 0] # {q ∈ I_N : s_q < 0}
    #display(Q)
    x_B = btab[2:end, 1];  T = btab[2:end,2:end];  I_B = BT.indices_basic

    set_of_pairs = filter( ((p,q),) -> T[𝓣(BT,p), 𝓣(BT,q)] > 0, collect(Iterators.product(I_B, Q))) # debug
    f = ((p,q),) -> x_B[𝓑⁻¹(BT,p)] / T[𝓣(BT,p), 𝓣(BT,q)]
    #display(collect(Iterators.product(I_B, Q)))
    #display(set_of_pairs)
    display(maximum(f, set_of_pairs))
    return argmin_set(
        ((p,q),) -> x_B[𝓑⁻¹(BT,p)] / T[𝓣(BT,p), 𝓣(BT,q)],
        filter( ((p,q),) -> T[𝓣(BT,p), 𝓣(BT,q)] > 0, collect(Iterators.product(I_B, Q))))
end
function pivotChoice_least_index_primal_simplex(BT::BasisTableau)    
    return (dual_admissible_pivots(BT)
        |> filter( ((i,j),) -> argmin_set(((i,j),) -> i, pairs)) # minimal i
        |> filter( ((i,j),) -> argmin_set(((i,j),) -> j, pairs)) # minimal j
        |> first) # (first, and really only, pair left)
end

# TO-DO: double-check condition on while-loop
# TO-DO: maybe factor out method? one method for primal simplex, dual simplex, criss-cross
function primal_simplex_method(c, A, b, initial_basis; pivotChoice = pivotChoice_least_index_dual_simplex)
    tableaux = Vector{BasisTableau}(undef, 0)
    push!(tableaux, BasisTableau(c, A, b, initial_basis))
    while !(is_optimal(tableaux[end]) || primal_unbounded_certificate(tableaux[end])
                                      || dual_unbounded_certificate(tableaux[end]))
        push!(tableaux, pivot(tableaux[end], pivotChoice(tableaux[end])...));
    end
    return tableaux
end

# TO-DO: Phase-I