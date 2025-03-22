function dual_admissible_pivots(BT::BasisTableau)
    @assert !dual_unbounded_certificate(BT) # too strong a check for learning purposes
    @assert !is_optimal(BT)
    (m, n) = size(BT.A);  btab = BT.btab;  x_B = btab[2:end, 1]
    P = [𝓑(BT,i) for i in 1:m if x_B[i] < 0] # {p ∈ I_B : x_p < 0}

    #s = -(btab[1, 2:end]);  T = btab[2:end,2:end];  I_N = BT.indices_nonbasic
    s_N = -btab[1, (1+m) .+ (1:n-m)];  T = btab[2:end,2:end];  I_N = BT.indices_nonbasic
    println(maximum(((p,q),) -> s_N[𝓝⁻¹(BT, q)] / -T[𝓣(BT,p), 𝓣(BT,q)], filter( ((p,q),) -> T[𝓣(BT,p), 𝓣(BT,q)] < 0, collect(Iterators.product(P, I_N)))))
    argmin_set(
        #((p,q),) -> s[𝓣(BT,q)] / -T[𝓣(BT,p), 𝓣(BT,q)],
        ### conceptually, T[𝓣(BT,p), 𝓣(BT,q)] = (A_B^-1 * A_N)[𝓑⁻¹(p), 𝓝⁻¹(q)] ###
        ((p,q),) -> s_N[𝓝⁻¹(BT, q)] / -T[𝓣(BT,p), 𝓣(BT,q)],
        filter( ((p,q),) -> T[𝓣(BT,p), 𝓣(BT,q)] < 0, collect(Iterators.product(P, I_N))))
#=
    # for ease of understanding
    θ = minimum([s[𝓣(BT,q)] / -T[𝓣(BT,p), 𝓣(BT,q)]
                    for s = -(btab[1, 2:end])
                    for p ∈ P, q ∈ I_N
                    if T[𝓣(BT,p), 𝓣(BT,q)] < 0])
=#
end
function pivotChoice_least_index_dual_simplex(BT::BasisTableau)    
    (dual_admissible_pivots(BT)
        |> filter( ((i,j),) -> argmin_set(((i,j),) -> i, pairs)) # minimal i
        |> filter( ((i,j),) -> argmin_set(((i,j),) -> j, pairs)) # minimal j
        |> first) # (first, and really only, pair left)
end
# TO-DO: double-check condition on while-loop
function dual_simplex_method(c, A, b, initial_basis; pivotChoice = pivotChoice_least_index_dual_simplex)
    tableaux = Vector{BasisTableau}(undef, 0)
    push!(tableaux, BasisTableau(c, A, b, initial_basis))
    while !(is_optimal(tableaux[end]) || primal_unbounded_certificate(tableaux[end])
                                      || dual_unbounded_certificate(tableaux[end]))
        push!(tableaux, pivot(tableaux[end], pivotChoice(tableaux[end])...));
    end
    return tableaux
end