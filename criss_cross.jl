Type_I_admissible_pivots(BT::BasisTableau) =
    [(𝓑(BT, i), 𝓝(BT, j))
        for (m, n) = [size(BT.A)], btab = [BT.btab]
        for i in 1:m, j in 1:(n-m)
        if btab[1+i, 1] < 0 && btab[1+i, 1+m+j] < 0]

Type_II_admissible_pivots(BT::BasisTableau) =
    [(𝓑(BT, i), 𝓝(BT, j))
        for (m, n) = [size(BT.A)], btab = [BT.btab]
        for i in 1:m, j in 1:(n-m)
        if btab[1, 1+m+j] > 0 && btab[1+i, 1+m+j] > 0]

# returns leaving variable p and entering variable q as specified by the criss-cross method 
function pivotChoice_least_index_criss_cross(BT::BasisTableau)
    m, n = size(BT.A);   I_B, I_N = BT.indices_basic, BT.indices_nonbasic
    p = filter(i -> x(BT)[i] < 0, I_B) |> set -> ~isempty(set) ? minimum(set) : typemax(Int)
    q = filter(j -> s(BT)[j] < 0, I_N) |> set -> ~isempty(set) ? minimum(set) : typemax(Int)
    if p ≤ q
        @assert !primal_inconsistent_certificate(BT::BasisTableau)
#=        
        q = Type_I_admissible_pivots(BT)
        q = q |> filter( ((i,j),) -> i==p)
        q = q |> pairs -> argmin_set(((i,j),) -> j, pairs) # Type-II admissible (p,j) with minimal j
        q = q |> ((i,j),) -> j # extract entering variable from the single remaining pair
=#            
        q = (Type_I_admissible_pivots(BT)
            |> filter( ((i,j),) -> i==p) # filter out pairs that don't have p as leaving variable
            |> pairs -> argmin_set(((i,j),) -> j, pairs) # Type-II admissible (p,j) with minimal j
            |> ((i,j),) -> j) # extract entering variable from the single remaining pair
    else
        @assert !dual_inconsistent_certificate(BT::BasisTableau)
        p = (Type_II_admissible_pivots(BT)
            |> filter( ((i,j),) -> j==q) # filter out pairs that don't have q as entering variable
            |> pairs -> argmin_set(((i,j),) -> i, pairs) # Type-II admissible (i,q) with minimal i
            |> ((i,j),) -> i) # extract leaving variable from the single remaining pair
    end
    return (p, q)
end

function criss_cross_method(c, A, b, initial_basis)
    tableaux = Vector{BasisTableau}(undef, 0)
    push!(tableaux, BasisTableau(c, A, b, initial_basis))
    while !(is_optimal(tableaux[end]) || primal_inconsistent_certificate(tableaux[end])
                                      || dual_inconsistent_certificate(tableaux[end]))
        push!(tableaux, pivot(tableaux[end], pivotChoice_least_index_criss_cross(tableaux[end])...));
    end
    return tableaux
end
criss_cross_method(c, A, b) = criss_cross_method(c, A, b, 1:length(b))