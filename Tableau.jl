struct BasisTableau
    c::Vector; A::Matrix; b::Vector # LP parameters
    #=
        indices_basic and indices_nonbasic are index sets, but can also be
        viewed as functions that that map a "raw" index to a variable index
    =#
    indices_basic::Array; indices_nonbasic::Array
    btab::Matrix # the entire tableau
    #T::Matrix # T := [I  (A_B)⁻¹∗A_N]
end
function BasisTableau(c, A, b, indices_basic)
    n = size(A)[2]; indices_nonbasic = setdiff(1:n, indices_basic)
    return BasisTableau(c, A, b, indices_basic, indices_nonbasic)
end
function BasisTableau(c, A, b, indices_basic, indices_nonbasic)
    m,n = size(A)
    I_B, I_N = copy.([indices_basic, indices_nonbasic])
    
    @assert length(b) == m;   @assert length(c) == n
    @assert length(I_B) == m "length(I_B) == $(length(I_B)) != $m == m"; @assert m ≤ n

    A_B = A[:, I_B];           A_N = A[:, I_N]
    c_B = c[I_B];              c_N = c[I_N]
    x_B = A_B^-1 * b;          x_N = zeros(n-m)
    s_B = zeros(size(A)[1]);   s_N = (c_N' - c_B' * A_B^-1 * A_N)'

    btab = [c_B'*(A_B^-1)*b   -s_B'   -s_N';
            x_B                I       A_B^-1*A_N]
    #T = [I  A_B^-1*A_N]
    return BasisTableau(c, A, b, sort(indices_basic), sort(indices_nonbasic), btab)
end
########################################################################################

#=
    We make a distinction between "variable" indices, "raw" indices, and "tableau" indices.
    Let [k] := {1, ..., k} for any k ∈ ℕ, and let B and N denote the set of basic and nonbasic
    variables respectively.
    Define 𝓑 : [m] ⟶ B and 𝓝 : [n-m] ⟶ N to be the bijective functions that determine them
    ordering of the indices within B and N. For the sake of readability, we will follow the
    convention that both 𝓑 and 𝓝 are increasing functions, which we will enforce and maintain
    as an invariant.
=#
# type `\bscrB`` for 𝓑; use similar commands for 𝓝, 𝓣
#, and 𝓡

# maps a variable index to its "raw" index w/ respect to said index set
function index_raw(BT::BasisTableau, index)
    @assert index ∈ BT.indices_basic ∪ BT.indices_nonbasic
    if index ∈ BT.indices_basic
        n = length(BT.indices_basic)
        return Dict([(BT.indices_basic[i],i) for i = 1:n])[index]
    else # index is nonbasic
        n = length(BT.indices_nonbasic)
        return Dict([(BT.indices_nonbasic[i],i) for i = 1:n])[index]
    end
end
𝓑(BT::BasisTableau, index) = BT.indices_basic[index]
function 𝓑⁻¹(BT::BasisTableau, index)
    @assert index ∈ BT.indices_basic; return index_raw(BT, index)
end
𝓝(BT::BasisTableau, variableIndex) = BT.indices_nonbasic[variableIndex]
function 𝓝⁻¹(BT::BasisTableau, index)
    @assert index ∈ BT.indices_nonbasic; return index_raw(BT, index)
end
# maps a variable index to its "T index", i.e., w/ respect to T := [I  (A_B)⁻¹∗A_N]
function index_T(BT::BasisTableau, index) # (also ordering of x and s)
    if index ∈ BT.indices_basic
        return 𝓑⁻¹(BT, index)
    else m = length(BT.b);
        return m + 𝓝⁻¹(BT, index)
    end
end
𝓣(BT::BasisTableau, index) = index_T(BT, index) # alias for index_T
# maps a variable index to its tableau index, i.e., w/ respect to the entire tableau
index_tableau(BT::BasisTableau, index) = 1 + index_T(BT, index)
#𝓡(BT::BasisTableau, index) = index_tableau(BT, index) # alias for index_tableau

############### convenience functions for relevant named entities ###############
#=
    In principle, one can avoid explicitly taking the inverse of A_B in the
    computation of x_B and s_N by simply reading them off of the tableau, which 
    is updated via row operations instead of matrix multiplication with inverses.
    These are mainly here for debugging purposes (and again, convenience).
=#
A_B(BT) = BT.A[:, BT.indices_basic];  A_N(BT) = BT.A[:, BT.indices_nonbasic]
c_B(BT) = BT.c[BT.indices_basic];     c_N(BT) = BT.c[BT.indices_nonbasic]
x_B(BT) = A_B(BT)^-1 * BT.b;          x_N(BT) = zeros([n-m for (m,n)=[size(BT.A)]][1])
s_B(BT) = zeros(size(BT.b));          s_N(BT) = (c_N(BT)' - c_B(BT)' * A_B(BT)^-1 * A_N(BT))'
function x(BT::BasisTableau)
    n = size(BT.c); x = zeros(n);
    [x[i] = x_B(BT)[𝓑⁻¹(BT, i)] for i ∈ BT.indices_basic]; return x
end
function s(BT::BasisTableau)
    n = size(BT.c); s = zeros(n);
    [s[i] = s_N(BT)[𝓝⁻¹(BT, i)] for i ∈ BT.indices_nonbasic]; return s
end
y(BT) = A_B(BT)'^-1 * c_B(BT)

#=
    If we cut off the 1st row of the tableau, the resulting augmented matrix is shorthand for:
   
        A_B^-1 * b = [I (A_B^-1*A_N)][x_B; x_N]

    (Here, [A B] := hcat(A,B), while [A; B] := vcat(A,B).)
    We can encode the cumulative row operations needed to turn a given column
    into a standard basis vector by a matrix R, and we can encode the row swaps
    needed to "untangle" x_B and x_N (as well as sort them in ascending order)
    by a permutation matrix P. Recall that the inverse of a permutation matrix is itself.
    We can put the tableau into "standard form" for our new x_B and x_N:

        [c_B()'*(A_B()^-1)*b  s_B()'  -s_N()';
         x_B()                I        A_B()^-1*A_N()]

    by multiplying both sides of the equation on the left by R, while inserting P*P^-1
    in-between [I (A_B^-1*A_N)] and [x_B; x_N] we get:
        
        R * A_B^-1 * b = R * [I (A_B^-1*A_N)] * (P^-1 * P) * [x_B; x_N]
                       = R * [I (A_B^-1*A_N)] * (P' * P) * [x_B; x_N]
                       = (R * [I (A_B^-1*A_N)] * P') * (P * [x_B; x_N])

    P is used to permute the rows of x, while P⁻¹ = Pᵀ permutes
    the columns of [I (A_B^-1*A_N)].

    TO-DO: update this description, particularly regarding cumulative row operation matrix R
=#

# TO-DO: fix pivot permutations; cheat and re-compute tableau from scratch for now
function pivot(BT::BasisTableau, leavingVariable, enteringVariable) # (inputs are "variable" indices)
    I_B, I_N = copy.([BT.indices_basic, BT.indices_nonbasic])
    @assert leavingVariable ∈ I_B
    @assert enteringVariable ∈ I_N
    btab = copy(BT.btab);  m, n = size(BT.A)
    #=
    Let x_p and x_q denote the leaving and entering variables respectively. We perform row
    operations so as to transform the column of the tableau corresponding to the tableau
    index of q into [0; e_(p; m)], where e_(p; m) is the p-th standard basis vector in R^m. =#
    
    transform_column_into_std_basis_vector!(btab, index_tableau(BT,leavingVariable),
                                                  index_tableau(BT,enteringVariable))
    #=
    In order to reorganize x into [x_B_new; x_N_new] for our new (ordered)
    sets of basic/nonbasic variables (and likewise for s), we first need to swap the rows of
    [x_B_old; x_N_old] corresponding to the leaving and entering variables. =#

    #println("B = $(BT.indices_basic)   N = $(BT.indices_nonbasic)") # debug
    #println("leaving: x_$leavingVariable") # debug
    #println("entering: x_$enteringVariable") # debug
    #println("𝓣(x_$leavingVariable) = $(𝓣(BT, leavingVariable))") # debug
    #println("𝓣(x_$enteringVariable) = $(𝓣(BT, enteringVariable))") # debug
    P_swap = row_swap!(I[1:n, 1:n], 𝓣(BT, leavingVariable), 𝓣(BT, enteringVariable))
    #=
    We then perform row swaps to sort x_B_intermediate and
    X_N_intermediate. At this point, row swaps should only be
    performed *within* the separate blocks, not *between* the two blocks. =#
    I_B[𝓑⁻¹(BT, leavingVariable)]  = enteringVariable
    I_N[𝓝⁻¹(BT, enteringVariable)] = leavingVariable
    #println("pivoting from $(BT.indices_basic) to $I_B") # debug
        # Comment: we do the updates for I_B and I_N directly instead of via
        #          unions and & set differences because we want the explicit
        #          matrix representation of the permutations.
    #P_B, P_N = ascending_sort_permutation_matrix.([I_B, I_N]) # note the broadcasting
    #P_B = ascending_sort_permutation_matrix(I_B) # note the broadcasting
    #P_N = ascending_sort_permutation_matrix(I_N) # note the broadcasting
    #P = direct_sum(P_B, P_N) * P_swap
    #P = P_swap
    #println("P_swap * [I_B; I_N] ="); display(P_swap * [BT.indices_basic; BT.indices_nonbasic])
    #println("P_B * (P_swap * [I_B; I_N])_B ="); display(P_B * (P_swap * [BT.indices_basic; BT.indices_nonbasic])[1:m])
    #println("P_N * (P_swap * [I_B; I_N])_N ="); display(P_N * (P_swap * [BT.indices_basic; BT.indices_nonbasic])[(m+1):end])
    #println("P_N * (P_swap * [I_B; I_N])_N ="); display(P * P_swap * [BT.indices_basic; BT.indices_nonbasic])

    #println("P_swap ="); display(P_swap); println()
    #println("P_B ="); display(P_B); println()
    #println("P_N ="); display(P_N); println()

    #println("before permuting:"); display(btab[:, 2:end]) # debug
    #btab[:, 2:end] = btab[:, 2:end] * P'
    #println("after permuting:"); display(btab[:, 2:end]); println() # debug
    sort!.([I_B, I_N]);
    BT_new = BasisTableau(BT.c, BT.A, BT.b, I_B, I_N, btab)

    ### debugging ###
    #println("pivoting from $(BT.indices_basic) to $I_B") # debug
    BT_new_naive = BasisTableau(BT.c, BT.A, BT.b, I_B)

    #BT_new.btab[1,:] = BT_new_naive.btab[1,:] # cheat on 1st row
    #BT_new.btab[2:end,1] = BT_new_naive.btab[2:end,1] # cheat on x_B

    #debug_equality("c_B' * A_B^-1 * b", BT_new.btab[1,1], "naive ver", BT_new_naive.btab[1,1])
    #debug_equality("x_B", BT_new.btab[2:end, 1], "naive ver", BT_new_naive.btab[2:end, 1])
    #debug_equality("-S_B", BT_new.btab[1,   1 .+ (1:m)], "naive ver", BT_new_naive.btab[1,   1 .+ (1:m)])
    #debug_equality("-S_N", BT_new.btab[1,   (1+m+1):end], "naive ver", BT_new_naive.btab[1,   (1+m+1):end])
    #debug_equality("I_m", BT_new.btab[2:end, 1 .+ (1:m)], "naive ver", BT_new_naive.btab[2:end, 1 .+ (1:m)])
    #debug_equality("A_B^-1 * A_N", BT_new.btab[2:end, (1+m+1):end], "naive ver", BT_new_naive.btab[2:end, (1+m+1):end])
    
    return BT_new_naive # debug
end

function debug_equality(name1, q1, name2, q2)
    error_message = "$name1 == $q1 != $q2 == $name2"
#    @assert isapprox(q1, q2) error_message
    
    if isapprox(q1, q2)
        println("$name1 correctly computed")
    else
        println("$name1")
        display(q1)
        println("naive ver:")
        display(q2)
        println()
    end
end

#### miscellaneous operations that aren't        ####
#### necessarily specific to linear programming: ####
function argmin_set(f, domain)
    f_min = minimum(f, domain)    
    indices = findall(x -> f(x) == f_min, domain)
    return domain[indices]
end

# returns a "block diagonal matrix" with A and B as submatrices along the diagonal
function direct_sum(A::Matrix,B::Matrix)
    m, n = size(A);   p, q = size(B);
    return [A zeros(m,q); zeros(p,n) B]
end
e(k,n) = I[1:n, k] # k-th standard basis vector in R^n
function row_swap!(A,i,j)
    A[[i,j],:] = A[[j,i],:]; return A
end
function ascending_sort_permutation_matrix(indices)
    indices = copy(indices)
    n = length(indices); P = I[1:n, 1:n]
    for i = 1:length(indices)
        j = (i-1) + argmin(indices[i:end])
        P = row_swap!(I[1:n, 1:n], i, j) * P
        indices[[i,j]] = indices[[j,i]]
    end
    return P
end
# performs row operations in order to transform A[:,j] into e(i,m), where m = size(A)[1]
function transform_column_into_std_basis_vector!(A::Matrix{T}, i, j) where T<:AbstractFloat
    #A = 1.0*A; # note: including this line causes the input matrix to not be mutated
    m,n = size(A);
    @assert 1 ≤ i ≤ m;  @assert 1 ≤ j ≤ n
    #print("i=$i, j=$j;  before:"); display(A); println()
    for k ∈ setdiff(1:m, i)
        if A[k,j] ≠ 0
            # R_k <-- R_k + (-a_ij/a_ij)R_k    for some k with a_kj ≠ 0
            A[k,:] += (-A[k,j] / A[i,j]) * A[i,:] ### no such k ###
        else
            # don't need to do anything
        end
    end
    A[i,:] *= 1/A[i,j]
    #print("After: "); display(A); println()
end

#=
# performs row operations in order to transform A[I,j] into v (or at least tries to?)
function replace_subcolumn_via_row_operations!(A::Matrix{T}, v, I::UnitRange, j) where T<:AbstractFloat
#    A = 1.0*A; # note: including this line causes the input matrix to not be mutated
    m,n = size(A);
    @assert length(I) == length(v); @assert typeof(v) <: Vector
    @assert 1 ≤ I[1] ≤ I[end] ≤ m;  @assert 1 ≤ j ≤ n
    for i = 1:length(v)
        k = findfirst(k -> k ≠ I[i] && A[k,j] ≠ 0, 1:m) # TO-DO: error handling?
        @assert k != nothing "i=$i, j=$j, k=$k,\n v = $v, I = $I,\nA = $A"    
        if A[I[i],j] == 0 && v[i] == 0
            # don't need to do anything
        elseif A[I[i],j] ≠ 0 && v[i] ≠ 0
            # R_i <-- (1/a_ij) * R_j
            A[I[i],:] *= (1 / A[I[i],j])
        elseif A[I[i],j] ≠ 0 && v[i] == 0
            # R_i <-- R_i + (-a_ij/a_kj)R_k    for some k with a_kj ≠ 0
            A[I[i],:] += (-A[I[i],j] / A[k,j]) * A[k,:] ### no such k ###
        else # A[I[i],j] == 0 && v[i] ≠ 0
            # R_i <-- R_i + (-v_i/a_kj)R_k    for some k with a_kj ≠ 0
            A[I[i],:] += (v[i] / A[k,j]) * A[k,:]
        end        
    end
    @assert A[I,j] == v
    #return A
end
=#