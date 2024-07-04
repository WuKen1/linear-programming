# naive implementation for converting a general form expression for a polyhedron (LP) to standard form
using LinearAlgebra
mutable struct Triplet; leq; geq; eq; end
struct StandardFormLP
    #c::Vector{<:Number} # cost vector
    A::Matrix{<:Number}
    b::Vector{<:Number}
    n::Integer # (original) number of variables
    m_leq::Integer # number of slack variables/less-than-equal constraints
    m_geq::Integer # number of surplus variables/greater-than-equal constraints
    m_eq::Integer # of equality constraints
    negative_var_indices#::Vector{<:Integer} # (with respect to the original variable indices)
    unconstrained_sign_indices#::Vector{<:Integer} # (with respect to the original variable indices)
end
# i.e. Ax = b and x ≥ 0 (both componentwise)
function to_std_form(A_leq::Matrix{<:Number}, b_leq::Vector{<:Number},
                     A_geq::Matrix{<:Number}, b_geq::Vector{<:Number},
                      A_eq::Matrix{<:Number},  b_eq::Vector{<:Number})

    A = Triplet(copy.([A_leq, A_geq, A_eq])...) # copy to avoid mutating user's input

    # For convenience. Relies on closure to ensure it remains up-to-date.
    all_A() = [A.leq, A.geq, A.eq] # Try switching to all_A = [A.leq, A.geq, A.eq] or all_A() = [A_leq, A_geq, A_eq] to see what can go wrong...
    function map_and_assign_to_all_A(f)
        A.leq, A.geq, A.eq = f.(all_A()); return nothing
    end

    @assert any(!isempty, all_A()) # make sure that there are constraints at all
    n = all_A() |> filter(!isempty) |> first |> size |> z->z[2] # number of variables
    [@assert isempty(A) || size(A)[2] == n for A in all_A()]

    # check that each A has consistent dimensions with its respective b
    A_b_pairs = ((A.leq,b_leq), (A.geq,b_geq), (A.eq,b_eq))
    [@assert size(A)[1] == length(b) for (A,b) in A_b_pairs]        
    m_leq = length(b_leq);  m_geq = length(b_geq);  m_eq = length(b_eq);

    ### identify variables that explicitly MUST be negative or nonnegative respectively ###
    #=
        If a given inequality only involves one variable (say, x_j), we have 8 cases:
        Case #1:  𝛜x ≥ 𝛅 ⟹ x ≥ 𝛅/𝛜 ≥ 0               Case #5:  𝛜x ≤  𝞭 (reduces to case #4)
        Case #2: -𝛜x ≥ 𝞭 ⟹ x ≤ -𝞭/𝛜 < 0              Case #6: -𝛜x ≤  𝞭 (reduces to case #3)
        Case #3:  𝛜x ≥ -𝞭 (inconclusive)              Case #7:  𝛜x ≤ -𝞭 (reduces to case #2)
        Case #4: -𝛜x ≥ -𝞭 ⟹ x ≤ 𝞭/𝛜 (inconclusive)   Case #8: -𝛜x ≤ -𝞭 (reduces to case #1)
        (𝛜 > 0 and  𝛅 ≥ 0)
    =# 

    function e(n, k) # k-th standard basis vector in ℝ^n
        @assert 1 ≤ k ≤ n
        return Int.((1:n) .== k)
    end
    function are_scalar_multiples(v,w)
        @assert length(v) == length(w);
        return rank([v w]) < 2
    end

    function 𝝋(i, j, ∧, ∨) # the use of symbols "∧", "∨" are to allow for infix notation
        @assert 1 ≤ j ≤ n
        # x_k uninvolved in ith [less-than/greater-than-or-equal-to] constraint for k ≠ j
        if i ≤ m_geq && are_scalar_multiples(A.geq[i,1:end], e(n, j))
            if A.geq[i,j] ∧ 0 && b_geq[i] ≥ 0 # ⟹ 𝛅/𝛜 ≥ 0 or -𝞭/𝛜 < 0
                return true
            end
        end
        if i ≤ m_leq && are_scalar_multiples(A.leq[i,1:end], e(n, j))
            if A.leq[i,j] ∨ 0 && b_leq[i] ≤ 0 # ⟹ 𝛅/𝛜 ≥ 0 or -𝞭/𝛜 < 0
                return true
            end        
        end
        return false
    end
    # xⱼ ≥ bᵢ/aᵢⱼ ≥ 0    Case #1: 𝛜x ≥ 𝛅, or Case #8: -𝛜x ≤ -𝞭, where -𝛜 = aᵢⱼ, -𝞭 = bᵢ
    xⱼ_geq_bᵢ➗aᵢⱼ_geq_0(i, j) = 𝝋(i, j, >, <)
    # xⱼ ≤ bᵢ/aᵢⱼ ≤ 0    Case #2: -𝛜x ≥ 𝞭, or Case #7: 𝛜x ≤ -𝞭
    xⱼ_leq_bᵢ➗aᵢⱼ_leq_0(i, j) = 𝝋(i, j, <, >)

    # identify variables that explicitly MUST be negative or nonnegative respectively
    nonnegative_indices = []; negative_indices = []
    for j = 1:n
        if any([xⱼ_geq_bᵢ➗aᵢⱼ_geq_0(i,j) for i=1:max(m_leq, m_geq, m_eq)])
            push!(nonnegative_indices, j)
        end
        if any([xⱼ_leq_bᵢ➗aᵢⱼ_leq_0(i, j) for i=1:max(m_leq, m_geq, m_eq)])
            push!(negative_indices, j)
        end
    end
    if ~isempty(nonnegative_indices ∩ negative_indices)
        msg = "the following variables have been identified as " *
                "having conflicting (non)negativity constraints:\n" *
                string(nonnegative_indices ∩ negative_indices)
        error(msg)
    end
    # remaining variables are presumably unconstrained in sign
    unconstrained_sign_indices = setdiff(1:n, nonnegative_indices ∪ negative_indices)

    # Fixes/corrects the number of columns in unused/empty matrices. This needs to be done so that:
    # i) we don't mistakenly "attempt to access 0x0 Matrix" potentially when scaling columns by -1, and
    # ii) to avoid "DimensionMismatch: number of columns of each array must match" error when using vcat.
    fix_dim_of_empty(M) = ~isempty(M) ? M : reshape(M, 0, n + length(unconstrained_sign_indices))
    map_and_assign_to_all_A(fix_dim_of_empty)
    
    # replace any negatively constrained variable xⱼ by xⱼ' ≔ -xⱼ, then multiply all of the
    # corresponding coeffecients aᵢⱼ by -1 accordingly to compensate
    [A[:,j] = -A[:,j] for A in all_A(), j in negative_indices]

    # Replaces each unconstrained (i.e. potentially negative or nonnegative) variable xⱼ
    # with the difference of two nonnegative artifical variables xⱼ⁺ - xⱼ⁻.
    #
    # Coeffecients are duplicated but with -1 since aᵢⱼxⱼ = aᵢⱼ(xⱼ⁺-xⱼ⁻) = aᵢⱼxⱼ⁺ + (-aᵢⱼ)xⱼ⁻.
    #
    # The column of coeffecients corresponding to xⱼ⁻ are inserted between the current jth column
    # (j+1)th column (if it even exists), starting with the largest j ∈ unconstrained_sign_indices,
    # and ending with the smallest. (done in this order for ease of insertion)
    function diff_of_nonnegative(M)
        if isempty(M); return M; end
        for j in sort(unconstrained_sign_indices, rev=true) # rev accounts for changing # of columns
            M = hcat(M[:,1:j-1], M[:,j], -M[:,j], M[:,j+1:end])
        end
        return M
    end
    map_and_assign_to_all_A(diff_of_nonnegative)
    
    # creates a block matrix that serves as a row of the end result block matrix
    block_A = [
                A.leq              1I                zeros(Integer, m_leq, m_geq);
                A.geq  zeros(Integer, m_geq, m_leq)             -I;
                A.eq   zeros(Integer, m_eq, m_leq)   zeros(Integer, m_eq, m_geq)
              ]
    b = vcat(b_leq, b_geq, b_eq)
    recovery_info = (n, m_leq, m_geq, m_eq, negative_indices, sort(unconstrained_sign_indices))
    return StandardFormLP(block_A, b, recovery_info...)
end

struct AugmentedMatrix
    A::Matrix{<:Number}
    b::Vector{<:Number}
end # TO-DO: maybe add related methods supporting row operations?
AugmentedMatrix() = AugmentedMatrix(Matrix{Number}(undef, 0, 0), Vector{Number}(undef, 0))
AugmentedMatrix(M) = isempty(M) ? AugmentedMatrix() : AugmentedMatrix(M[:,1:end-1], M[:,end])

function to_std_form(; aug_leq::AugmentedMatrix=AugmentedMatrix(),
                       aug_geq::AugmentedMatrix=AugmentedMatrix(),
                        aug_eq::AugmentedMatrix=AugmentedMatrix())
    return to_std_form(aug_leq.A, aug_leq.b, aug_geq.A, aug_geq.b, aug_eq.A, aug_eq.b)
end

### test example below ###
AUG_LEQ = AugmentedMatrix(
[ 1  0  0 -7;   #   x            ≤  -7   (note that x must be negative)
 -2  0  5  8;   # -2x      + 5z  ≤   8
  3 -4  6  9]); #  3x - 4y + 6z  ≤   9

AUG_GEQ = AugmentedMatrix(
[10 0 11 -13    # 10x + 11z  ≥ -13
  0 0 12  14]); #       12z  ≥  14   (note that z must be positive)

#AUG_EQ = AugmentedMatrix([0 0 1 -3]); # z = 3
AUG_EQ = AugmentedMatrix();

q = to_std_form(aug_leq=AUG_LEQ, aug_geq=AUG_GEQ, aug_eq=AUG_EQ);

println("AUG_LEQ.A:")
display(AUG_LEQ.A); println()
println("AUG_GEQ.A:")
display(AUG_GEQ.A); println()
println("q.A:")
display(q.A)

# debug_print(x) = printstyled("$x\n"; color = :red, bold=true, italic=true, underline=true);

#num_of_row_in_block_matrix = m_leq + m_geq + m_eq
#num_of_col_in_block_matrix = (n + length(unconstrained_sign_indices)) + (m_leq + m_geq)
#block_A = Matrix{Number}(undef, num_of_row_in_block_matrix, num_of_col_in_block_matrix)



### UNTESTED; need to rewrite ###
#=
function revert_sol_from_std_form(std::Std_form, sol)
    s = std; x = sol[2s.m]
    for j in s.unconstrained_sign_indices
        x[j] = x[j] - x[j+1]
        x = x[1:end .!= j+1]
    end
    x[s.negative_var_indices] = -x[s.negative_var_indices]
    x = x[1:(end-(s.m_leq + s.m_geq))]
end
=#
