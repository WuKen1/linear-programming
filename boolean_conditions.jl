primal_is_feasible(BT::BasisTableau) = all(BT.btab[2:end, 1]  .|>  ≥(0)) # i.e. x .≥ 0
dual_is_feasible(BT::BasisTableau) = all(BT.btab[1, 2:end]  .|>  ≤(0)) # i.e., s .≥ 0

is_optimal(BT::BasisTableau) = primal_is_feasible(BT)  &&  dual_is_feasible(BT)

# certificate: if it returns true, the primal is infeasible;
#              if it returns false, inconclusive
primal_inconsistent_certificate(BT::BasisTableau) =
    any([BT.btab[1+k,1] < 0 && all(BT.btab[1+k,2:end] .≥ 0)
        for m = [length(BT.b)]
        for k in 1:m])
dual_inconsistent_certificate(BT::BasisTableau) =
    any([BT.btab[1,1+ℓ] > 0 && all(BT.btab[2:end,1+ℓ] .≤ 0)
        for n = [length(BT.c)]
        for ℓ in 1:n]) # (ℓ is lowercase L)

primal_unbounded_certificate(BT::BasisTableau) =
    dual_inconsistent_certificate(BT) && primal_is_feasible(BT)
dual_unbounded_certificate(BT::BasisTableau) =
    primal_inconsistent_certificate(BT) && dual_is_feasible(BT)