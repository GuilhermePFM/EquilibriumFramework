
using JuMP
using Cbc

abstract type EquilibriumConstraint end


struct GreaterOrEqualThanEquilibriumConstraint <:EquilibriumConstraint
    vars::Vector{Any}
    coef::Vector{Float64}
    rhs::Float64
end

struct LowerOrEqualThanEquilibriumConstraint <:EquilibriumConstraint
    vars::Vector{Any}
    coef::Vector{Float64}
    rhs::Float64
end

struct FreeEquilibriumConstraint <:EquilibriumConstraint
    vars::Vector{Any}
    coef::Vector{Float64}
    rhs::Float64
end

struct IntervalEquilibriumConstraint <:EquilibriumConstraint
    vars::Vector{Any}
    coef::Vector{Float64}
    rhs::Float64
    LB::Float64
    UB::Float64
end

struct EqualEquilibriumConstraint <:EquilibriumConstraint
    vars::Vector{Any}
    coef::Vector{Float64}
    rhs::Float64
end

struct ComplementarityEquilibriumConstraint
    cstr1::EquilibriumConstraint
    cstr2::EquilibriumConstraint
    M::Float64
    ComplementarityEquilibriumConstraint(cstr1, cstr2) = new(cstr1, cstr2, 1e8)
end

function add_equilibrium_constraint(m, c::ComplementarityEquilibriumConstraint)
    slack = @variable(m, [1], Bin)
    
    # notacao:
    # - se a restricao ta desativada, a slack que a funcao recebe eh igual a 1
    # - se a restricao ta ativada, a slack que a funcao recebe eh igual a 0
    add_equilibrium_constraint(m, c.cstr1, slack[1], c.M)
    add_equilibrium_constraint(m, c.cstr2, (1-slack[1]), c.M)

    nothing
end

function add_equilibrium_constraint(m, c::GreaterOrEqualThanEquilibriumConstraint, slack, M::Float64)
    # complementarity greater or equal constraints can be modeled using Fortuny-Amat McCarl Linearization [Complementarity modeling in energy markets]: 
    # 0 ≤ x + RHS ⊥ y + RHS ≥ 0
    # equivalent to:
    # 0 ≤ x + RHS ≤ Mu
    # 0 ≤ y + RHS ≤ M(1 − u)
    @constraint(m,  0 <= sum(c.coef'*c.vars) + c.rhs)
    @constraint(m, sum(c.coef'*c.vars) + c.rhs <= M * slack)
    nothing
end

function add_equilibrium_constraint(m, c::LowerOrEqualThanEquilibriumConstraint, slack, M::Float64)
    # complementarity lower or equal constraints can be modeled using Fortuny-Amat McCarl Linearization [Complementarity modeling in energy markets]: 
    # 0 ≥ x + RHS ⊥ y + RHS ≥ 0
    # 0 ≤ -x -RHS ⊥ y + RHS ≥ 0
    # equivalent to:
    # 0 ≤ -x - RHS ≤ Mu
    # 0 ≤ y + RHS ≤ M(1 − u)
    @constraint(m, 0 <= -sum(c.coef'*c.vars) - c.rhs)
    @constraint(m, -sum(c.coef'*c.vars) - c.rhs <= M * slack )
    nothing
end

function add_equilibrium_constraint(m, c::FreeEquilibriumConstraint, slack, M::Float64)
    # complementarity free constraints have two options: 
    # 1. be binding at UB or LB (Inf in this case)
    # 2. be relaxed in the defined interval (be free)
    # free constraints have 2 ways to be binding (Modeling Mathematical Programs with Equilibrium Constraints in Pyomo pg 9:
    # 1. x = UB (Inf)
    # 2. x = LB (-Inf)
    active_aux = @variable(m, [1:2], Bin)

    @constraint(m, c.coef'*c.vars >= M *(active_aux[1]) - M *(active_aux[2]) - M * slack  )
    @constraint(m, c.coef'*c.vars <= M *(active_aux[1]) - M *(active_aux[2]) + M * slack )

    # condition for fee variable slacks
    @constraint(m, sum(active_aux) == 1 - slack)

    nothing
end

function add_equilibrium_constraint(m, c::IntervalEquilibriumConstraint, slack, M::Float64)
    # complementarity interval constraints have two options: 
    # 1. be binding at UB or LB
    # 2. be relaxed in the defined interval
    # free constraints have 2 ways to be binding (Modeling Mathematical Programs with Equilibrium Constraints in Pyomo pg 9:
    # 1. x = UB 
    # 2. x = LB 
    active_aux = @variable(m, [1:2], Bin)

    @constraint(m, c.coef'*c.vars >= c.UB *(active_aux[1]) + c.LB *(active_aux[2]) + c.LB * slack  )
    @constraint(m, c.coef'*c.vars <= c.UB *(active_aux[1]) + c.LB *(active_aux[2]) + c.UB * slack )

    # condition for fee variable slacks
    @constraint(m, sum(active_aux) == 1 - slack)

    nothing
end

function add_equilibrium_constraint(m, c::EqualEquilibriumConstraint, slack, M::Float64)
    # complementarity constraints have two options: 
    # 1. be binding at equality
    # 2. be relaxed in the defined interval
    # equality constraints have 2 ways to be relaxed (Modeling Mathematical Programs with Equilibrium Constraints in Pyomo pg 9:
    # 1. g(x) >=0
    # 2. g(x) <=0
    @variable(m, slack_free_eq[1:2], Bin)
    @constraint(m, sum(c.coef'*c.vars) >= c.rhs - M * slack_free_eq[1])
    @constraint(m, sum(c.coef'*c.vars) <= c.rhs + M * slack_free_eq[2])
    @constraint(m, sum(slack_free_eq) == slack)
    nothing
end

using JuMP, Cbc

function example_1()

    m = Model(solver = CbcSolver())
    
 
    @variable(m, q1>=0)
    @variable(m, q2>=0)
    @variable(m, d >=0)
    @variable(m, p)
    
    # Consumers
    c = LowerOrEqualThanEquilibriumConstraint([d, p], [-1, -1], 100)
    cdual = GreaterOrEqualThanEquilibriumConstraint([d], [1], 0)
    consumer = ComplementarityEquilibriumConstraint(c, cdual)
    add_equilibrium_constraint(m, consumer)
    
    # Firm 1
    f1 = LowerOrEqualThanEquilibriumConstraint([p, q1], [1, -2], 0)
    f1dual = GreaterOrEqualThanEquilibriumConstraint([q1], [1], 0)
    firm1 = ComplementarityEquilibriumConstraint(f1, f1dual)
    add_equilibrium_constraint(m, firm1)
    
    # Firm 2
    f2 = LowerOrEqualThanEquilibriumConstraint([p, q2], [1, -1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([q2], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m, firm2)
    
    # Market Clearing
    cl = EqualEquilibriumConstraint([d, q1, q2], [1, -1, -1], 0)
    cldual = FreeEquilibriumConstraint([p], [1], 0)
    clearing = ComplementarityEquilibriumConstraint(cl, cldual)
    add_equilibrium_constraint(m, clearing)
    # @constraint(m, d - q1 - q2 == 0)
    
    @objective(m, Min, d - q1 - q2 )
    println(m)
    status = solve(m)
    println(status)
    println("d=" , getvalue(d))
    println("q1=" ,getvalue(q1))
    println("q2=" ,getvalue(q2))
    println("p=" ,getvalue(p))
    writeLP(m, "lp", genericnames=false)
end

function example_2()

    m = Model(solver = CbcSolver())
 
    @variable(m, q1)
    @variable(m, q2)
    @variable(m, d)
    @variable(m, p)
    @variable(m, pc02)
    @variable(m, e1)
    @variable(m, lambda1)
    @variable(m, e2)
    @variable(m, lambda2)
    
    # Consumers
    c = LowerOrEqualThanEquilibriumConstraint([d, p], [-1, -1], 100)
    cdual = GreaterOrEqualThanEquilibriumConstraint([d], [1], 0)
    consumer = ComplementarityEquilibriumConstraint(c, cdual)
    add_equilibrium_constraint(m, consumer)
    
    # Firm 1
    # utility function
    f1 = LowerOrEqualThanEquilibriumConstraint([p, q1, lambda1], [1, -2, -1/2], 0)
    f1dual = GreaterOrEqualThanEquilibriumConstraint([q1], [1], 0)
    firm1 = ComplementarityEquilibriumConstraint(f1, f1dual)
    add_equilibrium_constraint(m, firm1)
    # co2 constraint
    e1_cstr_1 = LowerOrEqualThanEquilibriumConstraint([pc02, lambda1], [-1, 1], 0)
    e1_cstr_2 = GreaterOrEqualThanEquilibriumConstraint([e1], [1], 0)
    c02_1_contraint = ComplementarityEquilibriumConstraint(e1_cstr_1, e1_cstr_2)
    add_equilibrium_constraint(m, c02_1_contraint)
    # lambda constraint
    lambda1_cstr_1 = LowerOrEqualThanEquilibriumConstraint([q1, e1], [1/2, -1], 0)
    lambda1_cstr_2 = GreaterOrEqualThanEquilibriumConstraint([lambda1], [1], 0)
    c02_2_contraint = ComplementarityEquilibriumConstraint(lambda1_cstr_1, lambda1_cstr_2)
    add_equilibrium_constraint(m, c02_2_contraint)

    # Firm 2
    f2 = LowerOrEqualThanEquilibriumConstraint([p, q2, lambda2], [1, -1, -1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([q2], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m, firm2)
    # co2 constraint
    e2_cstr_1 = LowerOrEqualThanEquilibriumConstraint([pc02, lambda2], [-1, 1], 0)
    e2_cstr_2 = GreaterOrEqualThanEquilibriumConstraint([e2], [1], 0)
    c02_3_contraint = ComplementarityEquilibriumConstraint(e2_cstr_1, e2_cstr_2)
    add_equilibrium_constraint(m, c02_3_contraint)
    # lambda constraint
    lambda2_cstr_1 = LowerOrEqualThanEquilibriumConstraint([q2, e2], [1, -1], 0)
    lambda2_cstr_2 = GreaterOrEqualThanEquilibriumConstraint([lambda2], [1], 0)
    c02_4_contraint = ComplementarityEquilibriumConstraint(lambda2_cstr_1, lambda2_cstr_2)
    add_equilibrium_constraint(m, c02_4_contraint)
    
    # Market Clearing
    cl = EqualEquilibriumConstraint([d, q1, q2], [1, -1, -1], 0)
    cldual = FreeEquilibriumConstraint([p], [1], 0)
    clearing = ComplementarityEquilibriumConstraint(cl, cldual)
    add_equilibrium_constraint(m, clearing)
    
    # M.C.C. of CO2 permits
    co2_permit_1 = LowerOrEqualThanEquilibriumConstraint([e1, e2], [1, 1], -30)
    co2_permit_2 = GreaterOrEqualThanEquilibriumConstraint([pc02], [1], 0)
    c02_5_contraint = ComplementarityEquilibriumConstraint(co2_permit_1, co2_permit_2)
    add_equilibrium_constraint(m, c02_5_contraint)
    
    @objective(m, Min, d - q1 - q2 )
    println(m)
    status = solve(m)
    println(status)
    println("d=" , getvalue(d))
    println("q1=" ,getvalue(q1))
    println("q2=" ,getvalue(q2))
    println("p=" ,getvalue(p))
    println("e1=" ,getvalue(e1))
    println("e2=" ,getvalue(e2))
    println("lambda1=" ,getvalue(lambda1))
    println("lambda2=" ,getvalue(lambda2))
    println("pc02=" ,getvalue(pc02))
    writeLP(m, "lp", genericnames=false)
end


# run
# example_1()
example_2()