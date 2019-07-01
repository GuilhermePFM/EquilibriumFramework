using JuMP
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

struct marketEquilibriumModel
    m::JuMP.Model
    cstrs::Array{ComplementarityEquilibriumConstraint}
    marketEquilibriumModel() = new(Model(solver = CbcSolver()), ComplementarityEquilibriumConstraint[])
end

function add_equilibrium_constraint(m::marketEquilibriumModel, c::ComplementarityEquilibriumConstraint)
    push!(m.cstrs, c)
    nothing
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
    
    @constraint(m, sum( c.coef'*c.vars ) >= M *(active_aux[1]) - M *(active_aux[2]) - M * slack  )
    @constraint(m, sum( c.coef'*c.vars ) <= M *(active_aux[1]) - M *(active_aux[2]) + M * slack )

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