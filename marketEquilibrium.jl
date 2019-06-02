
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
    @constraint(m,  0 <= sum(c.coef'*c.vars) + c.rhs)
    @constraint(m, sum(c.coef'*c.vars) + c.rhs <= M * slack)
    nothing
end
function add_equilibrium_constraint(m, c::LowerOrEqualThanEquilibriumConstraint, slack, M::Float64)
    @constraint(m, 0 <= -sum(c.coef'*c.vars) - c.rhs)
    @constraint(m, -sum(c.coef'*c.vars) -c.rhs <= M * slack )
    nothing
end
function add_equilibrium_constraint(m, c::FreeEquilibriumConstraint, slack, M::Float64)
    @variable(m, slack_free[1:2], Bin)

    # is in lower limit or upper limit
    # @constraint(m, c.coef'*c.vars == M *(slack_free[1]) - M * slack_free[2])

    @constraint(m, c.coef'*c.vars >= M *(slack_free[1]) - M *(slack_free[2]) - M * (1- slack)  )
    @constraint(m, c.coef'*c.vars <= M *(slack_free[1]) - M *(slack_free[2])  + M * (1 - slack) )

    # is in between
    # @constraint(m, sum(c.coef'*c.vars) >= -M *  slack)
    # @constraint(m, sum(c.coef'*c.vars) <=  M *  slack)

    # condition for fee variable slacks
    @constraint(m, sum(slack_free) == 1 - slack)

    nothing
end
function add_equilibrium_constraint(m, c::EqualEquilibriumConstraint, slack, M::Float64)
    @constraint(m, sum(c.coef'*c.vars) >= c.rhs - M * (1-slack))
    @constraint(m, sum(c.coef'*c.vars) <= c.rhs + M * (1-slack))
    nothing
end

using JuMP, Cbc#, GLPK, Ipopt
# using Gurobi
# using Plots
# plotlyjs()

M = 10^6

function normal()
    
    M = 10^6
    
    m = Model(solver = CbcSolver())
    
    @variable(m, δ1, Bin)
    @variable(m, δ2, Bin)
    @variable(m, δ3, Bin)
    @variable(m, δ4, Bin)
    @variable(m, d >= 0)
    @variable(m, q1 >= 0)
    @variable(m, q2 >= 0)
    @variable(m, p)
    
    @constraint(m, d <= M*δ1)
    @constraint(m, q1 <= M*δ2)
    @constraint(m, q2 <= M*δ3)
    @constraint(m, p <= M*δ4)
    
    @constraint(m, (100 - d - p) >= M*(1 - δ1))
    @constraint(m, (p - 2*q1) >= M*(1 - δ2))
    @constraint(m, (p - q2) >= M*(1 - δ3))
    @constraint(m, (d - q1 - q2) >= M*(1 - δ4))
    
    @constraint(m, (100 - d - p) <= 0)
    @constraint(m, (p - 2*q1) <= 0)
    @constraint(m, (p - q2) <= 0)
    @constraint(m, (d - q1 - q2) == 0)
    
    @objective(m, Min, d - q1 - q2)
    print(m)
    status = solve(m)
    
    println(getvalue(d))
    println(getvalue(q1))
    println(getvalue(q2))
    println(getvalue(p))
end
function novo()

    m = Model(solver = CbcSolver())
    
    @variable(m, δ1, Bin)
    @variable(m, δ2, Bin)
    @variable(m, δ3, Bin)
    @variable(m, δ4, Bin)
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
    cl = EqualEquilibriumConstraint([d, q1, q2], [1, -1, 1], 0)
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
novo()
# normal()