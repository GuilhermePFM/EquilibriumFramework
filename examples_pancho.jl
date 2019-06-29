include("marketEquilibrium")

using JuMP
using Cbc

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

function example_3()
    T = 2
    g1_0 = 0
    g2_0 = 0
    d0 = 0
    ramp_g1 = 1
    ramp_g2 = 8

    m = Model(solver = CbcSolver())
 
    # demand
    @variable(m, d[1:T+1])
    @variable(m, p[1:T+1])
    # g1
    @variable(m, g1[1:T+1])
    @variable(m, lambda1[1:T+1])
    @variable(m, alpha1[1:T+1])
    @variable(m, beta1[1:T+1])
    # g2
    @variable(m, g2[1:T+1])
    @variable(m, lambda2[1:T+1])
    @variable(m, alpha2[1:T+1])
    @variable(m, beta2[1:T+1])
    
    # @constraint(m, [t=2:T+1],alpha2[t] == 0.0)
    # @constraint(m, [t=2:T+1],alpha1[t] == 0.0)
    # @constraint(m, [t=2:T+1],beta1[t] == 0.0)
    # @constraint(m, [t=2:T+1], beta2[t] == 0.0)

    @constraint(m, g1[1] == g1_0)
    @constraint(m, g2[1] == g2_0)
    @constraint(m, d[1] == d0)
    for t in 2:(T+1)
        # Consumers
        c = LowerOrEqualThanEquilibriumConstraint([p[t], d[t]], [-1, -10/3], 20)
        cdual = GreaterOrEqualThanEquilibriumConstraint([d[t]], [1], 0)
        consumer = ComplementarityEquilibriumConstraint(c, cdual)
        add_equilibrium_constraint(m, consumer)

        c2 = EqualEquilibriumConstraint([d[t], g1[t], g2[t]], [1, -1, -1], 0)
        cdual2 = FreeEquilibriumConstraint([p[t]], [1], 0)
        consumer2 = ComplementarityEquilibriumConstraint(c2, cdual2)
        add_equilibrium_constraint(m, consumer2)
        
        # Firm 1
        # utility function
        f1 = LowerOrEqualThanEquilibriumConstraint([g1[t], p[t], lambda1[t], alpha1[t], beta1[t]], [-2, 1, -1, 1, 1], 0)
        f1dual = GreaterOrEqualThanEquilibriumConstraint([g1[t]], [1], 0)
        firm1 = ComplementarityEquilibriumConstraint(f1, f1dual)
        add_equilibrium_constraint(m, firm1)

        f1_2 = LowerOrEqualThanEquilibriumConstraint([g1[t]], [1], -3)
        f1dual_2 = GreaterOrEqualThanEquilibriumConstraint([lambda1[t]], [1], 0)
        firm1_2 = ComplementarityEquilibriumConstraint(f1_2, f1dual_2)
        add_equilibrium_constraint(m, firm1_2)
        # subida
        f1_3 = LowerOrEqualThanEquilibriumConstraint([g1[t], g1[t-1]], [1, -1], -1)
        f1dual_3 = GreaterOrEqualThanEquilibriumConstraint([alpha1[t]], [1], 0)
        firm1_3 = ComplementarityEquilibriumConstraint(f1_3, f1dual_3)
        add_equilibrium_constraint(m, firm1_3)
        # descida
        f1_4 = GreaterOrEqualThanEquilibriumConstraint([g1[t], g1[t-1]], [1, -1], 1)
        f1dual_4 = GreaterOrEqualThanEquilibriumConstraint([beta1[t]], [1], 0)
        firm1_4 = ComplementarityEquilibriumConstraint(f1_4, f1dual_4)
        add_equilibrium_constraint(m, firm1_4)

        # Firm 2
        f2 = LowerOrEqualThanEquilibriumConstraint([g2[t], p[t], lambda2[t], alpha2[t], beta2[t]], [-2, 1, -1, 1, 1], 0)
        f2dual = GreaterOrEqualThanEquilibriumConstraint([g2[t]], [1], 0)
        firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
        add_equilibrium_constraint(m, firm2)

        f2_2 = LowerOrEqualThanEquilibriumConstraint([g2[t]], [1], -10)
        f2dual_2 = GreaterOrEqualThanEquilibriumConstraint([lambda2[t]], [1], 0)
        firm2_2 = ComplementarityEquilibriumConstraint(f2_2, f2dual_2)
        add_equilibrium_constraint(m, firm2_2)
        
        # # subida
        f2_3 = LowerOrEqualThanEquilibriumConstraint([g2[t], g2[t-1]], [1, -1], -8)
        f2dual_3 = GreaterOrEqualThanEquilibriumConstraint([alpha2[t]], [1], 0)
        firm2_3 = ComplementarityEquilibriumConstraint(f2_3, f2dual_3)
        add_equilibrium_constraint(m, firm2_3)
        # descida
        f2_4 = GreaterOrEqualThanEquilibriumConstraint([g2[t], g2[t-1]], [1, -1], 8)
        f2dual_4 = GreaterOrEqualThanEquilibriumConstraint([beta2[t]], [1], 0)
        firm2_4 = ComplementarityEquilibriumConstraint(f2_4, f2dual_4)
        add_equilibrium_constraint(m, firm2_4)
    end
    @variable(m, s[1:T+1] >=0)
    @constraint(m, [t=1:T], s[t] >= d[t] - g1[t] - g2[t])
    @constraint(m, [t=1:T], s[t] >= g1[t] + g2[t] - d[t] )
    @objective(m, Min, sum(s))
    println(m)
    status = solve(m)
    println(status)
    println("s=" , getvalue(s))
    println("d=" , getvalue(d))
    println("g1=", getvalue(g1))
    println("g2=", getvalue(g2))
    println("p=", getvalue(p))
    println("lambda1=", getvalue(lambda1))
    println("lambda2=", getvalue(lambda2))
    println("alpha1=", getvalue(alpha1))
    println("alpha2=", getvalue(alpha2))
    println("beta1=", getvalue(beta1))
    println("beta2=", getvalue(beta2))
    println("fobj=", getobjectivevalue(m))
    fob_g1 = -getvalue(g1[2])^2 + getvalue(g1[2]) * getvalue(p[2]) + 3
    fob_g2 = -getvalue(g2[2])^2 + getvalue(g2[2]) * getvalue(p[2]) + 5  
    fob_d  = (20 - getvalue(p[2]) ) * getvalue(d[2]) - 5/3 * getvalue(d[2])^2 
    println("fobj=" , fob_g1 + fob_g2 + fob_d)
    # writeLP(m, "lp", genericnames=false)
end

function example_4()

    m = Model(solver = CbcSolver())
    eq_m = marketEquilibriumModel()
 
    @variable(m, g1>=0)
    @variable(m, g2>=0)
    @variable(m, q1_l >=0)
    @variable(m, q2_l >=0)
    @variable(m, qc_l >=0)
    @variable(m, q1_e >=0)
    @variable(m, q2_e >=0)
    @variable(m, qc_e >=0)
    @variable(m, pe)
    @variable(m, pl)
    @variable(m, beta1>=0)
    @variable(m, beta2>=0)
    @variable(m, gamma1>=0)
    @variable(m, gamma2>=0)
    @variable(m, alpha>=0)
    @variable(m, pi)
    @variable(m, I1>=0)
    @variable(m, I2>=0)

    d = 100
    c1 = 200
    c2 = 300
    k1 = 6
    k2 = 5
    
    #g1
    c = LowerOrEqualThanEquilibriumConstraint([beta1, gamma1], [1, 1], -k1)
    cdual = GreaterOrEqualThanEquilibriumConstraint([I1], [1], 0)
    consumer = ComplementarityEquilibriumConstraint(c, cdual)
    add_equilibrium_constraint(m, consumer)
    add_equilibrium_constraint(eq_m, consumer)
    
    # g2
    c = LowerOrEqualThanEquilibriumConstraint([beta2, gamma2], [1, 1], -k2)
    cdual = GreaterOrEqualThanEquilibriumConstraint([I2], [1], 0)
    consumer = ComplementarityEquilibriumConstraint(c, cdual)
    add_equilibrium_constraint(m, consumer)
    add_equilibrium_constraint(eq_m, consumer)
    
    # g1 2
    f1 = LowerOrEqualThanEquilibriumConstraint([beta1, pi], [-1, +1], -c1)
    f1dual = GreaterOrEqualThanEquilibriumConstraint([g1], [1], 0)
    firm1 = ComplementarityEquilibriumConstraint(f1, f1dual)
    add_equilibrium_constraint(m, firm1)
    add_equilibrium_constraint(eq_m, firm1)
    # g2 2
    f1 = LowerOrEqualThanEquilibriumConstraint([beta2, pi], [-1, +1], -c2)
    f1dual = GreaterOrEqualThanEquilibriumConstraint([g2], [1], 0)
    firm1 = ComplementarityEquilibriumConstraint(f1, f1dual)
    add_equilibrium_constraint(m, firm1)
    add_equilibrium_constraint(eq_m, firm1)
    
    # g1 3
    f2 = LowerOrEqualThanEquilibriumConstraint([pl, gamma1], [1, -1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([q1_l], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m, firm2)
    add_equilibrium_constraint(eq_m, firm2)
    # g2 3
    f2 = LowerOrEqualThanEquilibriumConstraint([pl, gamma2], [1, -1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([q2_l], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m, firm2)
    add_equilibrium_constraint(eq_m, firm2)
    # g1 4 - energia
    f2 = LowerOrEqualThanEquilibriumConstraint([pe], [1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([q1_e], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m, firm2)
    add_equilibrium_constraint(eq_m, firm2)
    # g2 4 - energia
    f2 = LowerOrEqualThanEquilibriumConstraint([pe], [1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([q2_e], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m, firm2)
    add_equilibrium_constraint(eq_m, firm2)
    # consumidor
    f2 = LowerOrEqualThanEquilibriumConstraint([pe], [-1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([qc_e], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m, firm2)
    add_equilibrium_constraint(eq_m, firm2)
    # consumidor 2 
    f2 = LowerOrEqualThanEquilibriumConstraint([pl, alpha], [-1, 1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([qc_l], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m, firm2)
    add_equilibrium_constraint(eq_m, firm2)
    
    # Market Clearing 1
    cl = EqualEquilibriumConstraint([g1, g2], [1, 1], d)
    cldual = FreeEquilibriumConstraint([pi], [1], 0)
    clearing = ComplementarityEquilibriumConstraint(cl, cldual)
    add_equilibrium_constraint(m, clearing)
    add_equilibrium_constraint(eq_m, clearing)
    # Market Clearing 2
    cl = EqualEquilibriumConstraint([q1_e, q2_e, qc_e], [1, 1, -1], 0)
    cldual = FreeEquilibriumConstraint([pe], [1], 0)
    clearing = ComplementarityEquilibriumConstraint(cl, cldual)
    add_equilibrium_constraint(m, clearing)
    add_equilibrium_constraint(eq_m, clearing)
    # Market Clearing 3
    cl = EqualEquilibriumConstraint([q1_l, q2_l, qc_l], [1, 1, -1], 0)
    cldual = FreeEquilibriumConstraint([pl], [1], 0)
    clearing = ComplementarityEquilibriumConstraint(cl, cldual)
    add_equilibrium_constraint(m, clearing)
    add_equilibrium_constraint(eq_m, clearing)
    
    @objective(m, Min, g1+g2-d)
    println(m)
    status = solve(m)
    println(status)
    # println("d=" , getvalue(d))
    println("q1_l=" ,getvalue(q1_l))
    println("q2_l=" ,getvalue(q2_l))
    println("qc_l=" ,getvalue(qc_l))
    println("qc_l=" ,getvalue(qc_l))
    println("g1=",getvalue(g1))
    println("g2=",getvalue(g2))
    println("q1_l=",getvalue(q1_l))
    println("q2_l=",getvalue(q2_l))
    println("qc_l=",getvalue(qc_l))
    println("q1_e=",getvalue(q1_e))
    println("q2_e=",getvalue(q2_e))
    println("qc_e=",getvalue(qc_e))
    println("pe=",getvalue(pe))
    println("pl=",getvalue(pl))
    println("beta1=",getvalue(beta1))
    println("beta2=",getvalue(beta2))
    println("gamma1=",getvalue(gamma1))
    println("alpha=",getvalue(alpha))
    println("pi=",getvalue(pi))
    println("I1=",getvalue(I1))
    println("I2=",getvalue(I2))
    println("gamma2=",getvalue(gamma2))
    writeLP(m, "lp", genericnames=false)

    println("aaaaaaaaaaaaaaq")
    print(eq_m)
end

# run
# example_1()
# example_2()
# example_3()
example_4()