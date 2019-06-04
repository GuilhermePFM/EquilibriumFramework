using JuMP, Ipopt
using Plots
# plotlyjs()

M = 10^6

# 
# Modelo dos agentes
# ------------------

function firma1(p)
    m = Model(solver = IpoptSolver())
    @variable(m, q2 >= 0)
    @objective(m, Max, p*q2 - q2^2/2 )
    solve(m)
    # print(getvalue(q2))
    return getvalue(q2)
end

function firma1_rampa(p)
    T = length(p)
    m = Model(solver = IpoptSolver())
    @variable(m, q2[1:T] >= 0)
    @objective(m, Max, sum(p[t]*q2[t] - q2[t]^2/2 for t in 1:T) )
    @constraint(m, q2[1] - q2[2] <= 10)
    @constraint(m, q2[2] - q2[1] >= -10)
    solve(m)
    # print(getvalue(q2))
    return getvalue(q2)
end

function firma2(p)
    m = Model(solver = IpoptSolver())
    @variable(m, q1 >= 0)
    @objective(m, Max, p*q1 - q1^2 )
    solve(m)
    # print(getvalue(q1))
    return getvalue(q1)
end

function firma1_e_2_c02(p, pc02)
    m = Model(solver = IpoptSolver())
    @variable(m, q2 >= 0)
    @variable(m, q1 >= 0)
    @variable(m, e1 >= 0)
    @variable(m, e2 >= 0)
    
    @constraint(m, q1/2 <= e1)
    @constraint(m, q2 <= e2)
    @constraint(m, e1 + e2 <= 30)
    @objective(m, Max, p*q2 - q2^2/2 - pc02*e2 +  p*q1 - q1^2 - pc02*e1)
    solve(m)
    # print(getvalue(q2))
    return getvalue(q1), getvalue(q2) 
end

function firma2_rampa(p)
    T = length(p)
    m = Model(solver = IpoptSolver())
    @variable(m, q1[1:T] >= 0)
    @constraint(m, q1[1] - q1[2] <= 10)
    @constraint(m, q1[2] - q1[1] >= -10)
    @objective(m, Max, sum(p[t]*q1[t] - q1[t]^2 for t in 1:T))
    solve(m)
    # print(getvalue(q1))
    return getvalue(q1)
end

function consumidor(p)
    m = Model(solver = IpoptSolver())
    @variable(m, d >= 0)
    @objective(m, Max, (100 - d/2) * d - p * d )
    solve(m)
    # print(getvalue(d))
    return getvalue(d)
end

function consumidor_rampa(p)
    T = length(p)
    m = Model(solver = IpoptSolver())
    @variable(m, d[1:T] >= 0)
    @constraint(m, d[1] - d[2] <= 10)
    @constraint(m, d[2] - d[1] >= -10)
    @objective(m, Max, sum((100 - d[t]/2) * d[t] - p[t] * d[t] for t in 1:T))
    solve(m)
    # print(getvalue(d))
    return getvalue(d)
end

function eq_completo(p)
    m = Model(solver = IpoptSolver())
    @variable(m, d >= 0)
    @variable(m, q1 >= 0)
    @variable(m, q2 >= 0)
    @objective(m, Max, p*q2 - q2^2/2 + p*q1 - q1^2 + (100 - d/2) * d - p * d)
    solve(m)
    return getvalue(d), getvalue(q1), getvalue(q2)
end

#
# Metodos de solução
# ------------------

# iterativo basico
function iterativo()
    d_vec = []
    p_vec = []
    q1_vec = []
    q2_vec = []

    p = 0
    maxit=40
    for i in 1:maxit
        p+=5
        println(p)
        push!(p_vec, p)

        d = consumidor(p)
        push!(d_vec, d)
        q1 = firma1(p)
        push!(q1_vec, q1)
        q2 = firma2(p)
        push!(q2_vec, q2)
        # criterio de parada
        if false#d - q1 - q2 <= 1e-3
            println("q1*=",q1)
            println("q2*=",q2)
            println("p*=",p)
            println("d*=",q1+q2)
            break
        end
    end
    plot(p_vec, hcat(d_vec, q1_vec .+ q2_vec), title="Busca iterativa pelo equilibrio",label=["Oferta" "Demanda"],
    xlabel="Preço", ylabel="Quantidade")
end

# iterativo com bissecao
function solve_it!(p, p_vec, d_vec, q1_vec, q2_vec)
    d, q1, q2 = eq_completo(p)
    push!(p_vec, p)
    push!(d_vec, d)
    push!(q1_vec, q1)
    push!(q2_vec, q2)
    return d, q1 + q2
end

function iterativo_bissecao()
    d_vec = []
    p_vec = []
    q1_vec = []
    q2_vec = []

    pold = 100
    pnew = 100
    maxit= 40
    for i in 1:maxit
        # solve pmax
        dmax, qtmax = solve_it!(pnew, p_vec, d_vec, q1_vec, q2_vec)
        
        # criterio de parada
        if abs(dmax - qtmax) <= 1e-1
            println("q1*=",q1_vec[end])
            println("q2*=",q2_vec[end])
            println("p*=",pnew)
            println("d*=",qtmax)
            break
        end

        # faz a bissecao
        if dmax - qtmax <= 0
            # busco um novo preco na metade e mantenho o atual (mais tight)
            pold = pnew
            pnew = pnew/ 2
            @show pold
            @show pnew
        else
            @show pold
            @show pnew
            if pnew == pold
                # tenta ir mais pra direita
                pnew = pnew + 20
                pold = pold + 20
            else
                # mantenho o old (da direita) e tento buscar um mais pra direita do new
                pnew = ( pold + pnew) / 2
            end
        end        
    end
    plot(p_vec, hcat(d_vec, q1_vec .+ q2_vec), title="Busca iterativa pelo equilibrio",label=["Demanda" "Oferta"],
    xlabel="Preço", ylabel="Quantidade")
end

# iterativo ccom rampa
function solve_rampa!(p, d_vec, q1_vec, q2_vec)
    d = consumidor_rampa(p)
    q1 = firma1_rampa(p)
    q2 = firma2_rampa(p)
    push!(d_vec, d)
    push!(q1_vec, q1)
    push!(q2_vec, q2)
    return d, q1 + q2
end

function iterativo_rampa()
    d_vec = []
    p1_vec = []
    p2_vec = []
    q1_vec = []
    q2_vec = []
    
    p1 = 0
    p2 = 0
    maxit=40
    # loop para p1
    for i in 1:maxit
        p1+=5
        push!(p1_vec, p1)
        P2=0
        # loop para p2
        for i in 1:maxit
            p2+=5
            push!(p2_vec, p1)
            d, qt = solve_rampa!([p1,p2], d_vec, q1_vec, q2_vec)
            # criterio de parada interior
            if abs(d[1] - qt[1]) <= 1e-3
                println("interno")
                println("q1*=",q1_vec[end])
                println("q2*=",q2_vec[end])
                println("p2*=",p2)
                println("d*=",qt)
                println("------------")
                break
            end
        end
        # criterio de parada total
        if abs(d_vec[end][1] - q1_vec[end][1] - q2_vec[end][1]) <= 1e-3
            println("externo")
            println("q1*=",q1_vec[end])
            println("q2*=",q2_vec[end])
            println("p1*=",p1)
            println("p2*=",p2)
            println("d*=",q1_vec[end] .+ q2_vec[end])
            println("-------------")
            break
        end
    end
    # plot(p1_vec, hcat(d_vec, q1_vec .+ q2_vec), title="Busca iterativa pelo equilibrio",label=["Preço" "Quantidade"])
end

# iterativo com preco de emissao
function solve_it_c02(p, pc02, d_vec, q1_vec, q2_vec)
    d = consumidor(p)
    q1, q2 = firma1_e_2_c02(p, pc02)
    push!(d_vec, d)
    push!(q1_vec, q1)
    push!(q2_vec, q2)
    return d, q1 + q2
end

function iterativo_co2()
    
    d_vec = []
    p_vec = []
    pc02_vec = []
    q1_vec = []
    q2_vec = []
    
    p = 50
    pc02 = 30
    maxit=10
    
    # global iteration
    for i in 1:maxit
        # convergencia do p
        for i in 1:maxit
            # convergencia do pc02
            pc02 +=5
            push!(pc02_vec, pc02)
            p=0
            for i in 1:maxit
                p+=5
                println(p)
                push!(p_vec, p)
                d, qt = solve_it_c02(p, pc02, d_vec, q1_vec, q2_vec)
                
                # criterio de parada
            end
            if any(d_vec .- q1_vec .- q2_vec .<= 1e-3) 
                println("q1*=",q1_vec[end])
                println("q2*=",q2_vec[end])
                println("p*=",p)
                println("pc02*=",pc02)
                println("d*=",qt)
                break
            end

        end

    end
    plot(p_vec, hcat(d_vec, q1_vec .+ q2_vec), title="Busca iterativa pelo equilibrio",label=["Preço" "Quantidade"])
end 


# iterativo()
# iterativo_bissecao()
# iterativo_co2()
# iterativo_rampa()

# Modelo nao linear de rampa
# --------------------------
function f()
    T = 1
    g1_0 = 0
    g2_0 = 0
    d0 = 0
    ramp_g1 = 1
    ramp_g2 = 8

    m = Model(solver = IpoptSolver( acceptable_tol = 1e-3))
    # variables
    @variable(m, 0<=g1[1:T+1]<=3)
    @variable(m, 0<=g2[1:T+1]<=10)
    @variable(m, d[1:T+1]>=0)
    @variable(m, p[1:T+1])

    # initial condition
    @constraint(m, g1[1] == g1_0)
    @constraint(m, g2[1] == g2_0)
    @constraint(m, d[1] == d0)

    @constraint(m, demand[t=1:T+1],  d[t] == g1[t] + g2[t] )

    # ramp constraints
    @constraint(m, ramp1_down[t=1:T], g1[t+1] - g1[t] <= ramp_g1)
    @constraint(m, ramp1_up[t=2:T+1], g1[t-1] - g1[t] <= ramp_g1)

    @constraint(m, ramp2_down[t=1:T], g2[t+1] - g2[t] <= ramp_g2)
    @constraint(m, ramp2_up[t=2:T+1], g2[t-1] - g2[t] <= ramp_g2)

    fob_g1 = @expression(m, sum(-g1[t]^2 +  g1[t] * p[t] + 3 for t in 2:T+1))
    fob_g2 = @expression(m, sum(-g2[t]^2 + g2[t] * p[t] + 5   for t in 2:T+1))
    fob_d  = @expression(m, sum((20 - p[t] ) * d[t] - 5/3 * d[t]^2 for t in 2:T+1))
    @objective(m, Max, fob_g1 + fob_g2 + fob_d)
    solve(m)

    @show getvalue(g1)
    @show getvalue(g2)
    @show getvalue(d)
    @show getvalue(d) - getvalue(g1) - getvalue(g2)
    @show getvalue(p)
    @show getobjectivevalue(m)
end
f()