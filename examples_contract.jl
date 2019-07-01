include("marketEquilibrium.jl")
using JuMP
using Cbc

#=========================================================================#
#=====                  Social Welfare - Approach 2                  =====#
#=========================================================================#

function SocialWelfare_1( path , d , c1 , c2 , k1 , k2 )

    m = Model(solver = ClpSolver())
    
    @variable( m , g1   >=0   )
    @variable( m , g2   >=0   )
    @variable( m , I1   >=0   )
    @variable( m , I2   >=0   )
    @variable( m , q1_l >=0   )
    @variable( m , q2_l >=0   )
    @variable( m , qc_l >=0   )
    @variable( m , q1_e >=0   )
    @variable( m , q2_e >=0   )
    @variable( m , qc_e >=0   )
    @variable( m , pi         )
    @variable( m , pe         )
    @variable( m , pl         )
    @variable( m , beta1  >=0 )
    @variable( m , beta2  >=0 )
    @variable( m , gamma1 >=0 )
    @variable( m , gamma2 >=0 )
    @variable( m , alpha      )

    @constraintref max_gen_1   
    @constraintref max_gen_2   
    @constraintref max_contr_1 
    @constraintref max_contr_2 
    @constraintref min_contr 
    @constraintref balanco_spot 
    @constraintref balanco_lastro 
    @constraintref balanco_energia 

    max_gen_1 = @constraint( m , g1 <= I1 )
    max_gen_2 = @constraint( m , g2 <= I2 )

    max_contr_1 = @constraint( m , q1_l <= I1 )
    max_contr_2 = @constraint( m , q2_l <= I2 )

    min_contr = @constraint( m , d <= qc_l )
    
    balanco_spot    = @constraint( m ,  d == g1 + g2 )
    balanco_energia = @constraint( m ,  qc_e == q1_e + q2_e)
    balanco_lastro  = @constraint( m ,  qc_l == q1_l + q2_l )

    @objective( m , Max , -I1 * k1 -I2 * k2 - g1*c1 - g2*c2 )

    status = solve(m)
    #println(status)
    # println("d=" , getvalue(d))

    beta1 = getdual( max_gen_1 )
    beta2 = getdual( max_gen_2 )
    gamma1 = getdual( max_contr_1 )
    gamma2 = getdual( max_contr_2 )
    alpha  = getdual( min_contr )
    pi = getdual( balanco_spot )
    pl = getdual( balanco_lastro )
    pe = getdual( balanco_energia )

    println("#----------------------------------------------------#")
    println("#-----               Investimento               -----#")
    println("#----------------------------------------------------#\n")

    println( "Gerador 1 = " , round( getvalue( I1 ) , 2 ) , " MW")
    println( "Gerador 2 = " , round( getvalue( I2 ) , 2 ) , " MW")
    println("")

    println("#----------------------------------------------------#")
    println("#-----               Mercado Spot               -----#")
    println("#----------------------------------------------------#\n")
    
    println( "Gerador 1  = " , round( getvalue( g1 ) , 2 ) , " MWh")
    println( "Gerador 2  = " , round( getvalue( g2 ) , 2 ) , " MWh")
    println( "Preco spot = " , round( pi , 2 ) , " R\$/MWh" )
    println("")

    println("#----------------------------------------------------#")
    println("#-----             Mercado Energia              -----#")
    println("#----------------------------------------------------#\n")

    println("Quantidade - Gerador 1  = " , round( getvalue( q1_e ) , 2 ) , " MWm" )
    println("Quantidade - Gerador 2  = " , round( getvalue( q2_e ) , 2 ) , " MWm" )
    println("Quantidade - Consumidor = " , round( getvalue( qc_e ) , 2 ) , " MWm" )
    println("Preco contrato          = " , round( pe   , 2 ) , " R\$/MWh" )
    println("")

    println("#----------------------------------------------------#")
    println("#-----              Mercado Lastro              -----#")
    println("#----------------------------------------------------#\n")

    println("Quantidade - Gerador 1  = " , round( getvalue( q1_l ) , 2 ) , " MWm" )
    println("Quantidade - Gerador 2  = " , round( getvalue( q2_l ) , 2 ) , " MWm" )
    println("Quantidade - Consumidor = " , round( getvalue( qc_l ) , 2 ) , " MWm" )
    println("Preco contrato          = " , round( pl   , 2 ) , " R\$/MW" )
    println("")
    
    println("#----------------------------------------------------#")
    println("#-----         Variaveis Duais de Apoio         -----#")
    println("#----------------------------------------------------#\n")

    println("Beta 1  = " , round( beta1 , 3 ) )
    println("Beta 2  = " , round( beta2 , 3 ) )
    println("Gamma 1 = " , round( gamma1 , 3 ) )
    println("Gamma 2 = " , round( gamma2 , 3 ) )
    println("Alpha   = " , round( alpha , 3 ) )
    
    writeLP( m, joinpath( path , "debug.lp" ) , genericnames=false )

    return nothing
end

#=========================================================================#
#=====                  Social Welfare - Approach 1                  =====#
#=========================================================================#

function SocialWelfare_2( path , d , c1 , c2 , k1 , k2 )

    # m = Model(solver = CbcSolver())
    m = marketEquilibriumModel()
    # eq_m = marketEquilibriumModel()
    
    @variable( m.m , g1   >=0   )
    @variable( m.m , g2   >=0   )
    @variable( m.m , I1   >=0   )
    @variable( m.m , I2   >=0   )
    @variable( m.m , q1_l >=0   )
    @variable( m.m , q2_l >=0   )
    @variable( m.m , qc_l >=0   )
    @variable( m.m , q1_e >=0   )
    @variable( m.m , q2_e >=0   )
    @variable( m.m , qc_e >=0   )
    @variable( m.m , pi         )
    @variable( m.m , pe         )
    @variable( m.m , pl         )
    @variable( m.m , beta1  >=0 )
    @variable( m.m , beta2  >=0 )
    @variable( m.m , gamma1 >=0 )
    @variable( m.m , gamma2 >=0 )
    @variable( m.m , alpha      )

    #g1
    c = LowerOrEqualThanEquilibriumConstraint([beta1, gamma1], [1, 1], -k1)
    cdual = GreaterOrEqualThanEquilibriumConstraint([I1], [1], 0)
    consumer = ComplementarityEquilibriumConstraint(c, cdual)
    add_equilibrium_constraint(m.m, consumer)
    add_equilibrium_constraint(m, consumer)
    
    # g2
    c = LowerOrEqualThanEquilibriumConstraint([beta2, gamma2], [1, 1], -k2)
    cdual = GreaterOrEqualThanEquilibriumConstraint([I2], [1], 0)
    consumer = ComplementarityEquilibriumConstraint(c, cdual)
    add_equilibrium_constraint(m.m, consumer)
    add_equilibrium_constraint(m, consumer)
    
    # g1 2
    f1 = LowerOrEqualThanEquilibriumConstraint([beta1, pi], [-1, 1], -c1)
    f1dual = GreaterOrEqualThanEquilibriumConstraint([g1], [1], 0)
    firm1 = ComplementarityEquilibriumConstraint(f1, f1dual)
    add_equilibrium_constraint(m.m, firm1)
    add_equilibrium_constraint(m, firm1)
    # g2 2
    f1 = LowerOrEqualThanEquilibriumConstraint([beta2, pi], [-1, 1], -c2)
    f1dual = GreaterOrEqualThanEquilibriumConstraint([g2], [1], 0)
    firm1 = ComplementarityEquilibriumConstraint(f1, f1dual)
    add_equilibrium_constraint(m.m, firm1)
    add_equilibrium_constraint(m, firm1)
    
    # g1 3
    f2 = LowerOrEqualThanEquilibriumConstraint([pl, gamma1], [1, -1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([q1_l], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m.m, firm2)
    add_equilibrium_constraint(m, firm2)
    # g2 3
    f2 = LowerOrEqualThanEquilibriumConstraint([pl, gamma2], [1, -1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([q2_l], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m, firm2)

    # g1 4 - energia
    f2 = LowerOrEqualThanEquilibriumConstraint([pe], [1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([q1_e], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m.m, firm2)
    add_equilibrium_constraint(m, firm2)
    # g2 4 - energia
    f2 = LowerOrEqualThanEquilibriumConstraint([pe], [1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([q2_e], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m.m, firm2)
    add_equilibrium_constraint(m, firm2)
    # consumidor
    f2 = LowerOrEqualThanEquilibriumConstraint([pe], [-1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([qc_e], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m.m, firm2)
    add_equilibrium_constraint(m, firm2)

    #dual 1 - g1
    f2 = LowerOrEqualThanEquilibriumConstraint([I1 , g1], [-1,1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([beta1], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m.m, firm2)
    add_equilibrium_constraint(m, firm2)

    #dual 1 - g2
    f2 = LowerOrEqualThanEquilibriumConstraint([I2 , g2], [-1,1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([beta2], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m.m, firm2)
    add_equilibrium_constraint(m, firm2)
 

    #- Dual 2 - g1
    f2 = LowerOrEqualThanEquilibriumConstraint([I1 , q1_l], [-1,1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([gamma1], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m.m, firm2)
    add_equilibrium_constraint(m, firm2)

    #- Dual 2 - g2
    f2 = LowerOrEqualThanEquilibriumConstraint([I2 , q2_l], [-1,1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([gamma2], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m.m, firm2)
    add_equilibrium_constraint(m, firm2)

    
    # consumidor 2 
    f2 = LowerOrEqualThanEquilibriumConstraint([pl, alpha], [-1, 1], 0)
    f2dual = GreaterOrEqualThanEquilibriumConstraint([qc_l], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m.m, firm2)
    add_equilibrium_constraint(m, firm2)

    #-Consumidor - Dual 1
    f2 = LowerOrEqualThanEquilibriumConstraint([-qc_l], [1], d )
    f2dual = GreaterOrEqualThanEquilibriumConstraint([alpha], [1], 0)
    firm2 = ComplementarityEquilibriumConstraint(f2, f2dual)
    add_equilibrium_constraint(m.m, firm2)
    add_equilibrium_constraint(m, firm2)
    
    # Market Clearing 1
    cl = EqualEquilibriumConstraint([g1, g2], [1, 1], -d)
    cldual = FreeEquilibriumConstraint([pi], [1], 0)
    clearing = ComplementarityEquilibriumConstraint(cl, cldual)
    add_equilibrium_constraint(m.m, clearing)
    add_equilibrium_constraint(m, clearing)

    # Market Clearing 2
    cl = EqualEquilibriumConstraint([q1_e, q2_e, qc_e], [1, 1, -1], 0)
    cldual = FreeEquilibriumConstraint([pe], [1], 0)
    clearing = ComplementarityEquilibriumConstraint(cl, cldual)
    add_equilibrium_constraint(m.m, clearing)
    add_equilibrium_constraint(m, clearing)
    
    # Market Clearing 3
    cl = EqualEquilibriumConstraint([q1_l, q2_l, qc_l], [1, 1, -1], 0)
    cldual = FreeEquilibriumConstraint([pl], [1], 0)
    clearing = ComplementarityEquilibriumConstraint(cl, cldual)
    add_equilibrium_constraint(m.m, clearing)
    add_equilibrium_constraint(m, clearing)
    
    #@objective(m, Min, g1+g2-d)
    
    status = solve(m.m)
    #println(status)
    # println("d=" , getvalue(d))
    print(m)
    println(status)
    println("\n\n")

    println("#----------------------------------------------------#")
    println("#-----               Investimento               -----#")
    println("#----------------------------------------------------#\n")

    println( "Gerador 1 = " , round( getvalue( I1 ) , 2 ) , " MW")
    println( "Gerador 2 = " , round( getvalue( I2 ) , 2 ) , " MW")
    println("")

    println("#----------------------------------------------------#")
    println("#-----               Mercado Spot               -----#")
    println("#----------------------------------------------------#\n")
    
    println( "Gerador 1  = " , round( getvalue( g1 ) , 2 ) , " MWh")
    println( "Gerador 2  = " , round( getvalue( g2 ) , 2 ) , " MWh")
    println( "Preco spot = " , round( getvalue( pi ) , 2 ) , " R\$/MWh" )
    println("")

    println("#----------------------------------------------------#")
    println("#-----             Mercado Energia              -----#")
    println("#----------------------------------------------------#\n")

    println("Quantidade - Gerador 1  = " , round( getvalue( q1_e ) , 2 ) , " MWm" )
    println("Quantidade - Gerador 2  = " , round( getvalue( q2_e ) , 2 ) , " MWm" )
    println("Quantidade - Consumidor = " , round( getvalue( qc_e ) , 2 ) , " MWm" )
    println("Preco contrato          = " , round( getvalue( pe )   , 2 ) , " R\$/MWh" )
    println("")

    println("#----------------------------------------------------#")
    println("#-----              Mercado Lastro              -----#")
    println("#----------------------------------------------------#\n")

    println("Quantidade - Gerador 1  = " , round( getvalue( q1_l ) , 2 ) , " MWm" )
    println("Quantidade - Gerador 2  = " , round( getvalue( q2_l ) , 2 ) , " MWm" )
    println("Quantidade - Consumidor = " , round( getvalue( qc_l ) , 2 ) , " MWm" )
    println("Preco contrato          = " , round( getvalue( pl )   , 2 ) , " R\$/MW" )
    println("")
    
    println("#----------------------------------------------------#")
    println("#-----         Variaveis Duais de Apoio         -----#")
    println("#----------------------------------------------------#\n")

    println("Beta 1  = " , round( getvalue( beta1  ) , 3 ) )
    println("Beta 2  = " , round( getvalue( beta2  ) , 3 ) )
    println("Gamma 1 = " , round( getvalue( gamma1 ) , 3 ) )
    println("Gamma 2 = " , round( getvalue( gamma2 ) , 3 ) )
    println("Alpha   = " , round( getvalue( alpha  ) , 3 ) )
    
    # writeLP( m.m, joinpath( path , "debug.lp" ) , genericnames=false )

    println( "" )

    return 0
end

#=========================================================================#
#=====               Social Welfare - Approach 3 Complementarity package =#
#=========================================================================#

function SocialWelfare_3( path , d , c1 , c2 , k1 , k2 )

    m = Model(solver=IpoptSolver())
    
    @variable( m , g1   >= 0   )
    @variable( m , g2   >= 0   )
    @variable( m , I1   >= 0   )
    @variable( m , I2   >= 0   )
    @variable( m , q1_l >= 0   )
    @variable( m , q2_l >= 0   )
    @variable( m , qc_l >= 0   )
    @variable( m , q1_e >= 0   )
    @variable( m , q2_e >= 0   )
    @variable( m , qc_e >= 0   )
    @variable( m , pi          )
    @variable( m , pe          )
    @variable( m , pl          )
    @variable( m , beta1  >= 0 )
    @variable( m , beta2  >= 0 )
    @variable( m , gamma1 >= 0 )
    @variable( m , gamma2 >= 0 )
    @variable( m , alpha  >= 0 )

    #g1
    @complements(m, 0 <= beta1 + gamma1 -k1, I1 >= 0)
    
    # g2
    @complements(m, 0 <= beta2 + gamma2 -k2, I2 >= 0)
    
    # g1 2
    @complements(m, 0 <= -beta1 + pi -c1, g1 >= 0)
    
    # g2 2
    @complements(m, 0 <= -beta2+ pi -c2, g2 >= 0)
    
    # g1 3
    @complements(m, 0 <=pl - gamma1, q1_l >= 0)
    
    # g2 3
    @complements(m, 0 <= pl - gamma2, q2_l >= 0)
    
    # g1 4 - energia
    @complements(m, 0 <= pe, q1_e >= 0)
    
    # g2 4 - energia
    @complements(m, 0 <= pe, q2_e >= 0)
    
    #dual 1 - g1
    @complements(m, 0 <= I1 - g1, beta1 >= 0)
    
    #dual 1 - g2
    @complements(m, 0 <= I2 - g2, beta2 >= 0)
    
    #- Dual 2 - g1
    @complements(m, 0 <= I1 - q1_l, gamma1 >= 0)
    
    #- Dual 2 - g2
    @complements(m, 0 <= I2 - q2_l, gamma2 >= 0)
    
    # consumidor
    @complements(m, 0 <= pe, qc_e >= 0)
    
    # consumidor 2 
    @complements(m, 0 <= -pl+ alpha, qc_l >= 0)
    
    #-Consumidor - Dual 1
    @complements(m, 0 <= qc_l - d, alpha >= 0)
    
    # Market Clearing 1
    @complements(m, 0 <= g1 + g2 - d, pi >= 0)
    
    # Market Clearing 2
    @complements(m, 0 <= q1_e + q2_e - qc_e, pe >= 0)
    
    # Market Clearing 3
    @complements(m, 0 <= q1_l+ q2_l- qc_l, pl >= 0)
    
    #@objective(m, Min, g1+g2-d)
    
    status = solve(m)
    #println(status)
    # println("d=" , getvalue(d))

    println("#----------------------------------------------------#")
    println("#-----               Investimento               -----#")
    println("#----------------------------------------------------#\n")

    println( "Gerador 1 = " , round( getvalue( I1 ) , 2 ) , " MW")
    println( "Gerador 2 = " , round( getvalue( I2 ) , 2 ) , " MW")
    println("")

    println("#----------------------------------------------------#")
    println("#-----               Mercado Spot               -----#")
    println("#----------------------------------------------------#\n")
    
    println( "Gerador 1  = " , round( getvalue( g1 ) , 2 ) , " MWh")
    println( "Gerador 2  = " , round( getvalue( g2 ) , 2 ) , " MWh")
    println( "Preco spot = " , round( getvalue( pi ) , 2 ) , " R\$/MWh" )
    println("")

    println("#----------------------------------------------------#")
    println("#-----             Mercado Energia              -----#")
    println("#----------------------------------------------------#\n")

    println("Quantidade - Gerador 1  = " , round( getvalue( q1_e ) , 2 ) , " MWm" )
    println("Quantidade - Gerador 2  = " , round( getvalue( q2_e ) , 2 ) , " MWm" )
    println("Quantidade - Consumidor = " , round( getvalue( qc_e ) , 2 ) , " MWm" )
    println("Preco contrato          = " , round( getvalue( pe )   , 2 ) , " R\$/MWh" )
    println("")

    println("#----------------------------------------------------#")
    println("#-----              Mercado Lastro              -----#")
    println("#----------------------------------------------------#\n")

    println("Quantidade - Gerador 1  = " , round( getvalue( q1_l ) , 2 ) , " MWm" )
    println("Quantidade - Gerador 2  = " , round( getvalue( q2_l ) , 2 ) , " MWm" )
    println("Quantidade - Consumidor = " , round( getvalue( qc_l ) , 2 ) , " MWm" )
    println("Preco contrato          = " , round( getvalue( pl )   , 2 ) , " R\$/MW" )
    println("")
    
    println("#----------------------------------------------------#")
    println("#-----         Variaveis Duais de Apoio         -----#")
    println("#----------------------------------------------------#\n")

    println("Beta 1  = " , round( getvalue( beta1  ) , 3 ) )
    println("Beta 2  = " , round( getvalue( beta2  ) , 3 ) )
    println("Gamma 1 = " , round( getvalue( gamma1 ) , 3 ) )
    println("Gamma 2 = " , round( getvalue( gamma2 ) , 3 ) )
    println("Alpha   = " , round( getvalue( alpha  ) , 3 ) )
    
    writeLP( m, joinpath( path , "debug.lp" ) , genericnames=false )

    #println( "" )
    #println( m )

    return 0
end

#=========================================================================#
#=====               Social Welfare - Approach 4 Non linear multiplicat  =#
#=========================================================================#

function SocialWelfare_4( path , d , c1 , c2 , k1 , k2 )

    m = Model(solver=IpoptSolver())
    
    @variable( m , g1   >= 0   )
    @variable( m , g2   >= 0   )
    @variable( m , I1   >= 0   )
    @variable( m , I2   >= 0   )
    @variable( m , q1_l >= 0   )
    @variable( m , q2_l >= 0   )
    @variable( m , qc_l >= 0   )
    @variable( m , q1_e >= 0   )
    @variable( m , q2_e >= 0   )
    @variable( m , qc_e >= 0   )
    @variable( m , pi          )
    @variable( m , pe          )
    @variable( m , pl          )
    @variable( m , beta1  >= 0 )
    @variable( m , beta2  >= 0 )
    @variable( m , gamma1 >= 0 )
    @variable( m , gamma2 >= 0 )
    @variable( m , alpha  >= 0 )

    eps = 1e-3
    # Estacionariedade (derivadas)
    # investimento
    @constraint(m, ( beta1 + gamma1 -k1)*(I1)>=-eps)
    @constraint(m, ( beta1 + gamma1 -k1)*(I1)<=eps)
    
    @constraint(m,  (beta2+ gamma2 -k2)*(I2)>=-eps)
    @constraint(m,  (beta2+ gamma2 -k2)*(I2)<=eps)
    
    # Geração
    @constraint(m,  (-beta1+ pi -c1)*(g1)>=-eps)
    @constraint(m,  (-beta1+ pi -c1)*(g1)<=eps)

    @constraint(m,  (-beta2+ pi -c2)*(g2)>=-eps)
    @constraint(m,  (-beta2+ pi -c2)*(g2)<=eps)
    
    # Quantidade de Lastro
    @constraint(m,  (pl- gamma1)* q1_l>=-eps)
    @constraint(m,  (pl- gamma1)* q1_l<=eps)
   
    @constraint(m,  (pl - gamma2 )* q2_l>=-eps)
    @constraint(m,  (pl - gamma2 )* q2_l<=eps)
    
    # Quantidade de Energia
    @constraint(m,  (pe) * q1_e>=-eps)
    @constraint(m,  (pe) * q1_e<=eps)
    
    @constraint(m,  pe* q2_e>=-eps)
    @constraint(m,  pe* q2_e<=eps)

    # Quantidade de Energia - demanda
    @constraint(m,  (-pe)* qc_e>=-eps)
    @constraint(m,  (-pe)* qc_e<=eps)
    
    # Quantidade de Lastro - demanda
    @constraint(m,  (-pl + alpha)* qc_l>=-eps)
    @constraint(m,  (-pl + alpha)* qc_l<=eps)
    
    # Duais
    # -----
    # Restrição de demanda lastreada
    @constraint(m,  (qc_l - d)*(alpha) >=-eps)
    @constraint(m,  (qc_l - d)*(alpha) <=eps)

    # Geração limitada por investimento
    @constraint(m,  (-I1 + g1)* beta1>=-eps)
    @constraint(m,  (-I1 + g1)* beta1<=eps)
    @constraint(m,  (-I2 + g2)* beta2>=-eps)
    @constraint(m,  (-I2 + g2)* beta2<=eps)
    
    # Lastro limitado por investimento
    @constraint(m,  (I1 - q1_l)* gamma1>=-eps)
    @constraint(m,  (I1 - q1_l)* gamma1<=eps)
    @constraint(m,  (I2 - q2_l )* gamma2>=-eps)
    @constraint(m,  (I2 - q2_l )* gamma2<=eps)

    # Clearing
    # --------
    
    # Clearing demanda
    @constraint(m,  (g1 + g2 - d)*(pi)>=-eps)
    @constraint(m,  (g1 + g2 - d)*(pi)<=eps)
    
    # Clearing Energia
    @constraint(m,  (q1_e+ q2_e- qc_e)*(pe)>=-eps)
    @constraint(m,  (q1_e+ q2_e- qc_e)*(pe)<=eps)
    
    # Clearing Lastro
    @constraint(m,  (q1_l+ q2_l- qc_l)*(pl)>=-eps)
    @constraint(m,  (q1_l+ q2_l- qc_l)*(pl)<=eps)
    
    # @objective(m, Min, (g1+g2-d)^2)
    # @constraint(m, g1+g2-d)
    
    status = solve(m)
    #println(status)
    # println("d=" , getvalue(d))

    println("#----------------------------------------------------#")
    println("#-----               Investimento               -----#")
    println("#----------------------------------------------------#\n")

    println( "Gerador 1 = " , round( getvalue( I1 ) , 2 ) , " MW")
    println( "Gerador 2 = " , round( getvalue( I2 ) , 2 ) , " MW")
    println("")

    println("#----------------------------------------------------#")
    println("#-----               Mercado Spot               -----#")
    println("#----------------------------------------------------#\n")
    
    println( "Gerador 1  = " , round( getvalue( g1 ) , 2 ) , " MWh")
    println( "Gerador 2  = " , round( getvalue( g2 ) , 2 ) , " MWh")
    println( "Preco spot = " , round( getvalue( pi ) , 2 ) , " R\$/MWh" )
    println("")

    println("#----------------------------------------------------#")
    println("#-----             Mercado Energia              -----#")
    println("#----------------------------------------------------#\n")

    println("Quantidade - Gerador 1  = " , round( getvalue( q1_e ) , 2 ) , " MWm" )
    println("Quantidade - Gerador 2  = " , round( getvalue( q2_e ) , 2 ) , " MWm" )
    println("Quantidade - Consumidor = " , round( getvalue( qc_e ) , 2 ) , " MWm" )
    println("Preco contrato          = " , round( getvalue( pe )   , 2 ) , " R\$/MWh" )
    println("")

    println("#----------------------------------------------------#")
    println("#-----              Mercado Lastro              -----#")
    println("#----------------------------------------------------#\n")

    println("Quantidade - Gerador 1  = " , round( getvalue( q1_l ) , 2 ) , " MWm" )
    println("Quantidade - Gerador 2  = " , round( getvalue( q2_l ) , 2 ) , " MWm" )
    println("Quantidade - Consumidor = " , round( getvalue( qc_l ) , 2 ) , " MWm" )
    println("Preco contrato          = " , round( getvalue( pl )   , 2 ) , " R\$/MW" )
    println("")
    
    println("#----------------------------------------------------#")
    println("#-----         Variaveis Duais de Apoio         -----#")
    println("#----------------------------------------------------#\n")

    println("Beta 1  = " , round( getvalue( beta1  ) , 3 ) )
    println("Beta 2  = " , round( getvalue( beta2  ) , 3 ) )
    println("Gamma 1 = " , round( getvalue( gamma1 ) , 3 ) )
    println("Gamma 2 = " , round( getvalue( gamma2 ) , 3 ) )
    println("Alpha   = " , round( getvalue( alpha  ) , 3 ) )
    
    # writeLP( m, joinpath( path , "debug.lp" ) , genericnames=false )

    #println( "" )
    #println( m )

    return 0
end

#=========================================================================#
#=====                              Run                              =====#
#=========================================================================#

#--- Defining constants
path = "D:\\Dropbox (PSR)\\Mateus\\Diversos\\Lastro_Energia"
d = 100
c1 = 200
c2 = 300
k1 = 600
k2 = 601

#--- Calling functions

#SocialWelfare_1( path , d , c1 , c2 , k1 , k2 );
SocialWelfare_2( path , d , c1 , c2 , k1 , k2 );