#--- Julia version 1.1

using JuMP, Ipopt



d = 100
c1 = 200
c2 = 300
k1 = 600
k2 = 601

function test_kkt( k1 , k2 , d , c1 , c2 )

    erro = 1e-3
    m = Model()

    #--- Variaveis primais
    #@variable( m , g1   >= 0   )
    @variable( m , g1   >= 0   )
    @variable( m , I1   >= 0   )
    @variable( m , q1_l >= 0   )
    @variable( m , q1_e >= 0   )

    @variable( m , qc_l >= 0   )
    @variable( m , qc_e >= 0   )

    #--- Variaveis duais
    @variable( m , g_beta   >= 0 )
    @variable( m , g_gamma  >= 0 )
    @variable( m , g_teta   >= 0 )
    @variable( m , g_qci    >= 0 )
    @variable( m , g_delta  >= 0 )
    @variable( m , g_epslon >= 0 )

    @variable( m , d_alpha  >= 0 )
    @variable( m , d_fi     >= 0 )
    @variable( m , d_ro     >= 0 )
    
    @variable( m , pi          )
    @variable( m , pe          )
    @variable( m , pl          )

    #------------
    #--- KKTs ---
    #------------

    #- Condicoes do primal

    @constraint( m , g1 - I1 <= 0 )

    @constraint( m , q1_l - I1 <= 0 )

    @constraint( m , d - qc_l <= 0 )

    @constraint( m , q1_e - qc_e == 0 )
    @constraint( m , q1_l -  qc_l == 0 )
    @constraint( m , g1 - d == 0 )

    #- Derivadas do Lagrangeano

    @constraint( m , -k1 + g_beta + g_gamma + g_qci == 0 )

    @constraint( m , -c1 - g_beta + pi + g_teta == 0 )

    @constraint( m , -g_gamma - pl + g_delta == 0 )
    
    @constraint( m , -pe + g_epslon == 0 )
    
    @constraint( m , d_alpha + pl + d_fi == 0 )

    @constraint( m , pe + d_ro == 0 )

    #- Produto das restricoes pela duais

    #@constraint( m , ( g1 - I1 ) * g_beta == 0 )
    @constraint( m , ( I1 - g1 ) * g_beta == 0 )
    
    #@constraint( m , ( q1_l - I1 ) * g_gamma == 0 )
    @constraint( m , ( I1 - q1_l ) * g_gamma == 0 )

    #@constraint( m , ( d - qc_l ) * d_alpha == 0 )
    @constraint( m , ( qc_l - d ) * d_alpha == 0 )

    #@constraint( m , ( q1_e - qc_e ) * pe == 0 )
    #@constraint( m , ( qc_e - q1_e ) * pe == 0 )

    #@constraint( m , ( q1_l - qc_l ) * pl == 0 )
    #@constraint( m , ( qc_l - q1_l ) * pl == 0 )

    #@constraint( m , ( g1 - d ) * pi == 0 )
    #@constraint( m , ( d - g1 ) * pi == 0 )

    @constraint( m , -g1 * g_teta == 0 )

    @constraint( m , -I1 * g_qci >= 0 )

    @constraint( m , -q1_l * g_delta == 0 )

    @constraint( m , -q1_e * g_epslon == 0 )

    @constraint( m , -qc_l * d_fi == 0 )

    @constraint( m , -qc_e * d_ro >= 0 )
    
    #-----------------------
    #--- Funcao Objetivo ---
    #-----------------------

    @objective( m , Max , qc_l  )

    
    set_optimizer( m , with_optimizer( Ipopt.Optimizer ) )
    MOIU.attach_optimizer( m )

    optimize!( m )


    println("#----------------------------------------------------#")
    println("#-----               Investimento               -----#")
    println("#----------------------------------------------------#\n")

    println( "Gerador 1 = " , round( value.( I1 ) ; digits = 2 ) , " MW")
    println("")

    println("#----------------------------------------------------#")
    println("#-----               Mercado Spot               -----#")
    println("#----------------------------------------------------#\n")
    
    println( "Gerador 1  = " , round( value.( g1 ) ; digits = 2 ) , " MWh")
    println( "Preco spot = " , round( value.( pi ) ; digits = 2 ) , " R\$/MWh" )
    println("")

    println("#----------------------------------------------------#")
    println("#-----             Mercado Energia              -----#")
    println("#----------------------------------------------------#\n")

    println("Quantidade - Gerador 1  = " , round( value.( q1_e ) ; digits = 2 ) , " MWm" )
    println("Quantidade - Consumidor = " , round( value.( qc_e ) ; digits = 2 ) , " MWm" )
    println("Preco contrato          = " , round( value.( pe )   ; digits = 2 ) , " R\$/MWh" )
    println("")

    println("#----------------------------------------------------#")
    println("#-----              Mercado Lastro              -----#")
    println("#----------------------------------------------------#\n")

    println("Quantidade - Gerador 1  = " , round( value.( q1_l ) ; digits = 2 ) , " MWm" )
    println("Quantidade - Consumidor = " , round( value.( qc_l ) ; digits = 2 ) , " MWm" )
    println("Preco contrato          = " , round( value.( pl )   ; digits = 2 ) , " R\$/MW" )
    println("")
    
    println("#----------------------------------------------------#")
    println("#-----         Variaveis Duais de Apoio         -----#")
    println("#----------------------------------------------------#\n")

    println("Beta   = " , round( value.( g_beta  ) ; digits = 3 ) )
    println("Gamma  = " , round( value.( g_gamma ) ; digits = 3 ) )
    println("Theta  = " , round( value.( g_teta ) ; digits = 3 ) )
    println("Qci    = " , round( value.( g_qci ) ; digits = 3 ) )
    println("Delta  = " , round( value.( g_delta ) ; digits = 3 ) )
    println("Epslon = " , round( value.( g_epslon ) ; digits = 3 ) )

    
    
    
    println("Alpha   = " , round( value.( d_alpha  ) ; digits = 3 ) )
    println("Fi      = " , round( value.( d_fi ) ; digits = 3 ) )
    println("Ro      = " , round( value.( d_ro ) ; digits = 3 ) )
    
    #- Post processing

    # if termination_status( m ) == MOI.OPTIMAL 

    #     fda = calculate_fda( value.( R ) , p );
    #     var = get_VaR( fda , d_alpha )
        
    #     println( string( " Q*       = " , round( value.( Q ) ; digits = 4 ) ) )
    #     println( string( " E[ R ]   = " , round( objective_value( m ) ; digits = 4 ) ) )
    #     println( string( " VaR[ R ] = " , round( var ; digits = 4 ) ) )
    #     println( " FDA \n" )
        
    #     return fda
    # else
    #     println( " Optimal solution not found ")
    #     return nothing
    # end

    return nothing

end



function kkt_operacao( d , I1 , c1 )

    #--- Modelo
    m = Model()

    #--- Variaveis primais
    @variable( m , g1   >= 0   )
    
    #--- Variaveis duais
    @variable( m , g_beta   >= 0 )
    @variable( m , g_gamma  >= 0 )
    
    @variable( m , pi )

    #------------
    #--- KKTs ---
    #------------

    #- Viabilidade primal
    @constraint( m , g1 <= I1 )
    @constraint( m , d - g1 == 0 )
    
    #- Viabilidade dual
    @constraint( m , -1*c1 - g_beta - pi + g_gamma == 0 )

    #- Complementariedade
    @constraint( m , -g1 * g_gamma == 0 )
    @constraint( m , ( g1 - I1 ) * g_beta == 0 )
    
    #-----------------------
    #--- Funcao Objetivo ---
    #-----------------------

    @objective( m , Min , d - g1 )

    set_optimizer( m , with_optimizer( Ipopt.Optimizer ) )
    MOIU.attach_optimizer( m )

    optimize!( m )

    println("#-------------------------------------------------#")
    println("#-----               Resultado               -----#")
    println("#-------------------------------------------------#\n")

    println( "Gerador 1  = " , round( value.( g1 ) ; digits = 2 ) , " MWh")
    println( "Preco spot = " , round( value.( pi ) ; digits = 2 ) , " R\$/MWh" )
    println( "Beta       = " , round( value.( g_beta  ) ; digits = 3 ) )
    println( "Gamma      = " , round( value.( g_gamma ) ; digits = 3 ) )
    
    return nothing

end

function kkt_expansao( d , k1 , c1 )

    #--- Modelo
    m = Model()

    #--- Variaveis primais
    @variable( m , g1   >= 0   )
    @variable( m , I1   >= 0   )
    
    #--- Variaveis duais
    @variable( m , g_beta   >= 0 )
    @variable( m , g_gamma  >= 0 )
    @variable( m , g_alpha  >= 0 )
    
    @variable( m , pi )

    #------------
    #--- KKTs ---
    #------------

    #- Viabilidade primal
    @constraint( m , g1 - I1 <= 0 )
    @constraint( m , d - g1 == 0 )
    
    #- Viabilidade dual
    @constraint( m , -1*c1 - g_beta - pi + g_gamma == 0 )
    @constraint( m , -k1 + g_beta + g_alpha        == 0 )

    #- Complementariedade
    @constraint( m , g1 * g_gamma == 0 )
    @constraint( m , I1 * g_alpha == 0 )
    @constraint( m , ( g1 - I1 ) * g_beta == 0 )
    
    #-----------------------
    #--- Funcao Objetivo ---
    #-----------------------

    @objective( m , Min , d - g1 )

    set_optimizer( m , with_optimizer( Ipopt.Optimizer ) )
    MOIU.attach_optimizer( m )

    optimize!( m )

    println("#-------------------------------------------------#")
    println("#-----               Resultado               -----#")
    println("#-------------------------------------------------#\n")

    println( "Investimento G1  = " , round( value.( I1 ) ; digits = 2 ) , " MW")
    println( "Geracao G1       = " , round( value.( g1 ) ; digits = 2 ) , " MWh")
    println( "Preco spot       = " , round( value.( pi ) ; digits = 2 ) , " R\$/MWh" )
    println( "Beta             = " , round( value.( g_beta  ) ; digits = 3 ) )
    println( "Gamma            = " , round( value.( g_gamma ) ; digits = 3 ) )
    println( "Alpha            = " , round( value.( g_alpha ) ; digits = 3 ) )
    
    return nothing

end

function kkt_lastro( d , k1 , c1 )

    #--- Modelo
    m = Model()

    #--- Variaveis primais
    @variable( m , g1    >= 0 )
    @variable( m , I1    >= 0 )
    @variable( m , ql_g1 >= 0 )
    @variable( m , ql_c  >= 0 )
    
    #--- Variaveis duais
    @variable( m , g_beta   >= 0 )
    @variable( m , g_gamma  >= 0 )
    @variable( m , g_alpha  >= 0 )
    @variable( m , g_ro     >= 0 )
    @variable( m , g_qci    >= 0 )
    @variable( m , c_phi    >= 0 )
    @variable( m , c_theta  >= 0 )
    
    @variable( m , pi )
    @variable( m , p_lastro )

    #------------
    #--- KKTs ---
    #------------

    #--- Viabilidade primal

    #- Gerador
    @constraint( m , g1 - I1    <= 0 )      # Capacidade maxima de geracao
    @constraint( m , ql_g1 - I1 <= 0 )      # Maximo volume de lastro que pode ser vendido
    
    #- Consumidor
    @constraint( m , d - ql_c <= 0 )        # Volume minimo de lastro a ser comprado
    
    #- Balanco
    @constraint( m , g1 - d       == 0 )            # Mercado spot
    @constraint( m , ql_g1 - ql_c == 0 )            # Mercado de lastro
    
    #--- Viabilidade dual

    #- Gerador
    @constraint( m , -1*c1 - g_beta - pi + g_gamma == 0 )   # Lagrange de g
    @constraint( m , -k1 + g_beta + g_alpha + g_ro == 0 )   # Lagrange de I
    @constraint( m , -g_ro - p_lastro + g_qci      == 0 )   # Lagrande de ql_g1

    #- Consumidor
    @constraint( m , c_phi + p_lastro + c_theta == 0 )      # Lagrange de ql_c

    #- Complementariedade
    @constraint( m , g1 * g_gamma == 0 )
    @constraint( m , I1 * g_alpha == 0 )
    @constraint( m , ql_g1 * g_qci == 0 )
    @constraint( m , ql_c * c_theta == 0 )

    @constraint( m , ( g1 - I1 ) * g_beta == 0 )
    @constraint( m , ( ql_g1 - I1 ) * g_ro == 0 )
    @constraint( m , ( d - ql_c ) * c_phi == 0 )
    
    #-----------------------
    #--- Funcao Objetivo ---
    #-----------------------

    @objective( m , Min , g1 - d )

    set_optimizer( m , with_optimizer( Ipopt.Optimizer ) )
    MOIU.attach_optimizer( m )

    optimize!( m )

    println("#----------------------------------------------------#")
    println("#-----               Investimento               -----#")
    println("#----------------------------------------------------#\n")

    println( "Gerador 1 = " , round( value.( I1 ) ; digits = 2 ) , " MW")
    println("")

    println("#----------------------------------------------------#")
    println("#-----               Mercado Spot               -----#")
    println("#----------------------------------------------------#\n")
    
    println( "Gerador 1  = " , round( value.( g1 ) ; digits = 2 ) , " MWh")
    println( "Preco spot = " , round( value.( pi ) ; digits = 2 ) , " R\$/MWh" )
    println("")

    println("#----------------------------------------------------#")
    println("#-----              Mercado Lastro              -----#")
    println("#----------------------------------------------------#\n")

    println("Quantidade - Gerador 1  = " , round( value.( ql_g1 ) ; digits = 2 ) , " MWm" )
    println("Quantidade - Consumidor = " , round( value.( ql_c ) ; digits = 2 ) , " MWm" )
    println("Preco contrato          = " , round( value.( p_lastro ) ; digits = 2 ) , " R\$/MW" )
    println("")
    
    println("#----------------------------------------------------#")
    println("#-----         Variaveis Duais de Apoio         -----#")
    println("#----------------------------------------------------#\n")

    println("Beta   = " , round( value.( g_beta  ) ; digits = 3 ) )
    println("Gamma  = " , round( value.( g_gamma ) ; digits = 3 ) )
    println("Alpha  = " , round( value.( g_alpha ) ; digits = 3 ) )
    println("Ro     = " , round( value.( g_ro ) ; digits = 3 ) )
    println("Qci    = " , round( value.( g_qci ) ; digits = 3 ) )
    
    println("Phi    = " , round( value.( c_phi  ) ; digits = 3 ) )
    println("Theta  = " , round( value.( c_theta ) ; digits = 3 ) )
    
    return nothing

end

function kkt_completo( d , k1 , c1 )

    #--- Modelo
    m = Model()

    #--- Variaveis primais
    @variable( m , g1    >= 0 )
    @variable( m , I1    >= 0 )
    @variable( m , ql_g1 >= 0 )
    @variable( m , ql_c  >= 0 )
    @variable( m , qe_g1 >= 0 )
    @variable( m , qe_c  >= 0 )
    
    #--- Variaveis duais
    @variable( m , g_beta   >= 0 )
    @variable( m , g_gamma  >= 0 )
    @variable( m , g_alpha  >= 0 )
    @variable( m , g_ro     >= 0 )
    @variable( m , g_qci    >= 0 )
    @variable( m , g_mu     >= 0 )
    @variable( m , c_phi    >= 0 )
    @variable( m , c_theta  >= 0 )
    @variable( m , c_omega  >= 0 )
    
    @variable( m , pi )
    @variable( m , p_lastro  )
    @variable( m , p_energia )

    #------------
    #--- KKTs ---
    #------------

    #--- Viabilidade primal

    #- Gerador
    @constraint( m , g1 - I1    <= 0 )      # Capacidade maxima de geracao
    @constraint( m , ql_g1 - I1 <= 0 )      # Maximo volume de lastro que pode ser vendido
    
    #- Consumidor
    @constraint( m , d - ql_c <= 0 )        # Volume minimo de lastro a ser comprado
    
    #- Balanco
    @constraint( m , g1 - d       == 0 )            # Mercado spot
    @constraint( m , ql_g1 - ql_c == 0 )            # Mercado de lastro
    @constraint( m , qe_g1 - qe_c == 0 )            # Mercado de energia
    
    #--- Viabilidade dual

    #- Gerador
    @constraint( m , -1*c1 - g_beta - pi + g_gamma == 0 )   # Lagrange de g
    @constraint( m , -k1 + g_beta + g_alpha + g_ro == 0 )   # Lagrange de I
    @constraint( m , -g_ro - p_lastro + g_qci      == 0 )   # Lagrande de ql_g1
    @constraint( m , -p_energia + g_mu             == 0 )   # Lagrande de qe_g1

    #- Consumidor
    @constraint( m , c_phi + p_lastro + c_theta == 0 )      # Lagrange de ql_c
    @constraint( m , p_energia + c_omega        == 0 )      # Lagrange de qe_c

    #- Complementariedade
    @constraint( m , g1 * g_gamma    == 0  )
    @constraint( m , I1 * g_alpha    == 0  )
    @constraint( m , ql_g1 * g_qci   == 0  )
    @constraint( m , ql_c * c_theta  == 0  )
    @constraint( m , qe_g1 * g_mu    == 0  )
    @constraint( m , qe_c  * c_omega == 0  )

    @constraint( m , ( g1 - I1 ) * g_beta  == 0 )
    @constraint( m , ( ql_g1 - I1 ) * g_ro == 0 )
    @constraint( m , ( d - ql_c ) * c_phi  == 0 )
    
    #-----------------------
    #--- Funcao Objetivo ---
    #-----------------------
    @variable(m , a >=0)
    @constraint(m , a >=g1 - d)
    @constraint(m , a >= d- g1)
    # @constraint(m , g1 >= 100)
    @objective( m , Min , a )

    set_optimizer( m , with_optimizer( Ipopt.Optimizer ) )
    MOIU.attach_optimizer( m )

    optimize!( m )

    println("#----------------------------------------------------#")
    println("#-----               Investimento               -----#")
    println("#----------------------------------------------------#\n")

    println( "Gerador 1 = " , round( value.( I1 ) ; digits = 2 ) , " MW")
    println("")

    println("#----------------------------------------------------#")
    println("#-----               Mercado Spot               -----#")
    println("#----------------------------------------------------#\n")
    
    println( "Gerador 1  = " , round( value.( g1 ) ; digits = 2 ) , " MWh")
    println( "Preco spot = " , round( value.( pi ) ; digits = 2 ) , " R\$/MWh" )
    println("")

    println("#----------------------------------------------------#")
    println("#-----              Mercado Lastro              -----#")
    println("#----------------------------------------------------#\n")

    println("Quantidade - Gerador 1  = " , round( value.( ql_g1 ) ; digits = 2 ) , " MWm" )
    println("Quantidade - Consumidor = " , round( value.( ql_c ) ; digits = 2 ) , " MWm" )
    println("Preco contrato          = " , round( value.( p_lastro ) ; digits = 2 ) , " R\$/MW" )
    println("")

    println("#----------------------------------------------------#")
    println("#-----             Mercado Energia              -----#")
    println("#----------------------------------------------------#\n")

    println("Quantidade - Gerador 1  = " , round( value.( qe_g1 ) ; digits = 2 ) , " MWm" )
    println("Quantidade - Consumidor = " , round( value.( qe_c ) ; digits = 2 ) , " MWm" )
    println("Preco contrato          = " , round( value.( p_energia ) ; digits = 2 ) , " R\$/MWh" )
    println("")
    
    println("#----------------------------------------------------#")
    println("#-----         Variaveis Duais de Apoio         -----#")
    println("#----------------------------------------------------#\n")

    println("Beta   = " , round( value.( g_beta  ) ; digits = 3 ) )
    println("Gamma  = " , round( value.( g_gamma ) ; digits = 3 ) )
    println("Alpha  = " , round( value.( g_alpha ) ; digits = 3 ) )
    println("Ro     = " , round( value.( g_ro ) ; digits = 3 ) )
    println("Qci    = " , round( value.( g_qci ) ; digits = 3 ) )
    
    println("Phi    = " , round( value.( c_phi  ) ; digits = 3 ) )
    println("Theta  = " , round( value.( c_theta ) ; digits = 3 ) )
    
    return nothing

end
kkt_completo( 100, 300 , 150 )