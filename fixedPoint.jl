
# build recursive function for firm 1, returning q1 as function of q1 old value and q2
function b1(q1::Float64, q2::Float64)
    return (100 - q1 - q2 ) / 2
end
# build recursive function for firm2
function b2(q1::Float64, q2::Float64)
    return 100 - q1 - q2 
end

# example 1
# ---------
function fixed_point_1(q1::Float64=0.0, q2::Float64=0.0;  maxit::Int = 40, ϵ::Float64=1e-6)
    p = 0
    for i in 1:maxit
        q1_new = b1(q1,q2)
        q2_new = b2(q1_new,q2)
        
        # stop criteria
        if abs(q2 - q2_new) < ϵ && abs(q1 - q1_new) < ϵ || i == maxit
            q1 = q1_new
            q2 = q2_new
            p = 100 - q1_new - q2_new
            println("number of iterations =",i)
            println("q1* =", round(q1_new, 4))
            println("q2* =", round(q2_new, 4))
            println("p* =", round(p, 4))
            break
        end
        # update
        q1 = q1_new
        q2 = q2_new
    end

    return q1, q2, p
end

fixed_point_1()

#
# example 2
# ---------
function b1_2(p::Float64, pco2::Float64)
    q1 = p / 2 - pco2 / 4
   return  q1
end
# build recursive function for firm2
function b2_2(q1::Float64)
    e1 = q1 / 2
    q2 = 30 - e1
    return q2
end

function fixed_point_2(q1::Float64=0.0, q2::Float64=0.0;  maxit::Int = 40, ϵ::Float64=1e-6)
    p = 0.0
    pco2 = 0.0
    for i in 1:maxit
        q2_new = b2_2(q1)
        p = 100 - q1 - q2_new
        pco2 = p - q2_new 
        q1_new = b1_2(p, pco2)
        # stop criteria
        if  i == maxit || abs(q2 - q2_new) < ϵ || abs(q1 - q1_new) < ϵ
            println("number of iterations =",i)
            println("q1* =", round(q1_new, 4))
            println("q2* =", round(q2_new, 4))
            println("p* =", round(p, 4))
            println("pco2* =", round(pco2, 4))
            q1 = q1_new
            q2 = q2_new
            break
        end
        # update
        q1 = q1_new
        q2 = q2_new
    end
    return q1, q2, p, pco2
end


fixed_point_2()