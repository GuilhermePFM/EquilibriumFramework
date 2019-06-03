
# build recursive function for firm 1, returning q1 as function of q1 old value and q2
function b1(q1,q2)
    return (100 - q1 - q2 ) / 2
end
# build recursive function for firm2
function b2(q1,q2)
    return 100 - q1 - q2 
end

# example 1
function fixed_point_1()
    maxit = 20
    q1 = 0
    q2 = 0
    for i in 1:maxit
        q1_new = b1(q1,q2)
        q2_new = b2(q1_new,q2)
        
        # stop criteria
        if abs(q2 - q2_new) < 1e-2 && abs(q1 - q1_new) < 1e-2 || i == maxit
            println("number of iterations =",i)
            println("q1* =", q1_new)
            println("q2* =", q2_new)
            break
        end
        # update
        q1 = q1_new
        q2 = q2_new
    end
end

# example 2
function b1_2(q2)
    return 20 + q2 / 2
end
# build recursive function for firm2
function b2_2(q1)
    return 70 - 3q1/2 
end

function fixed_point_2()
    maxit = 40
    q1 = 0
    q2 = 0
    for i in 1:maxit
        q1_new = b1_2(q2)
        q2_new = b2_2(q1_new)
        println("q2_new = $q2_new" )
        # stop criteria
        if abs(q2 - q2_new) < 1e-3 || i == maxit
            println("number of iterations =",i)
            println("q1* =", q1_new)
            println("q2* =", q2_new)
            break
        end
        # update
        q1 = q1_new
        q2 = q2_new
    end
end