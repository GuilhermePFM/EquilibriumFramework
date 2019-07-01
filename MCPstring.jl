
# OVERLOAD
import Base.string
function string(exp::JuMP.GenericAffExpr{Float64,JuMP.Variable})
    io = IOBuffer()
    println(io,exp)
    str =String(take!(io))
    return str[1:end-1]
end

# Framework methods:
import Base.print
function print(m::marketEquilibriumModel)
    str = ""
    for c in m.cstrs
        str *= string(c) * "\n"
    end
    println("Equilibrium model:")
    println(str)
end

function string(c::ComplementarityEquilibriumConstraint)
    
    str = string(c.cstr1)
    str *= " _|_ "
    str *= string(c.cstr2)

    return str    
end

function string(c::GreaterOrEqualThanEquilibriumConstraint)
    return string(sum(c.coef'*c.vars) + c.rhs) * " >= 0"
end

function string(c::LowerOrEqualThanEquilibriumConstraint)
    return "0 >= " * string(sum(c.coef'*c.vars) + c.rhs)
end

function string(c::FreeEquilibriumConstraint)
    return string(sum(c.coef'*c.vars) + c.rhs) * " free"
end

function string(c::IntervalEquilibriumConstraint)
    return string(sum(c.coef'*c.vars) + c.rhs) * " in " * "[$(c.LB)], $(c.UB)]"
end

function string(c::EqualEquilibriumConstraint)
    return "0 = " * string(sum(c.coef'*c.vars) + c.rhs)
end