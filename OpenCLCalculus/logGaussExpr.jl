mutable struct logGaussExpr <:ScalarDataIndex
        id::Union{Int, Nothing}
        name::String
        pos::ScalarDataIndex
        p::ScalarDataIndex
        s::ScalarDataIndex
        hash::Union{Nothing,String}

        function logGaussExpr( pos::ScalarDataIndex,
                        p::ScalarDataIndex,
                        s::ScalarDataIndex,
                        id::Union{Int, Nothing} = nothing  )
                new(id,"logGaussExpr", pos, p, s, nothing)
        end
end

function simplifyExpr(lge::logGaussExpr)
        pos = simplifyExpr(lge.pos)
        minusP = simplifyExpr(neg(lge.p))
        s = simplifyExpr(lge.s)
        distance = divExpr( sumExpr([ pos,minusP ]), s )

        return sumExpr( [prodExpr([constIndex(0, "simplLogGaussConstMinusHalf", -0.5),distance, distance]), neg(logExpr(s)), constIndex(0, "simplLogGaussConstLogSqrt2pi", -0.5*log(2*π) )  ] )
end

function reorder(e::logGaussExpr)
        reorder(e.pos)
        reorder(e.p)
        reorder(e.s)
        e.hash = string(e)
        return nothing
end

function Base.:string(lge::logGaussExpr)
        #I call it lognormal, so it is visually different than the log(gauss()) hash
        if lge.id==nothing
                return string("logNormal(",string(lge.pos),",",string(lge.p),",",string(lge.s),")"  )
        else
                return string("logNormal", lge.id)
        end
end

function sameTypeIsLess(a::logGaussExpr, b::logGaussExpr)
        if a.id != b.id
                return a.id < b.id
        else
                return string(a) < string(b)
        end
end

exprOrder[logGaussExpr] = 14

function derive(lge::logGaussExpr, byX::X, byP::X)
        simplExpr = simplifyExpr(lge)
        return derive(simplExpr, byX, byP)
end

function juliaLiteral(e::logGaussExpr, hc::hiddenCounter = hiddenCounter())
        if e.id == nothing
                hc.counter+=1
                baseName = string(e.name, "Anon", hc.counter)
        else
                baseName = string(e.name, e.id)
        end

        pos = juliaLiteral(e.pos,hc)
        p = juliaLiteral(e.p, hc)
        s = juliaLiteral(e.s, hc)

        e.hash = string(e)
        if !haskey(hc.exprList, e.hash)
                eval = baseName
                hc.exprList[e.hash] = (hc.counter, eval, baseName)
                defs = string(baseName, "= -0.5*((",pos.expression,"-",p.expression,")/",s.expression,")^2 - log(",s.expression,") - 0.5*log(2.0π)\n")
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(string(pos.definitions, p.definitions, s.definitions, defs), eval, dataAttributes(hc))
end

function members(e::logGaussExpr)
        return [e.pos, e.p, e.s]
end

function deRef(e::logGaussExpr)
        pr = logGaussExpr( constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0) )
        pr.pos = deRef(e.pos)
        pr.p = deRef(e.p)
        pr.s = deRef(e.s)

        pr.id = deepcopy(e.id)
        pr.name = deepcopy(e.name)
        pr.hash = deepcopy(e.hash)
        return pr
end
