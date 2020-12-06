mutable struct logExpr <: ScalarDataIndex
        id::Union{Int, Nothing}
        name::String
        pos::ScalarDataIndex
        hash::Union{Nothing,String}

        function logExpr(pos::ScalarDataIndex, id::Union{Int, Nothing} = nothing)
                new(id, "logExpr", pos, nothing)
        end
end

function simplifyExpr(le::logExpr)
        pos = simplifyExpr(le.pos)
        if typeof(pos) == gaussExpr
                return logGaussExpr(pos.pos, pos.p, pos.s)
        elseif typeof(pos) == doubleLogisticExpr
                return logDoubleLogisticExpr(pos.pos, pos.p, pos.w,pos.s)
        elseif typeof(pos) == constIndex
                return constIndex(0, "simplLogExprConst", log(pos.value))
        else
                return logExpr(pos)
        end
end

function Base.:string(le::logExpr)
        if le.id == nothing
                return string("log(",string(le.pos),")")
        else
                return string("log",le.id)
        end
end

function sameTypeIsLess(a::logExpr, b::logExpr)
        if a.id != b.id
                return a.id < b.id
        else
                return string(a) < string(b)
        end
end

exprOrder[logExpr] = 12

function reorder(le::logExpr)
        reorder(le.pos)
        le.hash = string(le)
        return nothing
end

function derive(le::logExpr, byX::X, byP::X)
        function deriveLog(le::logExpr, xExpr::Vector{T}, dPosdx::Vector{U}) where T<:ScalarDataIndex where U<:ScalarDataIndex
                c0 = constIndex(0, "derLogConst", 0.0)
                for i in 1:length(xExpr)
                        if !(typeof(dPosdx[i]) == constIndex && dPosdx[i].value == 0.0)
                                xExpr[i] = divExpr(dPosdx[i], le.pos)
                        else
                                xExpr[i] = c0
                        end
                end
                return xExpr
        end

        pos = derive(le.pos, byX, byP)
        return deriveLog(le, createConstIndices(byX), pos[1]), deriveLog(le, createConstIndices(byP), pos[2])
end

function juliaLiteral(e::logExpr, hc::hiddenCounter = hiddenCounter())
        if e.id == nothing
                hc.counter+=1
                baseName = string(e.name, "Anon", hc.counter)
        else
                baseName = string(e.name, e.id)
        end

        pos = juliaLiteral(e.pos,hc)
        e.hash = string(e)
        if !haskey(hc.exprList, e.hash)
                eval = baseName
                hc.exprList[e.hash] = (hc.counter, eval, baseName)
                defs = string(baseName, "= log(",pos.expression,")\n" )

        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end

        literal(string(pos.definitions, defs), eval, dataAttributes(hc))
end

function members(e::logExpr)
        return [e.pos]
end

function deRef(e::logExpr)
        pr = logExpr(  constIndex(0, "dummy", 1.0) )
        pr.pos = deRef(e.pos)

        pr.id = deepcopy(e.id)
        pr.name = deepcopy(e.name)
        pr.hash = deepcopy(e.hash)
        return pr
end
