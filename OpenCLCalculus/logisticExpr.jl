mutable struct logisticExpr<:ScalarDataIndex
        id::Union{Int, Nothing}
        name::String
        pos::ScalarDataIndex
        p::ScalarDataIndex
        s::ScalarDataIndex
        hash::Union{Nothing,String}

        function logisticExpr( pos::ScalarDataIndex,
                        p::ScalarDataIndex,
                        s::ScalarDataIndex,
                        id::Union{Int, Nothing} = nothing  )
                new(id,"logisticExpr", pos, p, s, nothing)
        end
end

function simplifyExpr(le::logisticExpr)
        return logisticExpr( simplifyExpr(le.pos), simplifyExpr(le.p), simplifyExpr(le.s) )
end

function Base.:string(le::logisticExpr)
        if le.id==nothing
                return string("logistic(",string(le.pos),",",string(le.p),",",string(le.s),")"  )
        else
                return string("logistic", le.id)
        end
end

function sameTypeIsLess(a::logisticExpr, b::logisticExpr)
        if a.id != b.id
                return a.id < b.id
        else
                return string(a) < string(b)
        end
end

function reorder(le::logisticExpr)
        reorder(le.pos)
        reorder(le.p)
        reorder(le.s)
        le.hash = string(le)
        return nothing
end

exprOrder[logisticExpr] = 8

function derive(le::logisticExpr, byX::X, byP::X)
        function logisticDerive(le::logisticExpr,
                xExpr::Vector{T},
                dPosdx::Vector{U},
                dpdx::Vector{V},
                dsdx::Vector{W}) where T<:ScalarDataIndex where U<:ScalarDataIndex where V<:ScalarDataIndex where W<:ScalarDataIndex


                #f1 = prodExpr([logisticExpr( le.pos, le.p, le.s ), logisticExpr(le.p, le.pos, le.s)] )
                cm1 = constIndex(0, "derLogisticConst", -1.0)
                f1 = sumExpr([le, prodExpr([cm1, le, le] )])

                e1 = divExpr(sumExpr([le.p, prodExpr([le.pos, cm1] )]), prodExpr([le.s, le.s]))

                for i in 1:length(xExpr)
                        a1 = nothing
                        if typeof(dPosdx[i]) == constIndex
                                if dPosdx[i].value != 0.0
                                        a1 = divExpr(dPosdx[i], le.s)
                                end
                        else
                                a1 = divExpr(dPosdx[i], f1)
                        end

                        a2 = nothing
                        if typeof(dpdx[i]) == constIndex
                                if dpdx[i].value != 0.0
                                        a2 = divExpr(constIndex(0, "derLogisticNegConst", -dpdx[i].value), le.s)
                                end
                        else
                                a2 = divExpr(prodExpr([cm1, dpdx[i]]), le.s)
                        end

                        a3 = nothing
                        if typeof(dsdx[i]) ==  constIndex
                                if dsdx[i].value != 0
                                        if dsdx[i].value == 1.0
                                                a3 = e1
                                        else
                                                a3 = prodExpr([dsdx[i], e1])
                                        end
                                end
                        else
                                a3 = prodExpr([dsdx[i], e2])
                        end

                        factor = ScalarDataIndex[]
                        if a1!= nothing
                                push!(factor, a1)
                        end

                        if a2!= nothing
                                push!(factor, a2)
                        end

                        if a3!= nothing
                                push!(factor, a3)
                        end

                        if length(factor) == 0
                                xExpr[i] = constIndex(0, "derivConst", 0.0)
                        elseif length(factor) == 1
                                xExpr[i] = prodExpr([f1, factor[1]])
                        else
                                xExpr[i] = prodExpr([f1, sumExpr(factor)])
                        end
                end
                return xExpr
        end

        pos = derive(le.pos, byX, byP)
        p = derive(le.p, byX, byP)
        s = derive(le.s, byX, byP)

        return logisticDerive(le, createConstIndices(byX), pos[1], p[1], s[1]),
                logisticDerive(le, createConstIndices(byP), pos[2], p[2], s[2])

end

function juliaLiteral(e::logisticExpr, hc::hiddenCounter = hiddenCounter())
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
                logisticFunction = string("function logistic(x::Float64, p::Float64, s::Float64)\n",
                                          "     e1 = (-x+p)/s\n",
                                          "     return 1.0/(1.0+exp(e1))\n",
                                          "end\n")
                logisticFunctionHash = string("function logistic")
                #defs = string(baseName, "= 1.0/(1+exp((-",pos.expression,"+",p.expression,")/",s.expression,"))\n" )
                defs = string(baseName, "= logistic(",pos.expression,",",p.expression,",",s.expression,")\n" )

                if !haskey(hc.exprList, logisticFunctionHash)
                        hc.exprList[logisticFunctionHash] = (hc.counter, logisticFunction, "function")
                        defs = string(logisticFunction, defs)
                end

        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end

        literal(string(pos.definitions, p.definitions, s.definitions, defs), eval, dataAttributes(hc))
end

function members(e::logisticExpr)
        return [e.pos, e.p, e.s]
end

function deRef(e::logisticExpr)
        pr = logisticExpr( constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0) )
        pr.pos = deRef(e.pos)
        pr.p = deRef(e.p)
        pr.s = deRef(e.s)

        pr.id = deepcopy(e.id)
        pr.name = deepcopy(e.name)
        pr.hash = deepcopy(e.hash)
        return pr
end
