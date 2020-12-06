mutable struct erfExpr <: ScalarDataIndex
        id::Union{Int, Nothing}
        name::String
        pos::ScalarDataIndex
        p::ScalarDataIndex
        s::ScalarDataIndex
        hash::Union{Nothing,String}

        function erfExpr(pos::ScalarDataIndex,
                                p::ScalarDataIndex,
                                s::ScalarDataIndex, id::Union{Int, Nothing} = nothing)
                new(id, "erfExpr", pos, p, s, nothing)
        end
end

function simplifyExpr(ee::erfExpr)
        return erfExpr(simplifyExpr(ee.pos), simplifyExpr(ee.p), simplifyExpr(ee.s))
end

function Base.:string(ee::erfExpr)
        if ee.id==nothing
                return string("erf(",string(ee.pos),",",string(ee.p),",",string(ee.s),")"  )
        else
                return string("erf", ee.id)
        end
end

function sameTypeIsLess(a::erfExpr, b::erfExpr)
        if a.id != b.id
                return a.id < b.id
        else
                return string(a) < string(b)
        end
end

function reorder(ee::erfExpr)
        reorder(ee.pos)
        reorder(ee.p)
        reorder(ee.s)
        ee.hash = string(ee)
        return nothing
end

exprOrder[erfExpr] = 10

function derive(ee::erfExpr, byX::X, byP::X)
        function deriveErf(ee::erfExpr,
                xExpr::Vector{T},
                dPosdx::Vector{U},
                dpdx::Vector{V},
                dsdx::Vector{W})  where T<:ScalarDataIndex where U<:ScalarDataIndex where V<:ScalarDataIndex where W<:ScalarDataIndex

                f1 = gaussExpr(ee.pos, ee.p, ee.s)
                e1 = divExpr(sumExpr([ee.p, neg(ee.pos) ]), ee.s)

                for i in 1:length(xExpr)
                        a1 = nothing
                        if typeof(dPosdx[i]) == constIndex
                                if dPosdx[i].value != 0.0
                                        a1 = dPosdx[i]
                                end
                        else
                                a1 = dPosdx[i]
                        end


                        a2 = nothing
                        if typeof(dpdx[i]) == constIndex
                                if dpdx[i].value != 0.0
                                        a2 = constIndex(0, "derErfNegConst", -dpdx[i].value)
                                end
                        else
                                a2 = prodExpr( [constIndex(0, "derErfConst", -1.0), dpdx[i]] )
                        end

                        a3 = nothing
                        if typeof(dsdx[i]) == constIndex
                                if dsdx[i].value != 0.0
                                        if dsdx[i].value == 1.0
                                                a3 = e1
                                        else
                                                a3 = prodExpr([ dsdx[i], e1])
                                        end
                                end
                        else
                                a3 = prodExpr( [dsdx[i], e1])

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

        pos = derive(ee.pos, byX, byP)
        p = derive(ee.p, byX, byP)
        s = derive(ee.s, byX, byP)

        return deriveErf(ee, createConstIndices(byX), pos[1], p[1], s[1]),
                deriveErf(ee, createConstIndices(byP), pos[2], p[2], s[2])
end

function juliaLiteral(e::erfExpr, hc::hiddenCounter = hiddenCounter())
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
                defs = string(baseName,"= 0.5*erf( (",pos.expression,"-",p.expression,")/sqrt(2.0)/",s.expression,")\n"  )
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end

        literal(string(pos.definitions, p.definitions, s.definitions, defs ), eval, dataAttributes(hc))
end


function members(e::erfExpr)
        return [e.pos, e.p, e.s]
end

function deRef(e::erfExpr)
        pr = erfExpr( constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0) )
        pr.pos = deRef(e.pos)
        pr.p = deRef(e.p)
        pr.s = deRef(e.s)

        pr.id = deepcopy(e.id)
        pr.name = deepcopy(e.name)
        pr.hash = deepcopy(e.hash)
        return pr
end
