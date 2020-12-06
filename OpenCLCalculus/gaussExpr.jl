mutable struct gaussExpr <: ScalarDataIndex
        id::Union{Int,Nothing}
        name::String
        pos::ScalarDataIndex
        p::ScalarDataIndex
        s::ScalarDataIndex
        hash::Union{Nothing,String}
#        nDim::Int
        function gaussExpr(pos::ScalarDataIndex, p::ScalarDataIndex, s::ScalarDataIndex, id::Union{Int,Nothing} = nothing; nDim::Int = 1)
                new(id, "gaussExpr", pos, p, s, nothing)#, nDim)
        end

end

function simplifyExpr(e::gaussExpr)
        return gaussExpr(simplifyExpr(e.pos), simplifyExpr(e.p), simplifyExpr(e.s), e.id)
end

function sameTypeIsLess(a::gaussExpr,b::gaussExpr)
        if a.id!=b.id
                return a.id < b.id
        else
                return string(a) < string(b)
        end
end

function reorder(ge::gaussExpr)
        reorder(ge.pos)
        reorder(ge.p)
        reorder(ge.s)
        ge.hash = string(ge)
        return nothing
end

function Base.:string(ge::gaussExpr)
        if ge.id==nothing
                return string("gauss(",string(ge.pos),",", string(ge.p),",",string(ge.s), ")")
        else
                return string("gaussExpr",ge.id)
        end
end

function derive(ge::gaussExpr, byX::X, byP::X)
        function gaussDerive(ge::gaussExpr, xExpr::Vector{ScalarDataIndex}, dPosdx::Vector{ScalarDataIndex},dpdx::Vector{ScalarDataIndex}, dsdx::Vector{ScalarDataIndex} )

                e1 = divExpr(sumExpr([ge.p, prodExpr([ge.pos, constIndex(0, "derGaussConst", -1.0)] )]), ge.s)
                e2 = sumExpr(
                          [divExpr(
                                prodExpr([sumExpr( [ge.pos, prodExpr([ge.p, constIndex(0, "derGaussConst", -1.0)] ) ] ),
                                         sumExpr( [ge.pos, prodExpr([ge.p, constIndex(0, "derGaussConst", -1.0)] ) ])] )
                                        , prodExpr([ge.s,ge.s, ge.s])
                                 )
                           divExpr(
                               constIndex(0, "derGaussConst", -1.0),
                               ge.s
                           )
                          ])
                for i in 1:length(xExpr)
                        a1 = nothing
                        if typeof(dPosdx[i]) == constIndex
                                if dPosdx[i].value != 0
                                        if dPosdx[i].value == 1.0
                                                a1 = e1
                                        else
                                                a1 = prodExpr([dPosdx[i], e1])
                                        end
                                end
                        else
                                a1 = prodExpr([dPosdx[i], e1])
                        end

                        a2 = nothing
                        if typeof(dpdx[i]) == constIndex
                                if dpdx[i].value != 0
                                        if dpdx[i].value == 1.0
                                                a2 = prodExpr([constIndex(0, "derGaussConst", -1.0), e1])
                                        else
                                                a2 = prodExpr([constIndex(0, "derGaussNegConst",-dpdx[i].value), e1])
                                        end
                                end
                        else
                                a2 = prodExpr([constIndex(0, "derGaussConst", -1.0), dpdx[i], e1])
                        end


                        a3 = nothing
                        if typeof(dsdx[i]) ==  constIndex
                                if dsdx[i].value != 0
                                        if dsdx[i].value == 1.0
                                                a3 = e2
                                        else
                                                a3 = prodExpr([dsdx[i], e2])
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
                                xExpr[i] = prodExpr([ge, factor[1]])
                        else
                                xExpr[i] = prodExpr([ge, sumExpr(factor)])
                        end
                end
                return xExpr
        end

        pos = derive(ge.pos, byX, byP)
        p =   derive(ge.p, byX, byP)
        s =   derive(ge.s, byX, byP)

        return gaussDerive(ge, createConstIndices(byX), pos[1], p[1], s[1]),
                gaussDerive(ge, createConstIndices(byP), pos[2], p[2], s[2])

end

exprOrder[ gaussExpr] = 7


function juliaLiteral(e::gaussExpr,hc::hiddenCounter = hiddenCounter())
        if e.id == nothing
                hc.counter+=1
                baseName = string(e.name,"Anon",hc.counter)
        else
                baseName = string(e.name,e.id)
        end
        pos = juliaLiteral(e.pos,hc)
        p = juliaLiteral(e.p,hc)
        s = juliaLiteral(e.s,hc)

        e.hash = string(e)
        if !haskey(hc.exprList, e.hash)
                eval = baseName
                hc.exprList[e.hash] = (hc.counter, eval, baseName)
                defs = string(baseName,"= 1/sqrt(2π) * exp(-0.5*((",pos.expression,"-",p.expression,")/",s.expression,")^2)/",s.expression," \n" )

        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end

        literal(string(pos.definitions, p.definitions, s.definitions,defs), eval, dataAttributes(hc))
end

function members(e::gaussExpr)
        return [e.pos, e.p, e.s]
end

function deRef(e::gaussExpr)
        pr = gaussExpr( constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0) )
        pr.pos = deRef(e.pos)
        pr.p = deRef(e.p)
        pr.s = deRef(e.s)

        pr.id = deepcopy(e.id)
        pr.name = deepcopy(e.name)
        pr.hash = deepcopy(e.hash)
        return pr
end
