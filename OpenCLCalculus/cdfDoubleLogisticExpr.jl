mutable struct cdfDoubleLogisticExpr <: ScalarDataIndex
        id::Union{Int, Nothing}
        name::String
        pos::ScalarDataIndex
        p::ScalarDataIndex
        w::ScalarDataIndex
        s::ScalarDataIndex
        hash::Union{Nothing,String}

        function cdfDoubleLogisticExpr(pos::ScalarDataIndex,
                                p::ScalarDataIndex,
                                w::ScalarDataIndex,
                                s::ScalarDataIndex, id::Union{Int, Nothing} = nothing)
                new(id, "cdfDoubleLogisticExpr", pos, p, w, s, nothing)
        end
end

function simplifyExpr(dle::cdfDoubleLogisticExpr)
        return cdfDoubleLogisticExpr( simplifyExpr(dle.pos), simplifyExpr(dle.p), simplifyExpr(dle.w), simplifyExpr(dle.s) )
end

function Base.:string(dle::cdfDoubleLogisticExpr)
        if dle.id==nothing
                return string("cdfDoubleLogistic(",string(dle.pos),",",string(dle.p),",",string(dle.w),",",string(dle.s),")"  )
        else
                return string("cdfDoubleLogistic", dle.id)
        end
end

function sameTypeIsLess(a::cdfDoubleLogisticExpr, b::cdfDoubleLogisticExpr)
        if a.id != b.id
                return a.id < b.id
        else
                return string(a) < string(b)
        end
end

exprOrder[cdfDoubleLogisticExpr] = 11

function reorder(dle::cdfDoubleLogisticExpr)
        reorder(dle.pos)
        reorder(dle.p)
        reorder(dle.w)
        reorder(dle.s)
        dle.hash = string(dle)
        return nothing
end

function derive(dle::cdfDoubleLogisticExpr, byX::X, byP::X)
        function cdfDoubleLogisticDerive(dle::cdfDoubleLogisticExpr,
                xExpr::Vector{T},
                dPosdx::Vector{U},
                dpdx::Vector{V},
                dwdx::Vector{W},
                dsdx::Vector{Z}) where T<:ScalarDataIndex where U<:ScalarDataIndex where V<:ScalarDataIndex where W<:ScalarDataIndex where Z<:ScalarDataIndex

                cm1 = constIndex(0, "derCDFDoubleLogistic", -1.0)
                chalf = constIndex(0, "derHalfCDFDoubleLogistic", 0.5)



                cdf = cdfDoubleLogisticExpr(dle.pos, dle.p, dle.w, dle.s)
                pdf = doubleLogisticExpr(dle.pos, dle.p, dle.w, dle.s)
                l1 = logisticExpr(sumExpr([dle.pos, dle.w]), dle.p, dle.s)
                l2 = logisticExpr(dle.pos, sumExpr([dle.p, dle.w]), dle.s )
                e1 = divExpr(sumExpr(  [ neg(dle.pos), dle.p, neg(dle.w)   ]  ) , dle.s)
                e2 = divExpr(sumExpr(  [ dle.pos, neg(dle.p), neg(dle.w)   ]  ) , dle.s)

                e3 = divExpr( sumExpr([prodExpr([chalf, l1]), prodExpr([chalf,l2]), prodExpr([cm1,cdf]) ]),dle.w  )
                e4 = sumExpr( [
                        divExpr(cdf,dle.s),
                        divExpr(sumExpr([
                                prodExpr([chalf, l1, e1]),
                                prodExpr([chalf, l2, e2])
                                ]) , dle.w )
                        ] )

                for i in 1:length(xExpr)
                        a1 = nothing
                        if typeof(dPosdx[i]) == constIndex
                                if dPosdx[i].value != 0
                                        if dPosdx[i].value == 1.0
                                                a1 = pdf
                                        else
                                                a1 = prodExpr([dPosdx[i], pdf])
                                        end
                                end
                        else
                                a1 = prodExpr([dPosdx[i], f2])
                        end

                        a2 = nothing
                        if typeof(dpdx[i]) == constIndex
                                if dpdx[i].value != 0
                                        a2 = prodExpr([constIndex(0, "derNegCDFDoubleLogistic", -dpdx[i].value), pdf])
                                end
                        else
                                a2 = prodExpr([cm1, dpdx[i], pdf])
                        end

                        a3 = nothing
                        if typeof(dwdx[i]) == constIndex
                                if dwdx[i].value != 0
                                        if dwdx[i].value == 1
                                                a3 = e3
                                        else
                                                a3 = prodExpr([dwdx[i], e3])
                                        end
                                end
                        else
                                a3 = prodExpr([dwdx[i], e3])
                        end

                        a4 = nothing
                        if typeof(dsdx[i]) == constIndex
                                if dsdx[i].value != 0
                                        if dsdx[i].value == 1
                                                a4 = e4
                                        else
                                                a4 = prodExpr([dsdx[i], e4])
                                        end
                                end
                        else
                                a4 = prodExpr([dsdx[i], e4])
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

                        if a4!= nothing
                                push!(factor, a4)
                        end

                        if length(factor) == 0
                                xExpr[i] = constIndex(0, "derivConst", 0.0)
                        elseif length(factor) == 1
                                xExpr[i] = factor[1]
                        else
                                xExpr[i] = sumExpr(factor)
                        end
                end
                return xExpr

        end

        pos = derive(dle.pos, byX, byP)
        p = derive(dle.p, byX, byP)
        w = derive(dle.w, byX, byP)
        s = derive(dle.s, byX, byP)

        return cdfDoubleLogisticDerive(dle, createConstIndices(byX), pos[1], p[1], w[1], s[1] ),
                cdfDoubleLogisticDerive(dle, createConstIndices(byP), pos[2], p[2], w[2], s[2] )
end

function juliaLiteral(e::cdfDoubleLogisticExpr, hc::hiddenCounter = hiddenCounter())
        if e.id == nothing
                hc.counter+=1
                baseName = string(e.name, "Anon", hc.counter)
        else
                baseName = string(e.name, e.id)
        end

        pos = juliaLiteral(e.pos,hc)
        p = juliaLiteral(e.p, hc)
        w = juliaLiteral(e.w, hc)
        s = juliaLiteral(e.s, hc)

        e.hash = string(e)
        if !haskey(hc.exprList, e.hash)
#there are several equvivalent parameterizations
#eg.
#exp(2w)  + (1-exp(2w))/(1+exp(x-w)) == (1+exp(x+w))/(1+exp(x-w))
#https://www.wolframalpha.com/input/?i=exp%282w%29++%2B+%281-exp%282w%29%29%2F%281%2Bexp%28x-w%29%29+%3D+%281%2Bexp%28x%2Bw%29%29%2F%281%2Bexp%28x-w%29%29
                # cdfDLFunction = string("function cdfDoubleLogistic(x::Float64, p::Float64, w::Float64, s::Float64)\n",
                #                         "       e1 = (x-p+w)/s\n",
                #                         "       e2 = (x-p-w)/s\n",
                #                         "       if e1>700.0\n",
                #                         "            return 0.5\n",
                #                         "        else\n",
                #                         "            return 0.5*s*log( (1.0+exp(e1) )/(1.0+exp(e2) ) )/w - 0.5\n",
                #                         "#           return 0.5*s*log(exp(2*w) + (1-exp(2*w))/(1+exp((x-p-w)/s)))/w\n",
                #                         "        end\n",
                #                         "end\n")
                cdfDLFunction = string("function cdfDoubleLogistic(x::Float64, p::Float64, w::Float64, s::Float64)
                       e1 = (x-p+w)/s
                       e2 = (x-p-w)/s

                       if e1 < 10.0
                           return 0.5*s*log( (1.0+exp(e1) )/(1.0+exp(e2) ) )/w - 0.5
                       end
                       if e1>=10.0 && e2 <-10.0
                           return 0.5*(x-p)/w
                        end
                        if e2>=-10.0 && e2<10.0
                            return 0.5*s*log(   1.0 /(1.0+exp(e2) )  )/w  + 0.5*(x-p)/w
                        end
                        if e2>=10.0
                            return 0.5
                        end
                end\n")

                cdfDLFunctionHash = string("function cdfDoubleLogistic")

                eval = baseName
                hc.exprList[e.hash] = (hc.counter, eval, baseName)
                #defs = string(baseName,"= 0.5*",s.expression,"*log( (1+exp( (",pos.expression,"-",p.expression,"+",w.expression," )/",s.expression," ))/ ",
                #"(1+exp( (",pos.expression,"-",p.expression,"-",w.expression," )/",s.expression," ) ) )/(",w.expression,")\n")
                defs = string(baseName,"= cdfDoubleLogistic(",pos.expression,",",p.expression,",",w.expression,",",s.expression,")\n" )

                if !haskey(hc.exprList, cdfDLFunctionHash)
                        hc.exprList[cdfDLFunctionHash] = (hc.counter, cdfDLFunction, "function")
                        defs = string( cdfDLFunction, defs)
                end
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end

        literal(string(pos.definitions, p.definitions, w.definitions, s.definitions, defs ), eval, dataAttributes(hc))
end

function members(e::cdfDoubleLogisticExpr)
        return [e.pos, e.p, e.w, e.s]
end

function deRef(e::cdfDoubleLogisticExpr)
        pr = cdfDoubleLogisticExpr( constIndex(0, "dummy", 1.0),constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0) )
        pr.pos = deRef(e.pos)
        pr.p = deRef(e.p)
        pr.s = deRef(e.s)
        pr.w = deRef(e.w)

        pr.id = deepcopy(e.id)
        pr.name = deepcopy(e.name)
        pr.hash = deepcopy(e.hash)
        return pr
end
