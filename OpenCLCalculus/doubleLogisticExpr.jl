mutable struct doubleLogisticExpr <: ScalarDataIndex
        id::Union{Int, Nothing}
        name::String
        pos::ScalarDataIndex
        p::ScalarDataIndex
        w::ScalarDataIndex
        s::ScalarDataIndex
        hash::Union{Nothing,String}

        function doubleLogisticExpr(pos::ScalarDataIndex,
                                p::ScalarDataIndex,
                                w::ScalarDataIndex,
                                s::ScalarDataIndex, id::Union{Int, Nothing} = nothing)
                new(id, "doubleLogisticExpr", pos, p, w, s, nothing)
        end
end

function simplifyExpr(dle::doubleLogisticExpr)
        return doubleLogisticExpr( simplifyExpr(dle.pos), simplifyExpr(dle.p), simplifyExpr(dle.w), simplify(dle.s) )
end

function Base.:string(dle::doubleLogisticExpr)
        if dle.id==nothing
                return string("doubleLogistic(",string(dle.pos),",",string(dle.p),",",string(dle.w),",",string(dle.s),")"  )
        else
                return string("doubleLogistic", dle.id)
        end
end

function sameTypeIsLess(a::doubleLogisticExpr, b::doubleLogisticExpr)
        if a.id != b.id
                return a.id < b.id
        else
                return string(a) < string(b)
        end
end

exprOrder[doubleLogisticExpr] = 9

function reorder(dle::doubleLogisticExpr)
        reorder(dle.pos)
        reorder(dle.p)
        reorder(dle.w)
        reorder(dle.s)
        dle.hash = string(dle)
        return nothing
end



function derive(dle::doubleLogisticExpr, byX::X, byP::X)
        # function logisticDerive(dle::doubleLogisticExpr,
        #         xExpr::Vector{ScalarDataIndex},
        #         dPosdx::Vector{ScalarDataIndex},
        #         dpdx::Vector{ScalarDataIndex},
        #         dwdx::Vector{ScalarDataIndex},
        #         dsdx::Vector{ScalarDataIndex})
        #
        #
        # end

        tempExpr = prodExpr( [
                        divExpr(constIndex(0, "derDoubleLogisticConst", 0.5), dle.w),
                        sumExpr( [logisticExpr( sumExpr([dle.pos, dle.w]), dle.p , dle.s  ),
                                neg(logisticExpr( dle.pos, sumExpr([dle.p, dle.w]), dle.s  ))] )
                        ] )

        derive(tempExpr, byX, byP)
end

# function testDLE(dle::doubleLogisticExpr)
#         tempExpr = prodExpr( [
#                         divExpr(constIndex(0, "derDoubleLogisticConst", 0.5), dle.w),
#                         sumExpr( [logisticExpr( sumExpr([dle.pos, dle.w]), dle.p , dle.s  ),
#                                 neg(logisticExpr( dle.pos, sumExpr([dle.p, dle.w]), dle.s  ))] )
#                         ] )
# end


function juliaLiteral(e::doubleLogisticExpr, hc::hiddenCounter = hiddenCounter())
        if e.id == nothing
                hc.counter+=1
                baseName = string(e.name, "Anon", hc.counter)
        else
                baseName = string(e.name, e.id)
        end

        pos = juliaLiteral(e.pos,hc)
        p = juliaLiteral(e.p, hc)
        s = juliaLiteral(e.s, hc)
        w = juliaLiteral(e.w, hc)

        e.hash = string(e)
        if !haskey(hc.exprList, e.hash)

                dlFunction = string("function doubleLogistic(x::Float64, p::Float64, w::Float64, s::Float64)\n",
                                    "   e1 = (-x+p-w)/s\n",
                                    "   e2 = (-x+p+w)/s\n",
                                    "    if x-p<0\n",
                                    "       return 0.5 * (1.0/(1.0+exp(e1) )- 1.0/(1.0+exp(e2)) )/w\n",
                                    "    else\n",
                                    "        return 0.5 * (-1.0/(1.0+exp(-e1)) +  1.0/(1.0+exp(-e2)) )/w\n",
                                    "    end\n",
                                    "end\n")
                dlFunctionHash = string("function doubleLogistic")

                eval = baseName
                hc.exprList[e.hash] = (hc.counter, eval, baseName)
                #defs = string(baseName, "= ( 1.0/(1.0+exp((-",pos.expression,"-",p.expression,"-",w.expression  ,  " )/",s.expression," )) + ",
                #                        "-1.0/(1.0+exp((-",pos.expression,"-",p.expression,"+",w.expression  ,  " )/",s.expression," )) ) *0.5/",w.expression,"\n" )
                defs = string(baseName, "= doubleLogistic(",pos.expression,",",p.expression,",",w.expression,",",s.expression, ")\n")

                if !haskey(hc.exprList, dlFunctionHash)
                        hc.exprList[dlFunctionHash] = (hc.counter, dlFunction, "function")
                        defs = string(dlFunction, defs)
                end
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end

        literal(string(pos.definitions, p.definitions, w.definitions,s.definitions, defs ), eval, dataAttributes(hc))
end

function members(e::doubleLogisticExpr)
        return [e.pos, e.p, e.w, e.s]
end

function deRef(e::doubleLogisticExpr)
        pr = doubleLogisticExpr( constIndex(0, "dummy", 1.0),constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0) )
        pr.pos = deRef(e.pos)
        pr.p = deRef(e.p)
        pr.s = deRef(e.s)
        pr.w = deRef(e.w)

        pr.id = deepcopy(e.id)
        pr.name = deepcopy(e.name)
        pr.hash = deepcopy(e.hash)
        return pr
end
