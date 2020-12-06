mutable struct logDoubleLogisticExpr <: ScalarDataIndex
        id::Union{Int, Nothing}
        name::String
        pos::ScalarDataIndex
        p::ScalarDataIndex
        w::ScalarDataIndex
        s::ScalarDataIndex
        hash::Union{Nothing,String}

        function logDoubleLogisticExpr( pos::ScalarDataIndex,
                        p::ScalarDataIndex,
                        w::ScalarDataIndex,
                        s::ScalarDataIndex,
                        id::Union{Int, Nothing} = nothing  )
                new(id,"logDoubleLogisticExpr", pos, p, w, s, nothing)
        end
end

function simplifyExpr(ldle::logDoubleLogisticExpr)
        return logDoubleLogisticExpr( simplifyExpr(ldle.pos), simplifyExpr(ldle.p),simplifyExpr(ldle.w),simplifyExpr(ldle.s) )
end

function sameTypeIsLess(a::logDoubleLogisticExpr, b::logDoubleLogisticExpr)
        if a.id != b.id
                return a.id < b.id
        else
                return string(a) < string(b)
        end
end

function reorder(dle::logDoubleLogisticExpr)
        reorder(dle.pos)
        reorder(dle.p)
        reorder(dle.w)
        reorder(dle.s)
        dle.hash = string(dle)
        return nothing
end

exprOrder[logDoubleLogisticExpr] = 15

function derive(e::logDoubleLogisticExpr, byX::X, byP::X)
        tempExpr = logExpr( doubleLogisticExpr(e.pos, e.p, e.w, e.s ) )
        return derive(tempExpr, byX, byP)
end

function juliaLiteral(e::logDoubleLogisticExpr, hc::hiddenCounter = hiddenCounter())

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
                eval = baseName
                hc.exprList[e.hash] = (hc.counter, eval, baseName)

                logDLFunction = string("function logDoubleLogistic(x::Float64, p::Float64, w::Float64, s::Float64)\n",
                                        "       e1 = (-x+p-w)/s\n",
                                        "       e2 = (-x+p+w)/s\n",
                                        "       c = 20.0\n",
                                        "       if x-p<0\n",
                                        "           if e1>c\n",
                                        "               #&& e2 > 2.0*c\n",
                                        "              # return -e1-log(2.0*w) +\n",
                                        "               return   -e1-log(w) + \n",
                                        "               log(0.5)+log(1.0-exp(-2.0*w/s) ) - exp(-e1) -exp(-e2)\n",
                                         #+ c + log( 0.5* (1/(1+exp(c)) -  1/(1+exp(c+2w/s) ) ) )\n",
                                        "           else\n",
                                        "               return log(0.5 * (1.0/(1.0+exp(e1)) - 1.0/(1.0+exp(e2)) )/w)\n",
                                        "           end\n",
                                        "        else\n",
                                        "           if e2<-c\n",
                                        "               return e2-log(w) + log(0.5) + log(1.0 - exp(-2.0*w/s ) )   - exp(e2) -exp(e1) \n",
                                        "               #c + log( 0.5* (1/(1+exp(c)) -  1/(1+exp(c+2w/s) ) ) )\n",
                                        "           else\n",
                                        "               return log(0.5 * (-1.0/(1.0+exp(-e1)) +  1.0/(1.0+exp(-e2)) )/w)\n",
                                        "           end\n",
                                        "       end\n",
                                        "end\n")
                logDLFunctionHash = string("function logDoubleLogistic")

                defs = string(baseName, "= logDoubleLogistic(",pos.expression,",",p.expression,",",w.expression,",",s.expression,")\n")

                if !haskey(hc.exprList, logDLFunctionHash)
                        hc.exprList[logDLFunctionHash] = (hc.counter, logDLFunction, "function")
                        defs = string( logDLFunction, defs)
                end

        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(string(pos.definitions, p.definitions, w.definitions, s.definitions, defs), eval, dataAttributes(hc))
end

function members(e::logDoubleLogisticExpr)
        return [e.pos, e.p, e.w, e.s]
end

function deRef(e::logDoubleLogisticExpr)
        pr = logDoubleLogisticExpr( constIndex(0, "dummy", 1.0),constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0) )
        pr.pos = deRef(e.pos)
        pr.p = deRef(e.p)
        pr.s = deRef(e.s)
        pr.w = deRef(e.w)

        pr.id = deepcopy(e.id)
        pr.name = deepcopy(e.name)
        pr.hash = deepcopy(e.hash)
        return pr
end
