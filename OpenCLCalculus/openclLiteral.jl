function firstPass(c::constIndex, hc::hiddenCounter = hiddenCounter())
        c.hash = string(c)
        if !haskey(hc.exprList, c.hash)
                #const are not so important, maybe they don't need to be tracked
                hc.exprList[c.hash] = (0,c.name, string(c.value) )
        else
                hc.exprList[c.hash].count += 1
        end
        return nothing
end

function firstPass(d::dataIndex, hc::hiddenCounter = hiddenCounter())
        d.hash = string(d)
        if !haskey(hc.exprList, d.hash)
                hc.exprList[d.hash] = (-d.id,d.name, string(d.attribute)) #I want to track if a dataIndex was used or is obsolete
        else
                hc.exprList[d.hash].count += 1
        end
        return nothing
end

function firstPass(dv::Vector{T}, hc::hiddenCounter = hiddenCounter()) where T<:Union{ScalarDataIndex, derivedParametericAvg}
        for d in dv
                firstPass(d, hc)
        end
        return nothing
end


function firstPass(e::ScalarDataIndex, hc::hiddenCounter = hiddenCounter())
        e.hash = string(e)

        if !haskey(hc.exprList, e.hash)
                for p in members(e)
                        firstPass(p, hc)
                end
                if e.id == nothing
                        hc.counter+=1
                        baseName = string(e.name,"Anon",hc.counter)
                else
                        baseName = string(e.name,e.id)
                end
                hc.exprList[e.hash] = (hc.counter, baseName, "expr")

        else
                #no call for firstPass, as in principle everything under members(e) were described once already
                #be careful though, those are still in the tree, without hash
                hc.exprList[e.hash].count += 1
        end
        return nothing
end

function firstPass(e::prodExpr, hc::hiddenCounter = hiddenCounter())
        e.hash = string(e)

        if !haskey(hc.exprList, e.hash)
                for p in members(e)
                        firstPass(p, hc)
                end
                if e.id == nothing
                        hc.counter+=1
                        baseName = string(e.name,"Anon",hc.counter)
                else
                        baseName = string(e.name,e.id)
                end
                hc.exprList[e.hash] = (hc.counter, baseName, "expr")


                hc.prodList[e.hash] = [i.hash for i in e.px] #unique ?

        else
                #no call for firstPass, as in principle everything under members(e) were described once already
                #be careful though, those are still in the tree, without hash
                hc.exprList[e.hash].count += 1
        end
        return nothing
end


function firstPass(e::derivedParametericAvg, hc::hiddenCounter = hiddenCounter())
        e.hash = string(e)

        if !haskey(hc.exprList, e.hash)
                for p in members(e)
                        firstPass(p, hc)
                end
                if e.id == nothing
                        hc.counter+1
                        baseName = string(e.name, "Anon", hc.counter)
                else
                        baseName = string(e.name,e.id)
                end
                hc.exprList[e.hash] = (hc.counter, baseName, "parExpr")
        else
                hc.exprList[e.hash].count += 1
        end
        return nothing
end

function secondPass(e::ScalarDataIndex, groupingCommands::Dict{String,Array{Tuple{String,String},1}})
        e.hash = string(e)
        for p in members(e)
                secondPass(p,groupingCommands)
        end
        return nothing
end

function secondPass(e::derivedParametericAvg, groupingCommands::Dict{String,Array{Tuple{String,String},1}})
        for p in members(e)
                secondPass(p,groupingCommands)
        end
        return nothing
end

function secondPass(p::prodExpr, groupingCommands::Dict{String,Array{Tuple{String,String},1}})
        for pm in p.px
                #@show pm
                if pm.hash == nothing
                        pm.hash = string(pm)
                        println("px.hash == nothing -> ", pm.hash)
                end
        end
        if p.hash == nothing
                p.hash = string(p)
                println("p.hash == nothing -> ", p.hash)
        end

        println("-> before", p.hash)
        #p.px = deepcopy(p.px) # dereferencing is necessary
        applyRemove(p, groupingCommands)
        println("<- after ", p.hash)
        for pm in members(p)
                secondPass(pm,groupingCommands)
        end
        return nothing
end

function secondPass(pv::Vector{T}, groupingCommands::Dict{String,Array{Tuple{String,String},1}}) where T<:Union{ScalarDataIndex, derivedParametericAvg}
        for pm in members(pv)
                secondPass(pm,groupingCommands)
        end
        return nothing
end

function openclLiteral(e::constIndex, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(e, hc)
        end
        literal("", string(e.value,"f"), dataAttributes(hc) )
end

function openclLiteral(d::dataIndex, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(d, hc)
        end
        #literal("",string(d),dataAttributes(hc))
        literal("",string(d.name,"[",d.id-1,"]"),dataAttributes(hc))
end

function openclLiteral(e::sumExpr, hc::hiddenCounter = hiddenCounter())

        if hc.counter == 0
                firstPass(e,hc)
        end


        if hc.exprList[e.hash].count > 0
                # if length(e.px) == 1
                #         myExpression = pxv[1].expression
                # else
                #
                # end
                pxv = [openclLiteral(p, hc) for p in e.px] # calling openclLiteral for each element, not for the vector!
                pxvDefs = string([ l.definitions for l in pxv ]...)

                myExpression = string([string(l.expression, "+") for l in pxv[1:end-1]]...,pxv[end].expression )
                if hc.exprList[e.hash].count > 1
                        eval = hc.exprList[e.hash].eval
                        defs = string("float ", eval, "=", myExpression, ";\n" )
                else
                        defs = ""
                        if length(e.px) == 1
                                eval = myExpression
                        else
                                eval = string("(",myExpression, ")" )
                        end
                end
                hc.exprList[e.hash].count = -hc.exprList[e.hash].count
                defs = string(pxvDefs, defs)
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal( defs, eval, dataAttributes(hc) )
end

function openclLiteral(e::prodExpr, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(e,hc)
        end


        if hc.exprList[e.hash].count > 0
                # if length(e.px) == 1
                #         myExpression = pxv[1].expression
                # else
                #
                # end
                pxv = [openclLiteral(p, hc) for p in e.px] # calling openclLiteral for each element, not for the vector!
                pxvDefs = string([ l.definitions for l in pxv ]...)

                myExpression = string([string(l.expression, "*") for l in pxv[1:end-1]]...,pxv[end].expression )
                if hc.exprList[e.hash].count > 1
                        eval = hc.exprList[e.hash].eval
                        defs = string("float ",eval, "=", myExpression, ";\n" )
                else
                        defs = ""
                        if length(e.px) == 1
                                eval = myExpression
                        else
                                eval = string("(",myExpression, ")" )
                        end
                end
                defs = string(pxvDefs, defs)
                hc.exprList[e.hash].count = -hc.exprList[e.hash].count
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(defs, eval, dataAttributes(hc) )
end

function openclLiteral(e::dotProdExpr, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(e,hc)
        end


        if hc.exprList[e.hash].count > 0
                posv = [openclLiteral(p, hc) for p in e.pos] # calling openclLiteral for each element, not for the vector!
                posDefs = string([ l.definitions for l in posv ]...)
                parv = [openclLiteral(p, hc) for p in e.par] # calling openclLiteral for each element, not for the vector!
                parDefs = string([ l.definitions for l in parv ]...)
                nDim = length(e.pos)

                stringVec = [ string(posv[i].expression, "*", parv[i].expression)  for i in 1:length(posv)  ]
                myExpression = string([string(str, "+") for str in stringVec[1:end-1]]...,stringVec[end] )

                if hc.exprList[e.hash].count > 1
                        eval = hc.exprList[e.hash].eval
                        defs = string("float ",eval, "=", myExpression, ";\n" )
                else
                        defs = ""
                        if nDim == 1
                                eval = myExpression
                        else
                                eval = string("(",myExpression, ")" )
                        end
                end
                defs = string(posDefs, parDefs, defs)
                hc.exprList[e.hash].count = -hc.exprList[e.hash].count
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(defs, eval, dataAttributes(hc) )
end

function openclLiteral(e::gaussExpr, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(e,hc)
        end


        if hc.exprList[e.hash].count > 0
                pos = openclLiteral(e.pos, hc)
                p = openclLiteral(e.p, hc)
                s = openclLiteral(e.s, hc)

                gaussFunctionDeclaration = string("float gauss(float x, float p, float s);\n")
                gaussFunction = string(
                        "float gauss(float x, float p, float s) {\n",
                        "        float distance = (x-p)/s;\n",
                        "        return 0.3989422804014327f * exp((float)( -0.5f*distance*distance))/s;\n",
                        "}\n\n")

                gaussFunctionHash = string("gauss")
                hc.exprList[gaussFunctionHash] = (-5, gaussFunction,"function")
                hc.exprList[gaussFunctionHash].count = 0
                gaussFunctionDeclarationHash = string("gaussdeclaration")
                hc.exprList[gaussFunctionDeclarationHash] = (-137, gaussFunctionDeclaration,"function declaration")
                hc.exprList[gaussFunctionDeclarationHash].count = 0


                myExpression = string("gauss(",pos.expression,",",p.expression,",",s.expression,")" )

                if hc.exprList[e.hash].count > 1
                        eval = hc.exprList[e.hash].eval
                        defs = string("float ",eval, "=", myExpression, ";\n" )
                else
                        defs = ""
                        eval = string("(",myExpression, ")" )
                end
                defs = string(pos.definitions, p.definitions, s.definitions, defs)
                hc.exprList[e.hash].count = -hc.exprList[e.hash].count
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(defs, eval, dataAttributes(hc) )
end

function openclLiteral(e::logGaussExpr, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(e,hc)
        end

        if hc.exprList[e.hash].count > 0
                pos = openclLiteral(e.pos, hc)
                p = openclLiteral(e.p, hc)
                s = openclLiteral(e.s, hc)

                myExpression = string("-0.5*(",pos.expression,"-",p.expression,")/",s.expression,"*(",pos.expression,"-",p.expression,")/",s.expression," - log((float)(",s.expression,")) - 0.9189385332046727f")

                if hc.exprList[e.hash].count > 1
                        eval = hc.exprList[e.hash].eval
                        defs = string("float ",eval, "=", myExpression, ";\n" )
                else
                        defs = ""
                        eval = string("(",myExpression, ")" )
                end
                defs = string(pos.definitions, p.definitions, s.definitions, defs)
                hc.exprList[e.hash].count = -hc.exprList[e.hash].count
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(defs, eval, dataAttributes(hc) )
end

function openclLiteral(e::erfExpr, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(e,hc)
        end

        if hc.exprList[e.hash].count > 0
                pos = openclLiteral(e.pos, hc)
                p = openclLiteral(e.p, hc)
                s = openclLiteral(e.s, hc)

                myExpression = string("0.5*erf( (float)((",pos.expression,"-",p.expression,")*0.7071067811865475f/",s.expression,"))")

                if hc.exprList[e.hash].count > 1
                        eval = hc.exprList[e.hash].eval
                        defs = string("float ",eval, "=", myExpression, ";\n" )
                else
                        defs = ""
                        eval = string("(",myExpression, ")" )
                end
                defs = string(pos.definitions, p.definitions, s.definitions, defs)
                hc.exprList[e.hash].count = -hc.exprList[e.hash].count
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(defs, eval, dataAttributes(hc) )
end

function openclLiteral(e::divExpr, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(e,hc)
        end

        if hc.exprList[e.hash].count > 0
                nom = openclLiteral(e.nominator, hc)
                den = openclLiteral(e.denominator, hc)

                myExpression = string(nom.expression,"/",den.expression )

                if hc.exprList[e.hash].count > 1
                        eval = hc.exprList[e.hash].eval
                        defs = string("float ",eval, "=", myExpression, ";\n" )
                else
                        defs = ""
                        eval = string("(",myExpression, ")" )
                end
                defs = string(nom.definitions, den.definitions, defs)
                hc.exprList[e.hash].count = -hc.exprList[e.hash].count
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(defs, eval, dataAttributes(hc) )
end

function openclLiteral(e::logExpr, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(e,hc)
        end

        if hc.exprList[e.hash].count > 0
                pos = openclLiteral(e.pos, hc)

                myExpression = string("log((float)(",pos.expression,"))" )

                if hc.exprList[e.hash].count > 1
                        eval = hc.exprList[e.hash].eval
                        defs = string("float ",eval, "=", myExpression, ";\n" )
                else
                        defs = ""
                        eval = string(myExpression)
                end
                defs = string(pos.definitions, defs)
                hc.exprList[e.hash].count = -hc.exprList[e.hash].count
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(defs, eval, dataAttributes(hc) )
end

function openclLiteral(e::logisticExpr, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(e,hc)
        end


        if hc.exprList[e.hash].count > 0
                pos = openclLiteral(e.pos, hc)
                p = openclLiteral(e.p, hc)
                s = openclLiteral(e.s, hc)


                logisticFunctionDeclaration = string("float logistic(float x, float p, float s);\n")
                logisticFunction = string(
                                          "float logistic(float x, float p, float s){\n",
                                          "     float e1 = (-x+p)/s;\n",
                                          "     return 1.0/(1.0+exp((float)(e1)));\n",
                                          "}\n")
                logisticFunctionHash = string("logistic")
                hc.exprList[logisticFunctionHash] = (-1, logisticFunction, "function")
                hc.exprList[logisticFunctionHash].count = 0
                logisticFunctionDeclarationHash = string("logisticdeclaration")
                hc.exprList[logisticFunctionDeclarationHash] = (-137, logisticFunctionDeclaration,"function declaration")
                hc.exprList[logisticFunctionDeclarationHash].count = 0

                myExpression = string("logistic(",pos.expression,",",p.expression,",",s.expression,")")

                if hc.exprList[e.hash].count > 1
                        eval = hc.exprList[e.hash].eval
                        defs = string("float ",eval, "=", myExpression, ";\n" )
                else
                        defs = ""
                        eval = string(myExpression )
                end
                defs = string(pos.definitions, p.definitions, s.definitions, defs)
                hc.exprList[e.hash].count = -hc.exprList[e.hash].count
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(defs, eval, dataAttributes(hc) )
end

function openclLiteral(e::doubleLogisticExpr, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(e,hc)
        end



        if hc.exprList[e.hash].count > 0
                pos = openclLiteral(e.pos, hc)
                p = openclLiteral(e.p, hc)
                w = openclLiteral(e.w, hc)
                s = openclLiteral(e.s, hc)

                doubleLogisticFunctionDeclaration = string("float doubleLogistic(float x, float p, float w, float s);\n")
                doubleLogisticFunction = string(
                                          "float doubleLogistic(float x, float p, float w, float s){\n",
                                          "     float e1 = (-x+p-w)/s;\n",
                                          "     float e2 = (-x+p+w)/s;\n",
                                          "     if (x-p<0.0f) {\n",
                                          "        return 0.5f * (1.0f/(1.0f+exp((float)(e1)) )- 1.0f/(1.0f+exp((float)(e2))) )/w;\n",
                                          "     } else {\n",
                                          "         return 0.5f * (-1.0f/(1.0f+exp(-(float)(e1))) +  1.0f/(1.0f+exp(-(float)(e2))) )/w;\n",
                                          "     }\n",
                                          "}\n")
                doubleLogisticFunctionHash = string("doubleLogistic")
                hc.exprList[doubleLogisticFunctionHash] = (-2, doubleLogisticFunction, "function")
                hc.exprList[doubleLogisticFunctionHash].count = 0
                doubleLogisticFunctionDeclarationHash = string("doubleLogisticdeclaration")
                hc.exprList[doubleLogisticFunctionDeclarationHash] = (-137, doubleLogisticFunctionDeclaration,"function declaration")
                hc.exprList[doubleLogisticFunctionDeclarationHash].count = 0

                myExpression = string("doubleLogistic(",pos.expression,",",p.expression,",",w.expression,",",s.expression,")")

                if hc.exprList[e.hash].count > 1
                        eval = hc.exprList[e.hash].eval
                        defs = string("float ",eval, "=", myExpression, ";\n" )
                else
                        defs = ""
                        eval = string(myExpression )
                end
                defs = string(pos.definitions, p.definitions, w.definitions, s.definitions, defs)
                hc.exprList[e.hash].count = -hc.exprList[e.hash].count
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(defs, eval, dataAttributes(hc) )
end

function openclLiteral(e::cdfDoubleLogisticExpr, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(e,hc)
        end

        if hc.exprList[e.hash].count > 0
                pos = openclLiteral(e.pos, hc)
                p = openclLiteral(e.p, hc)
                w = openclLiteral(e.w, hc)
                s = openclLiteral(e.s, hc)

                # cdfDoubleLogisticFunction = string("float cdfDoubleLogistic(float x, float p, float w, float s);\n",
                #                           "float cdfDoubleLogistic(float x, float p, float w, float s){\n",
                #                           "     float e1 = (-x+p-w)/s;\n",
                #                           "     float e2 = (-x+p+w)/s;\n",
                #                           "     if (e1>88.0f) {\n",
                #                           "        return 0.5;\n",
                #                           "     } else {\n",
                #                           "         return 0.5f*s*log((float)( (1.0f+exp((float)(e1)) )/(1.0f+exp((float)(e2)) ) ))/w - 0.5f;\n",
                #                           "     }\n",
                #                           "}\n")

                cdfDoubleLogisticFunctionDeclaration = string("float cdfDoubleLogistic(float x, float p, float w, float s);\n")
                cdfDoubleLogisticFunction = string(
                                        "float cdfDoubleLogistic(float x, float p, float w, float s) {
                                                                            float e1 = (-x+p-w)/s;
                                                                            float e2 = (-x+p+w)/s;
                                                                            if (e1<10.0f) {
                                                                                    return 0.5*s*log( (float) (
                                                                                                (1.0+exp((float)(e1))) /( 1.0+exp((float)(e2)))
                                                                                                             )
                                                                                                          )/w -0.5;
                                                                            }
                                                                            if (e1>=10.0f && e2<-10.0f ) {
                                                                                    return 0.5*(x-p)/w;
                                                                            }
                                                                            if (e2>=-10.0f && e2<10.0f) {
                                                                                    return 0.5*s*log(   (float)(  1.0 /(1.0+exp((float)(e2)))  )    )/w  + 0.5*(x-p)/w;
                                                                            }
                                                                            if (e2>=10.0f)
                                                                                    {return 0.5;}
                                                                            return 0.0f;
                                                                    }\n")
                cdfDoubleLogisticFunctionHash = string("cdfDoubleLogistic")
                hc.exprList[cdfDoubleLogisticFunctionHash] = (-3, cdfDoubleLogisticFunction, "function")
                hc.exprList[cdfDoubleLogisticFunctionHash].count = 0
                cdfDoubleLogisticFunctionDeclarationHash = string("cdfDoubleLogisticdeclaration")
                hc.exprList[cdfDoubleLogisticFunctionDeclarationHash] = (-137, cdfDoubleLogisticFunctionDeclaration,"function declaration")
                hc.exprList[cdfDoubleLogisticFunctionDeclarationHash].count = 0


                myExpression = string("cdfDoubleLogistic(",pos.expression,",",p.expression,",",w.expression,",",s.expression,")")

                if hc.exprList[e.hash].count > 1
                        eval = hc.exprList[e.hash].eval
                        defs = string("float ",eval, "=", myExpression, ";\n" )
                else
                        defs = ""
                        eval = string(myExpression )
                end
                defs = string(pos.definitions, p.definitions, w.definitions, s.definitions, defs)
                hc.exprList[e.hash].count = -hc.exprList[e.hash].count
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(defs, eval, dataAttributes(hc) )
end

function openclLiteral(e::logDoubleLogisticExpr, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(e,hc)
        end

        if hc.exprList[e.hash].count > 0
                pos = openclLiteral(e.pos, hc)
                p = openclLiteral(e.p, hc)
                w = openclLiteral(e.w, hc)
                s = openclLiteral(e.s, hc)

                logDoubleLogisticFunctionDeclaration = string("float logDoubleLogistic(float x, float p, float w, float s);\n")
                logDoubleLogisticFunction = string(
                                          "float logDoubleLogistic(float x, float p, float w, float s){\n",
                                          "     float e1 = (-x+p-w)/s;\n",
                                          "     float e2 = (-x+p+w)/s;\n",
                                          "     float c = 20.0f;\n",
                                          "       if (x-p<0.0f){\n",
                                          "           if (e1>c){\n",
                                          "               return   -e1-log((float)(w)) + \n",
                                          "               log(0.5f)+log(1.0f-exp((float)(-2.0f*w/s)) ) - exp((float)(-e1)) -exp((float)(-e2));\n",
                                           #+ c + log( 0.5* (1/(1+exp(c)) -  1/(1+exp(c+2w/s) ) ) )\n",
                                          "           } else {\n",
                                          "               return log( (float)(0.5f * (1.0/(1.0+exp((float)(e1))) - 1.0f/(1.0f+exp((float)(e2))) )/w));\n",
                                          "           }\n",
                                          "        } else {\n",
                                          "           if (e2<-c) {\n",
                                          "               return e2-log((float)(w)) + log(0.5f) + log(1.0f - exp((float)(-2.0f*w/s) ) )   - exp((float)(e2)) -exp((float)(e1));\n",
                                          "           } else {\n",
                                          "               return log(0.5f * (-1.0f/(1.0f+exp((float)(-e1))) +  1.0f/(1.0f+exp((float)(-e2))) )/w);\n",
                                          "           }\n",
                                          "       }\n",
                                          "}\n")
                logDoubleLogisticFunctionHash = string("logDoubleLogistic")
                hc.exprList[logDoubleLogisticFunctionHash] = (-4, logDoubleLogisticFunction, "function")
                hc.exprList[logDoubleLogisticFunctionHash].count = 0

                logDoubleLogisticFunctionDeclarationHash = string("logDoubleLogisticdeclaration")
                hc.exprList[logDoubleLogisticFunctionDeclarationHash] = (-137, logDoubleLogisticFunctionDeclaration,"function declaration")
                hc.exprList[logDoubleLogisticFunctionDeclarationHash].count = 0


                myExpression = string("logDoubleLogistic(",pos.expression,",",p.expression,",",w.expression,",",s.expression,")")

                if hc.exprList[e.hash].count > 1
                        eval = hc.exprList[e.hash].eval
                        defs = string("float ",eval, "=", myExpression, ";\n" )
                else
                        defs = ""
                        eval = string(myExpression )
                end
                defs = string(pos.definitions, p.definitions, w.definitions, s.definitions, defs)
                hc.exprList[e.hash].count = -hc.exprList[e.hash].count
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(defs, eval, dataAttributes(hc) )
end

function openclLiteral(dv::Vector{T}, hc::hiddenCounter = hiddenCounter()) where T<:ScalarDataIndex


        if hc.counter == 0
                firstPass(dv,hc)

                if length(dv) >0
                         hc.exprList["vectorOutput"] = (length(dv), "vectorSize","output")
                end
        end

        if length(dv) == 1
                return openclLiteral(dv[1], hc)
        else
                oclVec = [openclLiteral(d,hc) for d in dv]
                defs = string([ l.definitions for l in oclVec ]...)

                eval = string( #"(float[",length(dv),"])",
                        "{",[ string(l.expression,",") for l in oclVec[1:end-1]]..., oclVec[end].expression ,"}" )
                return literal(defs, eval, dataAttributes(hc))
        end
end

function openclVectorFunction(seed::ScalarDataIndex, commonVariables::X, individualVariables::X, nRepeat::String = "nRepeat")
        hc = hiddenCounter()
        firstPass(seed, hc)
        groupingCommands = findSimilarities(hc.prodList)
        secondPass(seed, groupingCommands)
        hc = hiddenCounter()
        firstPass(seed,hc)

        e = openclLiteral(seed,hc)
        hc.exprList[nRepeat] = (-1, nRepeat, "repeat length", "int")
        hc.exprList[nRepeat].extraAttributes = ["extraInput", "int"]

        localMem = false

        nParam =individualVariables.counter #maximum([i[3] for i in filter(x->x[2]==individualVariables.name,e.input)])

        myFunctionName = string(seed.name, "VectorCall")

        myFunctionDeclaration = string( "float ", myFunctionName, "(__global float * ",commonVariables.name,", __global float * ",individualVariables.name,");\n")
        myFunction = string(
                                "float ", myFunctionName, "(__global float * ",commonVariables.name,", __global float * ",individualVariables.name,"){\n",
                                e.definitions,
                                "return ",e.expression,";\n",
                                "}\n")


        vectorFunctionHash = myFunctionName
        hc.exprList[vectorFunctionHash] = (-8,string(myFunction), "function")
        hc.exprList[vectorFunctionHash].count = 0
        myFunctionDeclarationHash = string(myFunctionName, "declaration")
        hc.exprList[myFunctionDeclarationHash] = (-137,myFunctionDeclaration, "function declaration")
        hc.exprList[myFunctionDeclarationHash].count = 0


        globalindicator = "__global "
        outPutOffset = "gid*nRepeat + "
        if localMem
                globalindicator = ""
                outPutOffset = ""
        end

        objName = string(seed.name,"Vec")
        defs = string("for(int i=0; i!= nRepeat; i++){
                        //kernel should assign gid = get_global_id(0);
                        ",objName,"[",outPutOffset," i] = ",myFunctionName,"(",commonVariables.name,",",individualVariables.name,"+i*",nParam,");
                }\n")

        eval = objName

        outPutHash = string(objName, "outputMemoryPointer")

        if !haskey(hc.exprList, outPutHash)
                hc.exprList[outPutHash] = (-1, objName, string("output memory"), "float *")
                hc.exprList[outPutHash].extraAttributes = ["extraOutput","float *"]
        end


        literal(defs, eval, dataAttributes(hc))
end

function openclVectorFunction(ve::Vector{derivedParametericAvg})

        nVecLength = length(ve)

        vhc = hiddenCounter()
        eval = ""
        defs = ""

        localMem = false
        globalindicator = "__global "
        if localMem
                globalindicator = ""
        end

        for (eIndex, eref) in enumerate(ve)
                e = deRef(eref)
                myFunctionName = string("vectorDerivedParamAverageInput",eIndex)
                myKernelName = string("vectorDerivedParamAvg",eIndex)
                nDim = e.commonVariables.counter
                nParam = e.individualVariables.counter

                hc = hiddenCounter()

                logWeights = false
                logWeightsCall = ""
                logWeightFunctionName = ""

                reorder(e)

                if typeof(e.wSeed) == prodExpr

                        logWeights = true
                        logWSeedExpr = sumExpr( simplifyExpr.(logExpr.( deRef.(e.wSeed.px) ) ) )
                        myhc = hiddenCounter()
                        reorder(logWSeedExpr)
                        firstPass(logWSeedExpr, myhc)
                        groupingCommandsLogW = findSimilarities(myhc.prodList)
                        println("logW secondPass")
                        secondPass(logWSeedExpr, groupingCommandsLogW)
                        myhc = hiddenCounter()
                        firstPass(logWSeedExpr, myhc)
                        logWSeed = openclLiteral(logWSeedExpr, myhc)

                        logWeightFunctionName = string("vectorDerivedParamAverageLogWeight", eIndex)
                        logWeightFunctionDeclaration = string("float ",logWeightFunctionName,"(__global float *",e.commonVariables.name,",",
                                                                "__global float *",e.individualVariables.name,");\n")

                        logWeightFunction = string("float ",logWeightFunctionName,"(__global float *",e.commonVariables.name,",",
                                                                "__global float *",e.individualVariables.name,"){\n",
                                                                logWSeed.definitions,
                                                                "return ",logWSeed.expression,";\n",
                                                                "}\n")
                        logWeightFunctionHash = logWeightFunctionName
                        logWeightFunctionDeclarationHash = string(logWeightFunctionName, "declaration")
                        myhc.exprList[logWeightFunctionDeclarationHash] = (-137, logWeightFunctionDeclaration, "function declaration")
                        myhc.exprList[logWeightFunctionDeclarationHash].count = 0

                        myhc.exprList[logWeightFunctionHash] = (-8, logWeightFunction,"function")
                        myhc.exprList[logWeightFunctionHash].count


                        for k in keys(myhc.exprList)
                                vhc.exprList[k] = myhc.exprList[k]
                        end
                        println("end of logW secondPass")
                end
                firstPass(e,hc)
                groupingCommands = findSimilarities(hc.prodList)
                secondPass(e, groupingCommands)

                hc = hiddenCounter()
                firstPass(e,hc)

#                hc.exprList[e.hash].eval = string(e.name,"Vec",eIndex)

                xSeed = openclLiteral(e.xSeed, hc)
                wSeed = openclLiteral(e.wSeed, hc)

                xSeedDerived = ( [openclLiteral( f,hc ) for f in e.xSeedDerived[1]], [openclLiteral( f,hc ) for f in e.xSeedDerived[2]] )
                wSeedDerived = ( [openclLiteral( f,hc ) for f in e.wSeedDerived[1]], [openclLiteral( f,hc ) for f in e.wSeedDerived[2]] )

                nRepeat = e.nRepeat

                e.nRepeat.hash = string(e.nRepeat)
                if !haskey(hc.exprList, e.nRepeat.hash)
                        hc.exprList[e.nRepeat.hash] = (-1, e.nRepeat.name, string(e.nRepeat.attribute), "int")
                        hc.exprList[e.nRepeat.hash].extraAttributes = ["extraInput","int"]
                end
                workMemHash = "workMemoryPointer"
                if !haskey(hc.exprList, workMemHash)
                        hc.exprList[workMemHash] = (-1, "workmemory", string("work memory"), "float *")
                        hc.exprList[workMemHash].extraAttributes = ["extraInput","float *"]
                end
                objName = string(e.name,"Obj")
                outPutHash = "outputMemoryPointer"
                if !haskey(hc.exprList, outPutHash)
                        hc.exprList[outPutHash] = (-1, objName, string("output memory"), "float *")
                        hc.exprList[outPutHash].extraAttributes = ["extraOutput","float *"]
                end


                # nDim = length(filter(x->x[2]== e.commonVariables.name, wSeed.input))
                # nParam =length(filter(x->x[2]==e.individualVariables.name, wSeed.input))


                #workMemoryOffset = string("+ gid*(2*",nDim,"+2*",nParam," + 2)*nRepeat*",nVecLength)
                workMemoryOffset = string("+ gid*(2*",nDim,"+2*",nParam," + 2)*nRepeat") #reusing work memory, when execution is sequential
                outputOffset = string("+ gid*",nVecLength,"*(",nDim,"+",nParam,"*nRepeat)")
                if localMem
                        workMemoryOffset = ""
                        outputOffset = ""
                end


                xIMin = -1
                xIMax = -1
                wIMin = -1
                wIMax = -1

                for i in 1:nParam
                        if !(typeof(e.xSeedDerived[2][i]) == constIndex && e.xSeedDerived[2][i].value == 0.0)
                                if xIMin == -1
                                        xIMin = i-1
                                end
                                xIMax = i-1
                        end

                        if !(typeof(e.wSeedDerived[2][i]) == constIndex && e.wSeedDerived[2][i].value == 0.0)
                                if wIMin == -1
                                        wIMin = i-1
                                end
                                wIMax = i-1
                        end
                end

                xCMin = -1
                xCMax = -1
                wCMin = -1
                wCMax = -1

                for i in 1:nDim
                        if !(typeof(e.xSeedDerived[1][i]) == constIndex && e.xSeedDerived[1][i].value == 0.0)
                                if xCMin == -1
                                        xCMin = i-1
                                end
                                xCMax = i-1
                        end

                        if !(typeof(e.wSeedDerived[1][i]) == constIndex && e.wSeedDerived[1][i].value == 0.0)
                                if wCMin == -1
                                        wCMin = i-1
                                end
                                wCMax = i-1
                        end
                end

                if logWeights
                        logWeightsCall = string("
                                float maxLogW = -INFINITY;
                                float penultimateLogW = -INFINITY;
                                int maxWindex = -1;
                                for(int i=0; i<nRepeat; i++){
                                                __global float * myParams = &",e.individualVariables.name,"[i*",nParam,"];
                                                float logW = ",logWeightFunctionName,"(ourCoord,myParams);
                                                if(logW > maxLogW)
                                                {
                                                        penultimateLogW = maxLogW;
                                                        maxLogW = logW; maxWindex = i;

                                                }
                                                else {
                                                        if(logW > penultimateLogW)
                                                                {penultimateLogW = logW;}
                                                }
                                        }
                                if(maxLogW<-25.0f || maxLogW-penultimateLogW > 15.0f) //-15.0f and 7.0f
                                {
                                        //call myfunction that returns derivatives, only for maxWindex;
                                        //
                                        __global float * myParams = &",e.individualVariables.name,"[maxWindex*",nParam,"];
                                        ",myFunctionName,"(ourCoord,myParams,",
                                        "xSeed+maxWindex,wSeed+maxWindex,",
                                        "xSeedDerivedC+maxWindex*",nDim,", xSeedDerivedP + maxWindex*",nParam,",",
                                        "wSeedDerivedC+maxWindex*",nDim,", wSeedDerivedP + maxWindex*",nParam,");\n",

                                        "for(int j=",xIMin,";j<=",xIMax,";j++) {
                                                derivedP[ (maxWindex*",nParam,"+j)*",nVecLength,"+",eIndex-1," ] = xSeedDerivedP[j+",nParam,"*maxWindex];
                                        }
                                        for(int j=",xCMin,";j<=",xCMax,";j++){
                                                derivedC[(j*",nVecLength,"+",eIndex-1,")] += xSeedDerivedC[j+maxWindex*",nDim,"];
                                        }
                                        return;
                                }
                                ")
                end

                derivedParametericAvgFunctionDeclaration = string("void ",myFunctionName,"(__global float * ",e.commonVariables.name,",",
                                                                "__global float * ",e.individualVariables.name,",",
                                                                "",globalindicator," float * xSeed, ",globalindicator," float * wSeed,",
                                                                "",globalindicator," float *xSeedDiffC, ",globalindicator," float * xSeedDiffP,",
                                                                "",globalindicator," float *wSeedDiffC, ",globalindicator," float * wSeedDiffP",
                                                                ");\n")
                derivedParametericAvgFunction = string(
                                                       "void ",myFunctionName,"(__global float * ",e.commonVariables.name,",",
                                                                "__global float *",e.individualVariables.name,",",
                                                                "",globalindicator," float * xSeed, ",globalindicator," float * wSeed,",
                                                                "",globalindicator," float *xSeedDiffC, ",globalindicator," float * xSeedDiffP,",
                                                                "",globalindicator," float *wSeedDiffC, ",globalindicator," float * wSeedDiffP",
                                                                ") {\n",
                                                                xSeed.definitions,
                                                                wSeed.definitions,
                                                                [xS.definitions for xS in xSeedDerived[1]]...,
                                                                [xS.definitions for xS in xSeedDerived[2]]...,
                                                                [wS.definitions for wS in wSeedDerived[1]]...,
                                                                [wS.definitions for wS in wSeedDerived[2]]...,

                                                                "\n",
                                                                "//****return values***//\n",
                                                                "xSeed[0] = ", xSeed.expression,";\n",
                                                                "wSeed[0] = ", wSeed.expression,";\n",
                                                                [ string("xSeedDiffC[",index-1,"] = " ,xSeedDerived[1][index].expression,";\n"  ) for index in 1:nDim ]...,
                                                                [ string("xSeedDiffP[",index-1,"] = " ,xSeedDerived[2][index].expression,";\n"  ) for index in 1:nParam ]...,
                                                                [ string("wSeedDiffC[",index-1,"] = " ,wSeedDerived[1][index].expression,";\n"  ) for index in 1:nDim ]...,
                                                                [ string("wSeedDiffP[",index-1,"] = " ,wSeedDerived[2][index].expression,";\n"  ) for index in 1:nParam ]...,
                                                                "}\n")




                derivedParAvgFunc2Declaration = string("void ",myKernelName,"(__global float * ",e.commonVariables.name,",",
                                            "__global float * ",e.individualVariables.name,", int nRepeat, ",globalindicator," float *workmemory, ",globalindicator," float * output);\n")
                derivedParAvgFunc2 = string( "void ",myKernelName,"(__global float * ",e.commonVariables.name,",",
                                            "__global float * ",e.individualVariables.name,", int nRepeat, ",globalindicator," float *workmemory, ",globalindicator," float * output){
                                                const int gid = get_global_id(0);

                                                //__global float temp[(2*",nDim,"+2*",nParam," + 2)*nRepeat];
                                                ",globalindicator," float *xSeed = workmemory ",workMemoryOffset,";//+2*(",nDim+nParam+1,")*",eIndex-1,"*nRepeat;
                                                ",globalindicator," float *wSeed = xSeed+nRepeat;
                                                ",globalindicator," float *xSeedDerivedC = wSeed+nRepeat;
                                                ",globalindicator," float *xSeedDerivedP = xSeedDerivedC+",nDim,"*nRepeat;
                                                ",globalindicator," float *wSeedDerivedC = xSeedDerivedP+",nParam,"*nRepeat;
                                                ",globalindicator," float *wSeedDerivedP = wSeedDerivedC+",nDim,"*nRepeat;

                                                ",globalindicator," float *derivedC = output ",outputOffset,";//+",nDim,"*",eIndex-1,";
                                                ",globalindicator," float *derivedP = output ",outputOffset,"+",nDim,"*",nVecLength,";// + ",nParam,"*",eIndex-1,"*nRepeat;
                                                for(int i=0;i<",nDim,";i++)
                                                        {derivedC[(i*",nVecLength,"+",eIndex-1,")] = 0.0f;}
                                                        //{derivedC[i] = 0.0f;}
                                                for(int i=0;i<",nParam,"*nRepeat;i++)
                                                        {derivedP[(i*",nVecLength,"+",eIndex-1,")] = 0.0f;}


                                                //__global float *ourCoord = &",e.commonVariables.name,"[gid*",nDim,"];
                                                __global float *ourCoord = ",e.commonVariables.name,";
                                                float sumW = 0.0f;
                                                float avg = 0.0f;
                                                //float maxAbsW = 0.0f;
                                                //int maxWindex = -1;
                                                ",logWeightsCall,"

                                                for(int i=0; i<nRepeat; i++) {

                                                                __global float * myParams = &",e.individualVariables.name,"[i*",nParam,"];
                                                                ",myFunctionName,"(ourCoord,myParams,",
                                                                "xSeed+i,wSeed+i,",
                                                                "xSeedDerivedC+i*",nDim,", xSeedDerivedP + i*",nParam,",",
                                                                "wSeedDerivedC+i*",nDim,", wSeedDerivedP + i*",nParam,");\n",
                                                                #"if(fabs((float)wSeed[i])>maxAbsW){maxAbsW = fabs((float)wSeed[i]); maxWindex = i; }",
                                                                "sumW += wSeed[i];\n",
                                                                "avg += wSeed[i]*xSeed[i];\n",
                                                "}\n",

                                                # "if(maxWindex!=0 && fabs((float)(maxAbsW-sumW))<1e-6f )\n",
                                                # "{
                                                #         for(int j=",xIMin,";j<=",xIMax,";j++) {
                                                #                 derivedP[ (maxWindex*",nParam,"+j)*",nVecLength,"+",eIndex-1," ] = xSeedDerivedP[j+",nParam,"*maxWindex];
                                                #         }
                                                #         for(int j=",xCMin,";j<=",xCMax,";j++){
                                                #                 derivedC[(j*",nVecLength,"+",eIndex-1,")] += xSeedDerivedC[j+maxWindex*",nDim,"];
                                                #         }
                                                #
                                                #         return;
                                                # }",


                                                "avg = avg/sumW;\n",
                                                "for(int i=0; i<nRepeat; i++){\n",
                                                                "for(int j=",xIMin,";j<=",xIMax,";j++) {
                                                                        derivedP[(i*",nParam,"+j)*",nVecLength,"+",eIndex-1,"] = wSeed[i]*xSeedDerivedP[j+",nParam,"*i]/sumW;
                                                                }\n",

                                                                "for(int j=",wIMin,";j<=",wIMax,";j++) {
                                                                        /*the indices of non-zero derivatives of w and x
                                                                                don't overlap, but it is important to use +=,
                                                                                instead of =, so the zero valued derivatives
                                                                                don't overwrite the values given in the line above */
                                                                        derivedP[(i*",nParam,"+j)*",nVecLength,"+",eIndex-1,"] += wSeedDerivedP[i*",nParam,"+j]*(xSeed[i] - avg)/sumW;
                                                                }\n",
                                                                "for(int j=",xCMin,";j<=",xCMax,";j++){
                                                                        derivedC[(j*",nVecLength,"+",eIndex-1,")] += wSeed[i]*xSeedDerivedC[j+i*",nDim,"];
                                                                }\n",
                                                                "for(int j=",wCMin,";j<=",wCMax,";j++){
                                                                        derivedC[(j*",nVecLength,"+",eIndex-1,")] += wSeedDerivedC[j+i*",nDim,"]*(xSeed[i] - avg ) ;
                                                                        //
                                                                }\n",
                                                        "}\n",
                                                "\n",
                                                "for(int i=0; i<",nDim,"; i++ ) {derivedC[i*",nVecLength,"+",eIndex-1,"]/=sumW;}\n",

                                            "}\n",


                                            )
                                            #for(int j=0;j<=",nParam-1,";j++) {
                                            #"


                    # derivedParametericAvgFunctionHash = myFunctionName
                    # hc.exprList[derivedParametericAvgFunctionHash] = (-7,derivedParametericAvgFunction, "function")
                    # hc.exprList[derivedParametericAvgFunctionHash].count = 0

                    derivedParametericAvgKernelHash = myKernelName
                    hc.exprList[derivedParametericAvgKernelHash] = (-8,string(derivedParametericAvgFunction,derivedParAvgFunc2), "function")
                    hc.exprList[derivedParametericAvgKernelHash].count = 0
                    #hc.exprList[e.hash].count = -hc.exprList[e.hash].count

                    derivedParAvgFunc2DeclarationHash = string(myKernelName,"declaration")
                    hc.exprList[derivedParAvgFunc2DeclarationHash] = (-137,string(derivedParametericAvgFunctionDeclaration,derivedParAvgFunc2Declaration), "function declaration")
                    hc.exprList[derivedParAvgFunc2DeclarationHash].count = 0

                    for k in keys(hc.exprList)
                            vhc.exprList[k] = hc.exprList[k]
                    end

                    defs = string(defs, myKernelName,"(",e.commonVariables.name,",",e.individualVariables.name,",",nRepeat.name,",workmemory,",objName,");\n")

        end

        if false
                #should change workmemory offset, if execution is not sequenctial, using separate kernels
                eval = string( "__kernel void diffAvg(__global float * ",ve[1].commonVariables.name,",",
                               "__global float * ",ve[1].individualVariables.name,", int nRepeat, ",globalindicator," float *workmemory, ",globalindicator," float * ",objName,"){\n",
                               eval,
                               "}\n")
        end
        literal(defs, eval, dataAttributes(vhc) )

end

function exprFactMap(prodList::Dict{String,Vector{String}})
        factorCountMap = Dict{String, Int}()
        factorCountMapMap = Dict{String, Dict{String,Int}}()
        for (k,v) in prodList
                vcm = countmap( v )
                factorCountMapMap[k] = vcm
                for c in vcm
                        if !haskey(factorCountMap,c.first)
                                push!(factorCountMap, c)
                        else
                                factorCountMap[c.first] = max(c.second, factorCountMap[c.first])
                        end
                end
        end


        factorCumulativeSum = cumsum(collect(values(factorCountMap)))
        if length(factorCumulativeSum)!=0
                nKeys = factorCumulativeSum[end]
        else
                nKeys = 0
        end
        factorReverse = Vector{String}(undef, nKeys)
        counter = 1
        for (k,v) in factorCountMap
                for i in 1:v
                        factorReverse[counter] = k
                        counter += 1
                end

        end
        factorCumulativeSum = vcat(0,factorCumulativeSum[1:end-1])
        facturMapCumulativeSum = Dict( keys(factorCountMap) .=> factorCumulativeSum )
        expressionFactorVector = Dict{String, Vector{Int}}()
        for k1 in keys(factorCountMapMap)
                cvec = zeros(nKeys)
                for k2 in keys(factorCountMapMap[k1])
                        #println(  facturMapCumulativeSum[k2], " ", k2, " ",  factorCountMapMap[k1][k2])
                        for i in 1:factorCountMapMap[k1][k2]
                                cvec[facturMapCumulativeSum[k2]+i ] = 1
                        end
                end
                expressionFactorVector[k1] = cvec
        end
        return expressionFactorVector, factorReverse, factorCountMapMap, factorCumulativeSum
end


function removeElementO( exprList::Array, overlapArray::Array, coord::CartesianIndex )
        overlap = exprList[coord[1],:] .* exprList[coord[2],:]
        NexprList = deepcopy(exprList)
        NexprList[ coord[1],:] .-= overlap
        NexprList[ coord[2],:] .-= overlap
        newExpr = zeros(size(exprList)[1]+1)
        newExpr[coord[1]] = 1
        newExpr[coord[2]] = 1
        NNexprList = hcat(vcat(NexprList, overlap') , newExpr)

        NoverlapArray = hcat( vcat(overlapArray, zeros(size(overlapArray)[1])' ), zeros(size(overlapArray)[1]+1) )

        #max is at coord[1],coord[2]
         #i=1:N,  j=1:i-1
         #NNoverlapArray[i,j] =dot(  h[ i,: ],h[j,:] )
        #rewrite every coord[1] -fixed, and coord[2] fixed
        # when i = coord[1]
        # j=1:coord[1]-1
        for j=1:coord[1]-1
                NoverlapArray[ coord[1],j ] = dot(NNexprList[coord[1],:], NNexprList[j,:])
        end

        #when j=coord[2]
        # j<i-1 -> i= j+1:end
        for i=coord[2]+1:size(NNexprList)[1]
                NoverlapArray[i, coord[2] ] = dot(NNexprList[i,:], NNexprList[coord[2],:])
        end


        return NNexprList, NoverlapArray
end

function removeElement( exprList::Array, overlapArray::Array, hashVector::Vector, overlapIndices::Vector{Tuple{Vector{Int},CartesianIndex{2}}}, coord::CartesianIndex , noBreakStr = "£prod")
        overlap = exprList[coord[1],:] .* exprList[coord[2],:]
        exprList = deepcopy(exprList)
        exprList[ coord[1],:] .-= overlap
        exprList[ coord[2],:] .-= overlap
        # println(coord)
        # println(overlap)
        push!(overlapIndices,(findall(x->x==1, overlap),coord) )
        newExpr = zeros(size(exprList)[1]+1)
        newExpr[coord[1]] = 1
        newExpr[coord[2]] = 1
        exprList = hcat(vcat(exprList, overlap') , newExpr)
        overlapArray = hcat( vcat(overlapArray, zeros(size(overlapArray)[1])' ), zeros(size(overlapArray)[1]+1) )


        #max is at coord[1],coord[2]
         #i=1:N,  j=1:i-1
         #overlapArray[i,j] =dot(  h[ i,: ],h[j,:] )
        #rewrite every coord[1] -fixed, and coord[2] fixed
        # when i = coord[1]
        # j=1:coord[1]-1
        for j=1:coord[1]-1
                overlapArray[ coord[1],j ] = dot(exprList[coord[1],:], exprList[j,:])
        end
        for j=1:coord[2]-1
                overlapArray[ coord[2],j ] = dot(exprList[coord[2],:], exprList[j,:])
        end

        #when j=coord[2]
        # j<i-1 -> i= j+1:end
        for i=coord[2]+1:size(exprList)[1]
                overlapArray[i, coord[2] ] = dot(exprList[i,:], exprList[coord[2],:])

        end
        for i=coord[1]+1:size(exprList)[1]
                overlapArray[i, coord[1] ] = dot(exprList[i,:], exprList[coord[1],:])

        end
        nExprList = size(exprList)[1]
        for j=1:nExprList-1
                overlapArray[nExprList,j] = dot(exprList[end,:], exprList[j,:])
        end

        push!(hashVector, string(noBreakStr,"(",hashVector[coord[1]],",",hashVector[coord[2]],")") )
        return exprList, overlapArray
end

function applyRemove(pe::prodExpr, removeCommands::Dict{String, Vector{ Tuple{String,String} } } )
        if haskey(removeCommands,pe.hash)
                # @show removeCommands[pe.hash]
                # @show pe.hash
                println("=> applyRemove \n")
                #println(" id is ",pe.id)
                #println(removeCommands[pe.hash])
                #pe.px = deepcopy(pe.px)

                for block in removeCommands[pe.hash]
                        #println(pe.hash)
                        #println(block)
                        for pex in pe.px
                                println("  -  ",pex.hash)
                        end
                        ff1 = findfirst( x->x.hash==block[1], pe.px )
                        if block[1] == block[2]
                                println(ff1)

                                ff2 = findnext(x->x.hash==block[2], pe.px,ff1+1)
                        else

                                ff2 = findfirst( x->x.hash==block[2], pe.px )
                        end
                        #println("ff = ", ff1, "    ", ff2)
                        if length(pe.px) > 2
                                blockProd = prodExpr([pe.px[ff1], pe.px[ff2]], nothing, true )
                                blockProd.hash = string("£prod(",block[1],",",block[2],")")
                                deleteat!(pe.px, sort([ff1, ff2]))
                                push!(pe.px, blockProd)
                        else
                                pe.keepAsBlock = true
                                pe.hash = string("£prod(",block[1],",",block[2],")")
                        end
                end
                #pe.id = rand(1:999)
                #println(" end id is ", pe.id)
                pe.hash = string(pe)
                println("<= applyRemove end")
        end
end

function findSimilarities( exprList::Dict{String,Vector{String}}, noBreakHash = "£prod" )
        b,bfm, bfmm, bfcs = exprFactMap(exprList)

        h = hcat(values(b)...)
        d=zeros(size(h)[1],size(h)[1])
        for i=1:size(h)[1]
                for j=1:i-1
                        d[i,j] = dot(h[i,:], h[j,:])
                end
        end

        overlapIndices = Vector{Tuple{Vector{Int},CartesianIndex{2}}}(undef,0)
        println("overlapIndices = ",overlapIndices)
        removeCommands = Dict{String, Vector{ Tuple{String,String} } }()

        if length(d) == 0
                return removeCommands
        end
        m = findmax(d)

        while m[1] > 1
                h, d = removeElement(h, d, bfm, overlapIndices, m[2], noBreakHash)

                m = findmax(d)
                #println(m, " " ,size(h), " ", size(d))
        end

        println("overlapIndices = ",overlapIndices)

        exprKeys = collect(keys(b))
        #overlapHashes = Vector{Tuple{Vector{String}, Tuple{String,String}}}(undef,0)
        for oi in overlapIndices
                for o in oi[1]
                        if !haskey(removeCommands,exprKeys[o])
                                removeCommands[exprKeys[o]] = Tuple{String,String}[]
                        end
                        push!( removeCommands[exprKeys[o]], (bfm[oi[2][1]], bfm[oi[2][2]] ) )
                end
        end

        return removeCommands
end
