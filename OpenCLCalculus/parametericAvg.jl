mutable struct parametericAvg
        id::Union{Int,Nothing}
        name::String

        nRepeat::ScalarDataIndex
        xSeed::ScalarDataIndex
        wSeed::ScalarDataIndex

        commonVariables::X
        individualVariables::X

        hash::Union{Nothing, String}

        function parametericAvg(xSeed::ScalarDataIndex, wSeed::ScalarDataIndex, nRepeat::ScalarDataIndex, commonVariables::X, individualVariables::X, id::Union{Int,Nothing} = nothing)
                new(id, "parametericAvg", nRepeat, xSeed, wSeed, commonVariables, individualVariables, nothing )
        end
end

function Base.:string(pa::parametericAvg)
        if pa.id == nothing
                return string("parametericAvg(",string(pa.nRepeat),",",string(pa.xSeed),",",string(pa.wSeed),")")
        else
                return string("parametericAvg",pa.id)
        end
end

function juliaLiteral(e::parametericAvg, hc::hiddenCounter = hiddenCounter())
        if e.id == nothing
                hc.counter +=1
                baseName =string(e.name,"Anon", hc.counter)
        else
                baseName = string(e.name,e.eid)
        end

        xSeed = juliaLiteral(e.xSeed, hc)
        hc2 = hiddenCounter()
        wSeed = juliaLiteral(e.wSeed, hc2)
        #nRepeat = juliaLiteral(e.nRepeat, hc)
        nRepeat = e.nRepeat
        e.nRepeat.hash = string(e.nRepeat.name)
        if !haskey(hc.exprList, e.nRepeat.hash)
                hc.exprList[e.nRepeat.hash] = (-1, e.nRepeat.name, string(e.nRepeat.attribute), "Int")
                hc.exprList[e.nRepeat.hash].extraAttributes = ["extraInput","Int"]
        end



        nDim = e.commonVariables.counter #maximum([i[3] for i in filter(x->x[2]== e.commonVariables.name, w.input)]) # w is the latest juliaLiteral, with updated hiddenCounter
        nParam =e.individualVariables.counter# maximum([i[3] for i in filter(x->x[2]== e.individualVariables.name, w.input)])

        e.hash = string(e)
        if !haskey(hc.exprList, e.hash)
                eval = baseName

                myFunctionNameX = string("calculateParametericAverageInputX",hc.counter)
                parametericAvgFunctionX = string("function ",myFunctionNameX,"( ",e.commonVariables.name,"::Vector{Float64},",
                                        e.individualVariables.name,"::Vector{Float64} )\n",
                                        xSeed.definitions,
                                        "return ",xSeed.expression,"\n",
                                        "end\n")
                myFunctionNameW = string("calculateParametericAverageInputW",hc.counter)
                parametericAvgFunctionW = string("function ",myFunctionNameW,"( ",e.commonVariables.name,"::Vector{Float64},",
                                        e.individualVariables.name,"::Vector{Float64} )\n",
                                        wSeed.definitions,
                                        "return ",wSeed.expression,"\n",
                                        "end\n")

                parametericAvgFunctionXHash = string("function ", myFunctionNameX)
                parametericAvgFunctionWHash = string("function ", myFunctionNameW)

                myXName = string("paramAvgX",hc.counter)
                myWName = string("paramAvgW",hc.counter)

                defs = string( myXName, "= Vector{Float64}(undef, ", nRepeat.name ,")\n",
                                myWName, "= Vector{Float64}(undef, ", nRepeat.name ,")\n",
                                "for i in 1:", nRepeat.name,"\n",
                                "    ",myXName,"[i] = ",  myFunctionNameX, "(",e.commonVariables.name,"[1:",nDim,"],\n",
                                "          ",e.individualVariables.name,"[(i-1)*",nParam,"+1:i*",nParam,"]",")\n",
                                "    ",myWName,"[i] = ",  myFunctionNameW, "(",e.commonVariables.name,"[1:",nDim,"],\n",
                                "          ",e.individualVariables.name,"[(i-1)*",nParam,"+1:i*",nParam,"]",")\n",
                                "end\n",

                                baseName,"= dot(",myXName,",",myWName,")/sum(",myWName,")\n"
                 )


                if !haskey(hc.exprList, parametericAvgFunctionXHash)
                        hc.exprList[parametericAvgFunctionXHash] = (hc.counter, parametericAvgFunctionX, "function")
                        defs = string( parametericAvgFunctionX, defs  )
                end

                if !haskey(hc.exprList, parametericAvgFunctionWHash)
                        hc.exprList[parametericAvgFunctionWHash] = (hc.counter, parametericAvgFunctionW, "function")
                        defs = string( parametericAvgFunctionW, defs  )
                end

                for k in keys(hc2.exprList)
                        hc.exprList[k] = hc2.exprList[k]
                end
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end

        literal(defs, eval, dataAttributes(hc)) #updata data attributes nRepeat times
end

mutable struct derivedParametericAvg
        id::Union{Int,Nothing}
        name::String

        nRepeat::ScalarDataIndex #I don't pass this in members(), and won't scan with firstPass
        xSeed::ScalarDataIndex
        wSeed::ScalarDataIndex
        xSeedDerived::Tuple{Vector{ScalarDataIndex},Vector{ScalarDataIndex}  }
        wSeedDerived::Tuple{Vector{ScalarDataIndex},Vector{ScalarDataIndex}  }

        commonVariables::X
        individualVariables::X

        hash::Union{Nothing, String}

        function derivedParametericAvg(xSeed::ScalarDataIndex, wSeed::ScalarDataIndex, nRepeat::ScalarDataIndex, commonVariables::X, individualVariables::X, id::Union{Int,Nothing} = nothing)
                new(id, "derivedParametericAvg", nRepeat, xSeed, wSeed, derive(simplifyExpr(xSeed), commonVariables, individualVariables),derive(simplifyExpr(wSeed), commonVariables, individualVariables),  commonVariables, individualVariables, nothing )
        end
end

function reorder(dpa::derivedParametericAvg)
        reorder(dpa.xSeed)
        reorder(dpa.wSeed)
        reorder.(dpa.xSeedDerived[1])
        reorder.(dpa.xSeedDerived[2])

        reorder.(dpa.wSeedDerived[1])
        reorder.(dpa.wSeedDerived[2])

        dpa.hash = string(dpa)
        return nothing
end

function Base.:string(pa::derivedParametericAvg)
        if pa.id == nothing
                return string("derivedParametericAvg(",string(pa.nRepeat),",",string(pa.xSeed),",",string(pa.wSeed),")")
        else
                return string("derivedParametericAvg",pa.id)
        end
end

function members(e::derivedParametericAvg)
        return [e.xSeed, e.wSeed, e.xSeedDerived[1][:]..., e.xSeedDerived[2][:]..., e.wSeedDerived[1][:]..., e.wSeedDerived[2][:]...  ]
end

function juliaLiteral(e::derivedParametericAvg, hc::hiddenCounter = hiddenCounter())
        if e.id == nothing
                hc.counter +=1
                baseName =string(e.name,"Anon", hc.counter)
        else
                baseName = string(e.name,e.eid)
        end

        xSeed = juliaLiteral(e.xSeed,hc)
        wSeed = juliaLiteral(e.wSeed,hc)

        xSeedDerived = (juliaLiteral(e.xSeedDerived[1],hc), juliaLiteral(e.xSeedDerived[2],hc) )
        wSeedDerived = (juliaLiteral(e.wSeedDerived[1],hc), juliaLiteral(e.wSeedDerived[2],hc) )

        nRepeat = e.nRepeat

        e.nRepeat.hash = string(e.nRepeat)
        if !haskey(hc.exprList, e.nRepeat.hash)
                hc.exprList[e.nRepeat.hash] = (-1, e.nRepeat.name, string(e.nRepeat.attribute), "Int")
                hc.exprList[e.nRepeat.hash].extraAttributes = ["extraInput","Int"]
        end

        nDim = length(filter(x->x[2]== e.commonVariables.name, wSeed.input))
        nParam =length(filter(x->x[2]==e.individualVariables.name, wSeed.input))

        e.hash = string(e)

        if !haskey(hc.exprList, e.hash)
                eval = baseName
                myFunctionName = string("calculateDerivedParametericAverageInput",hc.counter)

                xIMin = -1
                xIMax = -1
                wIMin = -1
                wIMax = -1

                for i in 1:nParam
                        if !(typeof(e.xSeedDerived[2][i]) == constIndex && e.xSeedDerived[2][i].value == 0.0)
                                if xIMin == -1
                                        xIMin = i
                                end
                                xIMax = i
                        end

                        if !(typeof(e.wSeedDerived[2][i]) == constIndex && e.wSeedDerived[2][i].value == 0.0)
                                if wIMin == -1
                                        wIMin = i
                                end
                                wIMax = i
                        end
                end

                xCMin = -1
                xCMax = -1
                wCMin = -1
                wCMax = -1

                for i in 1:nDim
                        if !(typeof(e.xSeedDerived[1][i]) == constIndex && e.xSeedDerived[1][i].value == 0.0)
                                if xCMin == -1
                                        xCMin = i
                                end
                                xCMax = i
                        end

                        if !(typeof(e.wSeedDerived[1][i]) == constIndex && e.wSeedDerived[1][i].value == 0.0)
                                if wCMin == -1
                                        wCMin = i
                                end
                                wCMax = i
                        end
                end


                derivedParametericAvgFunction = string("function ",myFunctionName,"(",e.commonVariables.name,"::Vector{Float64},",
                                                                e.individualVariables.name,"::Vector{Float64} )\n",
                                                                xSeed.definitions,
                                                                wSeed.definitions,
                                                                xSeedDerived[1].definitions,
                                                                xSeedDerived[2].definitions,
                                                                wSeedDerived[1].definitions,
                                                                wSeedDerived[2].definitions,

                                                                "return ",
                                                                xSeed.expression,", ",
                                                                wSeed.expression,", ",
                                                                "(",xSeedDerived[1].expression,", ",xSeedDerived[2].expression,"), ",
                                                                "(",wSeedDerived[1].expression,", ",xSeedDerived[2].expression,")\n",
                                                                "end\n")

                derivedParametericAvgFunctionHash = string("function ", myFunctionName)

                myXName = string("diffParamAvgX", hc.counter)
                myWName = string("diffParamAvgW", hc.counter)
                myDiffXName = string("diffParamAvgDiffX", hc.counter)
                myDiffWName = string("diffParamAvgDiffW", hc.counter)

                sumWName = string("sumW", hc.counter)
                avgName = string("pAvg", hc.counter)

                defs = string(  eval ," = (zeros(",nDim,"),zeros(",nParam,"*",nRepeat.name,")  )  \n",#"= (Vector{Float64}(undef, ",nDim,"),Vector{Float64}(undef, ",nParam,"))\n"
                                myXName, "= Vector{Float64}(undef, ",nRepeat.name ,")\n",
                                myWName, "= Vector{Float64}(undef, ",nRepeat.name ,")\n",
                                myDiffXName, " = Vector{Tuple{ Vector{Float64}, Vector{Float64} }}(undef,",nRepeat.name,")\n",
                                myDiffWName, " = Vector{Tuple{ Vector{Float64}, Vector{Float64} }}(undef,",nRepeat.name,")\n",
                                "for i in 1:", nRepeat.name,"\n",
                                        myXName, "[i], ", myWName, "[i], ", myDiffXName, "[i], ", myDiffWName, "[i] = ",
                                                myFunctionName,"(",e.commonVariables.name,"[1:",nDim,"],\n",
                                                "          ",e.individualVariables.name,"[(i-1)*",nParam,"+1:i*",nParam,"]",")\n",

                                "end\n",

                                sumWName, " = sum(",myWName,")\n",
                                avgName, " = dot(",myWName,",",myXName,")/",sumWName,"\n",

                                "for j in 1:", nRepeat.name,"\n",
                                "       for i in ",xIMin,":",xIMax,"\n",
                                "               ",eval,"[2][i + (j-1)*",nParam,"] = ", myWName,"[j]*",myDiffXName,"[j][2][i]/",sumWName, "\n",
                                "       end\n",
                                "       for i in ",wIMin,":",wIMax,"\n",
                                "               ",eval,"[2][i + (j-1)*",nParam,"] = ", myDiffWName,"[j][2][i]*(",myXName,"[j] - ",avgName," )/",sumWName,"\n",
                                "       end\n",
                                "       for i in ",xCMin,":",xCMax,"\n",
                                "               ",eval,"[1][i] += ",myWName,"[j]*",myDiffXName,"[j][1][i]","\n",
                                "       end\n",
                                "       for i in ",wCMin,":",wCMax,"\n",
                                "               ",eval,"[1][i] += ",myDiffWName,"[j][1][i]*(",myXName,"[j] - ",avgName,")\n",
                                "       end\n",
                                "end\n",
                                eval,"[1]./=",sumWName,"\n"

                                )

                if !haskey(hc.exprList, derivedParametericAvgFunctionHash)
                        hc.exprList[derivedParametericAvgFunctionHash] = (hc.counter, derivedParametericAvgFunction, "function")
                        defs = string( derivedParametericAvgFunction, defs  )
                end
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end

        literal(defs, eval, dataAttributes(hc)) #updata data attributes nRepeat times
end

function openclLiteral(e::derivedParametericAvg, hc::hiddenCounter = hiddenCounter())
        if hc.counter == 0
                firstPass(e, hc)
        end

        localMem = false
        globalindicator = "__global "
        if localMem
                globalindicator = ""
        end
        if hc.exprList[e.hash].count > 0
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


                # nDim = length(filter(x->x[2]== e.commonVariables.name, wSeed.input))
                # nParam =length(filter(x->x[2]==e.individualVariables.name, wSeed.input))
                nDim = e.commonVariables.counter
                nParam = e.individualVariables.counter

                workMemoryOffset = string("+ gid*(2*",nDim,"+2*",nParam," + 2)*nRepeat")
                outPutOffset = string("+ gid*(",nDim,"+",nParam,"*nRepeat)")
                if localMem
                        workMemoryOffset = ""
                        outPutOffset = ""
                end

                myFunctionName = string("calculateDerivedParametericAverageInput",hc.counter)

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



                derivedParametericAvgFunctionHash = myFunctionName
                hc.exprList[derivedParametericAvgFunctionHash] = (-7,derivedParametericAvgFunction, "function")
                hc.exprList[derivedParametericAvgFunctionHash].count = 0

                derivedParAvgFunc2Declaration = string("void diffAvg(__global float * ",e.commonVariables.name,",",
                                            "__global float * ",e.individualVariables.name,", int nRepeat, ",globalindicator," float *workmemory, ",globalindicator," float * output);\n")
                derivedParAvgFunc2 = string("void diffAvg(__global float * ",e.commonVariables.name,",",
                                            "__global float * ",e.individualVariables.name,", int nRepeat, ",globalindicator," float *workmemory, ",globalindicator," float * output){
                                                const int gid = get_global_id(0);

                                                //__global float temp[(2*",nDim,"+2*",nParam," + 2)*nRepeat];
                                                ",globalindicator," float *xSeed = workmemory ",workMemoryOffset,";
                                                ",globalindicator," float *wSeed = xSeed+nRepeat;
                                                ",globalindicator," float *xSeedDerivedC = wSeed+nRepeat;
                                                ",globalindicator," float *xSeedDerivedP = xSeedDerivedC+",nDim,"*nRepeat;
                                                ",globalindicator," float *wSeedDerivedC = xSeedDerivedP+",nParam,"*nRepeat;
                                                ",globalindicator," float *wSeedDerivedP = wSeedDerivedC+",nDim,"*nRepeat;

                                                ",globalindicator," float *derivedC = output ",outPutOffset,";
                                                ",globalindicator," float *derivedP = derivedC + ",nDim,";
                                                for(int i=0;i<(",nDim,"+",nParam,"*nRepeat);i++)
                                                        {output[i] = 0;}


                                                __global float *ourCoord = &",e.commonVariables.name,"[gid*",nDim,"];
                                                float sumW = 0.0f;
                                                float avg = 0.0f;
                                                for(int i=0; i<nRepeat; i++) {

                                                                __global float * myParams = &",e.individualVariables.name,"[i*",nParam,"];
                                                                ",myFunctionName,"(ourCoord,myParams,",
                                                                "xSeed+i,wSeed+i,",
                                                                "xSeedDerivedC+i*",nDim,", xSeedDerivedP + i*",nParam,",",
                                                                "wSeedDerivedC+i*",nDim,", wSeedDerivedP + i*",nParam,");\n",
                                                                "sumW += wSeed[i];\n",
                                                                "avg += wSeed[i]*xSeed[i];\n",
                                                "}\n",
                                                "avg = avg/sumW;\n",
                                                "for(int i=0; i<nRepeat; i++){\n",
                                                                "for(int j=",xIMin,";j<=",xIMax,";j++) {
                                                                        derivedP[i*",nParam,"+j] = wSeed[i]*xSeedDerivedP[j+",nParam,"*i]/sumW;
                                                                }\n",
                                                                "for(int j=",wIMin,";j<=",wIMax,";j++) {
                                                                        derivedP[i*",nParam,"+j] = wSeedDerivedP[i*",nParam,"+j]*(xSeed[i] - avg)/sumW;
                                                                }\n",
                                                                "for(int j=",xCMin,";j<=",xCMax,";j++){
                                                                        derivedC[j] += wSeed[i]*xSeedDerivedC[j+i*",nDim,"];
                                                                }\n",
                                                                "for(int j=",wCMin,";j<=",wCMax,";j++){
                                                                        derivedC[j] += wSeedDerivedC[j+i*",nDim,"]*(xSeed[i] - avg ) ;
                                                                }\n",
                                                        "}\n",
                                                "\n",
                                                "for(int i=0; i<",nDim,"; i++ ) {derivedC[i]/=sumW;}\n",

                                            "}\n",


                                            )


                    derivedParametericAvgFunctionHash = myFunctionName
                    hc.exprList[derivedParametericAvgFunctionHash] = (-7,string(derivedParametericAvgFunction,derivedParAvgFunc2), "function")
                    hc.exprList[derivedParametericAvgFunctionHash].count = 0

                    derivedParametericAvgFunctionDeclarationHash = string(myFunctionName,"declaration")
                    hc.exprList[derivedParametericAvgFunctionDeclarationHash] = (-137,string(derivedParametericAvgFunctionDeclaration,derivedParAvgFunc2Declaration),"function declaration")
                    hc.exprList[derivedParametericAvgFunctionDeclarationHash].count =0
                myExpression = ""

                if hc.exprList[e.hash].count >1
                        eval = hc.exprList[e.hash].eval
                        defs = string("float ", eval, "=", myExpression,";\n")
                else
                        defs = ""
                        eval = string(myExpression )
                end

                hc.exprList[e.hash].count = -hc.exprList[e.hash].count

        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end
        literal(defs, eval, dataAttributes(hc) )

end

function deRef(e::parametericAvg)
        pr = parametericAvg(constIndex(0, "dummy", 1.0),
                        constIndex(0, "dummy", 1.0),
                        constIndex(0, "dummy", 1.0),
                        e.commonVariables,
                        e.individualVariables)

        pr.nRepeat = deRef(e.nRepeat)
        pr.xSeed = deRef(e.xSeed)
        pr.wSeed = deRef(e.wSeed)
        pr.id = deepcopy(e.id)
        pr.name = deepcopy(e.name)
        pr.hash = deepcopy(e.hash)
        return pr
end

function deRef(e::derivedParametericAvg)
        pr = derivedParametericAvg(constIndex(0, "dummy", 1.0),
                        constIndex(0, "dummy", 1.0),
                        constIndex(0, "dummy", 1.0),
                        X(),
                        X())

        pr.commonVariables = e.commonVariables
        pr.individualVariables = e.individualVariables
        pr.wSeedDerived = (deRef.(e.wSeedDerived[1]), deRef.(e.wSeedDerived[2]))
        pr.xSeedDerived = (deRef.(e.xSeedDerived[1]), deRef.(e.xSeedDerived[2]))
        pr.nRepeat = deRef(e.nRepeat)
        pr.xSeed = deRef(e.xSeed)
        pr.wSeed = deRef(e.wSeed)
        pr.id = deepcopy(e.id)
        pr.name = deepcopy(e.name)
        pr.hash = deepcopy(e.hash)
        return pr
end
