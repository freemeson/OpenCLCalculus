function writeFunctionFile(jl::literal, functionname::String, filename::String)
        extraInputVars = filter(x->x[3] == -1,  jl.input)
        if length(extraInputVars) != 0
                extraInputVarString = string([ string(",",e[2],"::",e[5]) for e in extraInputVars]...)
        else
                extraInputVarString = ""
        end

        open(filename, "w") do f
                write(f, string("function ",functionname,"(x::Vector,p::Vector",extraInputVarString,")\n"))
                write(f, jl.definitions)
                write(f, string("return ", jl.expression,"\n") )
                write(f, "end\n")
        end
        return nothing
end

function writeOpenCLPrg(jl::literal, kernelName::String)
        nDim = maximum([i[3] for i in filter(x->x[2]=="x",jl.input)])
        nParam = maximum([i[3] for i in filter(x->x[2]=="p",jl.input)])
        println("nDim = ", nDim)
        println("nParam = ", nParam)

        ffoutput = findfirst(i->i[4] == "output", jl.input)
        if ffoutput!=nothing
                nOutPut = -jl.input[ffoutput][3]
                println("nOutPut = ", nOutPut)
                println("nPeak = ", nOutPut - nParam*nDim-nDim*nDim)
                resultStr = string("float outPut[",nOutPut,"] = ",jl.expression,";\n",
                                        "for(int i=0; i<",nOutPut,"; i++){ result[",nOutPut,"*gid + i] = outPut[i]; } \n"  )
        else
                resultStr = string("result[gid] = ", jl.expression,";\n")
        end
        string("
            #pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
            #pragma OPENCL EXTENSION cl_khr_fp64 : disable
        ",
                [j[2] for j in filter(x->x[4]=="function", jl.input)]...,
                "__kernel void ",kernelName,"(__global float* data,
                                                __global float *p,
                                                __global float *result){
                        const int gid = get_global_id(0);
                        const size_t nDim=",nDim,";
                        __global float* x = &data[gid*nDim];\n;
                                                \n",
                                                jl.definitions,
                                                resultStr,
                                                "}\n"
         )
end

function writeOpenCLFunction(jl::literal, functionName::String)
        nDim = maximum([i[3] for i in filter(x->x[2]=="x",jl.input)])
        nParam = maximum([i[3] for i in filter(x->x[2]=="p",jl.input)])

        println("nDim = ", nDim)
        println("nParam = ", nParam)

        ffoutput = findfirst(i->i[4] == "output", jl.input)
        if ffoutput!=nothing
                nOutPut = -jl.input[ffoutput][3]
                println("nOutPut = ", nOutPut)
                println("nPeak = ", nOutPut - nParam*nDim-nDim*nDim)

                resultStr = string("float outPut[",nOutPut,"] = ",jl.expression,";\n",
                                   "for(int i=0; i<",nOutPut,"; i++){ result[i] = outPut[i]; } \n"  )
        else
                nOutPut = 1
                resultStr = string("result[gid] = ", jl.expression,";\n")
        end
        string("
            #pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable
            #pragma OPENCL EXTENSION cl_khr_fp64 : disable\n",

                [j[2] for j in filter(x->x[4]=="function", jl.input)]...,
                "void ",functionName,"(__global float* x,__global float *p, float *result);\n",
                "void ",functionName,"(__global float* x,
                                       __global float *p,
                                        float *result
                                        ){\n",
                                        # " const size_t nDim=",nDim,";\n",
                                        # " const size_t nParam=", nParam, ";\n",
                                        # " const size_t nOutPut=", nOutPut,";\n",
                                        # " const size_t nPeak=",nOutPut-nParam*nDim-nDim*nDim,";\n\n",

                                                jl.definitions,
                                                resultStr,
                                                "}\n"
         )
end

### example code
# function juliaLiteral(f::sigmaDistance)
#         #if the parameters are not expressions
#         xStr = string("x[1:",f.nDimX,"]")
#         mStr = string("p[1:",f.nDimX,"]")
#         sStr = string("p[",f.nDimX+1,":",f.nDimP,"]")
#         #preProcessingStr
#         codeStr = string("sum(((",xStr,"-",mStr,")./",sStr,").^2)")
# end
#
#
# expr1 = sigmaDistance(4)
# juliaLiteral(expr1)
#
#
#
#
# function juliaFunction(f::sigmaDistance)
#         functionHeaderStr = string("function ", typeof(f), "(x::Vector,p::Vector)")
#         codeStr = juliaLiteral(f)
#         functionStr = string(functionHeaderStr,"\n",codeStr,"\n","end\n")
# end
#
# function sigmaDistance(x::Vector,p::Vector)
#         sum(((x[1:4]-p[1:4])./p[5:8]).^2)
# end

#sigmaDistance(collect(1.0:4.0), collect(9.0:16.0))









#need to add comments, attributes to dataIndices, so later I can determine their min,max properties and generate boundary values, boundary checks, random generators for them
# types like
# - amplitude [0, 1]
# - polynome or affine transformation factors [-inf, +inf]
# - width, deviation, sigma paramaters [+0., inf]

#need more PDF and CDF
# - erf
# - integral of logistic

function marginalCDF(nGauss, nDL, nDim)

        #array for pdf = dotprod(prod(pdf), amplitudes)
        #each pdf is univariate, prod(pdf) = multivariate
        #marginalCDFs are from marginal pdf mpdf[i] = dotprod(prod(1..i, pdf[k,i]), amplitude[k])
        #  the cdf is produced by the integral of the ith dimension mcdf[i] = dotprod( prod(1..i-1,pdf[k,i]*cdf(k,i), amplitude[k]  )
        #  the conditionalMarginalCDF = mcdf[i]/mpdf[i-1]
        #for the ghost entropy, we need the derivatives d/dp mcd[i]/mpdf[i-1]
        x = X("x", "data")
        p = X("p", "parameter")
        cOne = constIndex(0, "mCDFgeneratorOne", 1.0)
        cF = constIndex(0, "mCDFgeneratorFactor", 1.0/3.0)

        function affine(nDim::Int)
                return dotProdExpr( vcat(cOne, x[1:nDim-1]), newVars(p, nDim, "affinePolynome") )
        end
        amplitudes = newVars(p, nGauss, "amplitude")


        function integral(ge::gaussExpr)
                return erfExpr(ge.pos, ge.p, ge.s)
        end

        function integral(de::doubleLogisticExpr)
                return cdfDoubleLogisticExpr(de.pos, de.p, de.w, de.s)
        end

        peaks = Array{ScalarDataIndex,2}(undef, (nGauss, nDim))
        peakProb = Vector{ScalarDataIndex}(undef,nGauss)
        logPeakProb = Vector{ScalarDataIndex}(undef,nGauss)

        for i in 1:nGauss
                for j in 1:nDim
                        #if true
                        if false #j<nDim #&& false
                                w = newVar(p, "width")
                                peaks[i,j] = doubleLogisticExpr( x[j], affine(j), w, prodExpr([w,cF]) )
                                #peaks[i,j] = doubleLogisticExpr( x[j], affine(j), newVar(p, "width"), newVar(p, "smear") )
                        else
                                peaks[i,j] = gaussExpr( x[j], affine(j), newVar(p,"sigma") )
                        end
                end
                peakProb[i] = prodExpr(vcat(amplitudes[i], peaks[i,:] ))
#                logPeakProb[i] = sumExpr( [logExpr(a) for a in  peaks[i,:]] ) #  vcat(amplitudes[i], peaks[i,:] ) ]  )
                logPeakProb[i] = sumExpr( logExpr.(vcat(amplitudes[i], peaks[i,:])))
         end


        probability = sumExpr(peakProb)

        #return juliaLiteral(simplifyExpr(probability))
        cdfs = integral.(peaks)

        marginalCDFs = Vector{ScalarDataIndex}(undef, nDim)
        derivedMarginalCDF = Vector{Tuple{Vector{ScalarDataIndex}, Vector{ScalarDataIndex}}}(undef, nDim)

        for i in 1:nDim
                mCDF = Vector{ScalarDataIndex}(undef, nGauss)
                mPDF = Vector{ScalarDataIndex}(undef, nGauss)
                for j in 1:nGauss
                        mPDF[j] = prodExpr(vcat(amplitudes[j], peaks[j, 1:i-1] ))
                        mCDF[j] = prodExpr(vcat(amplitudes[j], cdfs[j,i], peaks[j, 1:i-1] ))
                end
                marginalCDFs[i] = divExpr(sumExpr(mCDF), sumExpr(mPDF))
                derivedMarginalCDF[i] = derive(marginalCDFs[i],x,p)
        end

        nParams = p.counter
        dFdx = Array{ScalarDataIndex}(undef,(nDim,nDim))
        dFdp = Array{ScalarDataIndex}(undef,(nDim,nParams))
        for i in 1:nDim
                dFdx[i,:] .= derivedMarginalCDF[i][1]
                dFdp[i,:] .= derivedMarginalCDF[i][2]
        end

        vectorExpression = simplifyExpr.(vcat( dFdx[:], dFdp[:], logPeakProb))
        #return juliaLiteral(dFdx[2])
        return juliaLiteral(vectorExpression), juliaLiteral(simplifyExpr(probability))

end

function marginalCDF_opencl(nPeaks, nDim, nonGauss = false, smoothFactor = 1.0/3.0 )

        #array for pdf = dotprod(prod(pdf), amplitudes)
        #each pdf is univariate, prod(pdf) = multivariate
        #marginalCDFs are from marginal pdf mpdf[i] = dotprod(prod(1..i, pdf[k,i]), amplitude[k])
        #  the cdf is produced by the integral of the ith dimension mcdf[i] = dotprod( prod(1..i-1,pdf[k,i]*cdf(k,i), amplitude[k]  )
        #  the conditionalMarginalCDF = mcdf[i]/mpdf[i-1]
        #for the ghost entropy, we need the derivatives d/dp mcd[i]/mpdf[i-1]
        x = X("x", "data")
        p = X("p", "parameter")
        cOne = constIndex(0, "mCDFgeneratorOne", 1.0)
        cF = constIndex(0, "mCDFgeneratorFactor", smoothFactor)

        function affine(nDim::Int)
                return dotProdExpr( vcat(cOne, x[1:nDim-1]), newVars(p, nDim, "affinePolynome") )
        end
        amplitudes = newVars(p, nPeaks, "amplitude")


        function integral(ge::gaussExpr)
                return erfExpr(ge.pos, ge.p, ge.s)
        end

        function integral(de::doubleLogisticExpr)
                return cdfDoubleLogisticExpr(de.pos, de.p, de.w, de.s)
        end

        peaks = Array{ScalarDataIndex,2}(undef, (nPeaks, nDim))
        peakProb = Vector{ScalarDataIndex}(undef,nPeaks)
        logPeakProb = Vector{ScalarDataIndex}(undef,nPeaks)

        for i in 1:nPeaks
                for j in 1:nDim
                        if  nonGauss
                                w = newVar(p, "width")
                                peaks[i,j] = doubleLogisticExpr( x[j], affine(j), w, prodExpr([w,cF]) )
                        else
                                peaks[i,j] = gaussExpr( x[j], affine(j), newVar(p,"sigma") )
                        end
                end
                peakProb[i] = prodExpr(vcat(amplitudes[i], peaks[i,:] ))
                logPeakProb[i] = sumExpr( logExpr.(vcat(amplitudes[i], peaks[i,:])))#must manually normalize aplitudes later
         end


        probability = divExpr(sumExpr(peakProb), sumExpr(amplitudes))
        #probability = sumExpr(peakProb)
        cdfs = integral.(peaks)

        marginalCDFs = Vector{ScalarDataIndex}(undef, nDim)
        derivedMarginalCDF = Vector{Tuple{Vector{ScalarDataIndex}, Vector{ScalarDataIndex}}}(undef, nDim)


        for i in 1:nDim
                mCDF = Vector{ScalarDataIndex}(undef, nPeaks)
                mPDF = Vector{ScalarDataIndex}(undef, nPeaks)
                for j in 1:nPeaks
                        mPDF[j] = prodExpr(vcat(amplitudes[j], peaks[j, 1:i-1] ))
                        mCDF[j] = prodExpr(vcat(amplitudes[j], cdfs[j,i], peaks[j, 1:i-1] ))
                end
                marginalCDFs[i] = divExpr(sumExpr(mCDF), sumExpr(mPDF))
                derivedMarginalCDF[i] = derive(marginalCDFs[i],x,p)

        end

        lastPDF = Vector{ScalarDataIndex}(undef, nPeaks)
        penultimatePDF = Vector{ScalarDataIndex}(undef, nPeaks)

        for i in 1:nPeaks
                lastPDF[i] = prodExpr(vcat( amplitudes[i], peaks[i,1:end] )  )
                penultimatePDF[i] = prodExpr(vcat( amplitudes[i], peaks[i,1:end-1] ) )
        end
        mPDFLast = divExpr( sumExpr(lastPDF), sumExpr(penultimatePDF ))
#ß        return mPDFLast
        nParams = p.counter
        dFdx = Array{ScalarDataIndex}(undef,(nDim,nDim))
        dFdp = Array{ScalarDataIndex}(undef,(nDim,nParams))
        for i in 1:nDim
                dFdx[i,:] .= derivedMarginalCDF[i][1]
                dFdp[i,:] .= derivedMarginalCDF[i][2]
        end

        vectorExpression = simplifyExpr.(vcat( dFdx[:], dFdp[:], logPeakProb))
        #return vectorExpression, probability
        return openclLiteral(vectorExpression), openclLiteral(simplifyExpr(probability)), openclLiteral(simplifyExpr(mPDFLast)), juliaLiteral(simplifyExpr(probability))
end

function regression_opencl(nDim::Int)
        x = X("x", "data")
        p = X("p", "parameter")
        cOne = constIndex(0, "mCDFgeneratorOne", 1.0)
        cF = constIndex(0, "mCDFgeneratorFactor", 1.0/3.0)

        function affine(nDim::Int)
                #return dotProdExpr( vcat(cOne, x[1:nDim-1]), newVars(p, nDim, "affinePolynome") )
                if nDim != 1
                        #xPow2 = Vector{ScalarDataIndex}(undef, nDim -1)
                        #xPow3 = Vector{ScalarDataIndex}(undef, nDim -1)
                        #xPow4 = Vector{ScalarDataIndex}(undef, nDim -1)
                        for i=1:nDim-1
                                #xPow2[i] = prodExpr([x[i],x[i]])
                                #xPow3[i] = prodExpr([x[i],x[i],x[i]])
                                #xPow4[i] = prodExpr([x[i],x[i],x[i],x[i]])
                        end
                        return sumExpr([newVar(p,"mean"), dotProdExpr( x[1:nDim-1], newVars(p, nDim-1, "affinePolynome"))
                                        #,dotProdExpr(xPow2,newVars(p, nDim-1, "affinePolynome"))
                                        #,dotProdExpr(xPow3,newVars(p, nDim-1, "affinePolynome"))
                                        #,dotProdExpr(xPow4,newVars(p, nDim-1, "affinePolynome"))
                                        ])

                        #return prodExpr([newVar(p,"affinePolynome"),logisticExpr(x[1], newVar(p,"affinePolynome"),newVar(p,"affinePolynome"))])
                else
                        return newVar(p,"mean")
                end
        end

        function integral(ge::gaussExpr)
                return erfExpr(ge.pos, ge.p, ge.s)
        end

        function integral(de::doubleLogisticExpr)
                return cdfDoubleLogisticExpr(de.pos, de.p, de.w, de.s)
        end

        peaks = Vector{ScalarDataIndex}(undef, nDim)
        amplitude = newVar(p, "amplitude")
        #amplitude = prodExpr([a,a])
        for j in 1:nDim
                if true
                        #affine(j)
                        peaks[j] = gaussExpr( x[j], affine(j), newVar(p,"sigma") )
                else
                # w = newVar(p, "width")
                # peaks[j] = doubleLogisticExpr( x[j], affine(j), w, prodExpr([w,cF]) )
                        #peaks[j] = doubleLogisticExpr( x[j], affine(j), newVar(p,"sigma"), newVar(p,"smear") )
                        me = affine(j)#newVar(p,"mean")
                        si = newVar(p,"sigma")
                        #sm = newVar(p,"smear")#
                        sm = prodExpr([si, cF])
                        peaks[j] = doubleLogisticExpr( x[j], me , si , sm)
                end
        end
        logPeakProb  = simplifyExpr(sumExpr(logExpr.(vcat(amplitude, peaks) ) )) #must manually normalize aplitudes
        logPeakProb.name = "logProbability"
        logPeakProbCond = simplifyExpr(sumExpr(logExpr.(vcat(amplitude, peaks[1:end-1]) ) ))
        logPeakProbCond.name = "logProbabilityCond"

        nRepeat = dataIndex(0,"nRepeat", "nPeaks")

        lastConditionalPDF = parametericAvg(peaks[nDim], simplifyExpr(prodExpr(vcat(amplitude, peaks[1:nDim-1] ) ) ), nRepeat, x,p )

        cdfs = sumExpr([integral(peaks[end]), newVar(p,"offset") ])
        marginalCDFs = Vector{derivedParametericAvg}(undef, 1)

        marginalCDFs[1] = derivedParametericAvg(cdfs,prodExpr(vcat(amplitude,peaks[1:end-1])),nRepeat,x,p,0)

        return openclVectorFunction(marginalCDFs), openclVectorFunction( logPeakProb, x, p ), openclVectorFunction( logPeakProbCond, x, p ), juliaLiteral(lastConditionalPDF)
end


function vector_marginalCDF_opencl(nDim::Int, probabilityOutputSelection::Symbol = :FullProbability )

        x = X("x", "data")
        p = X("p", "parameter")
        cOne = constIndex(0, "mCDFgeneratorOne", 1.0)
        cF = constIndex(0, "mCDFgeneratorFactor", 1.0/3.0)

        function affine(nDim::Int)
                #return dotProdExpr( vcat(cOne, x[1:nDim-1]), newVars(p, nDim, "affinePolynome") )
                if nDim != 1
                        #xPow2 = Vector{ScalarDataIndex}(undef, nDim -1)
                        #xPow3 = Vector{ScalarDataIndex}(undef, nDim -1)
                        #xPow4 = Vector{ScalarDataIndex}(undef, nDim -1)
                        for i=1:nDim-1
                                #xPow2[i] = prodExpr([x[i],x[i]])
                                #xPow3[i] = prodExpr([x[i],x[i],x[i]])
                                #xPow4[i] = prodExpr([x[i],x[i],x[i],x[i]])
                        end
                        return sumExpr([newVar(p,"mean"), dotProdExpr( x[1:nDim-1], newVars(p, nDim-1, "affinePolynome"))
                                        #,dotProdExpr(xPow2,newVars(p, nDim-1, "affinePolynome"))
                                        #,dotProdExpr(xPow3,newVars(p, nDim-1, "affinePolynome"))
                                        #,dotProdExpr(xPow4,newVars(p, nDim-1, "affinePolynome"))
                                        ])

                        #return prodExpr([newVar(p,"affinePolynome"),logisticExpr(x[1], newVar(p,"affinePolynome"),newVar(p,"affinePolynome"))])
                else
                        return newVar(p,"mean")
                end
        end

        function integral(ge::gaussExpr)
                return erfExpr(ge.pos, ge.p, ge.s)
        end

        function integral(de::doubleLogisticExpr)
                return cdfDoubleLogisticExpr(de.pos, de.p, de.w, de.s)
        end

        peaks = Vector{ScalarDataIndex}(undef, nDim)
        amplitude = newVar(p, "amplitude")
        #amplitude = prodExpr([a,a])
        for j in 1:nDim
                if  true
                        peaks[j] = gaussExpr( x[j], affine(j), newVar(p,"sigma") )
                else
                # w = newVar(p, "width")
                # peaks[j] = doubleLogisticExpr( x[j], affine(j), w, prodExpr([w,cF]) )
                        #peaks[j] = doubleLogisticExpr( x[j], affine(j), newVar(p,"sigma"), newVar(p,"smear") )
                        me = newVar(p,"mean")
                        si = newVar(p,"sigma")
                        sm = prodExpr([si, cF])
                        peaks[j] = doubleLogisticExpr( x[j], me , si , sm)
                end
        end
        logPeakProb  = simplifyExpr(sumExpr(logExpr.(vcat(amplitude, peaks) ) )) #must manually normalize aplitudes

        nRepeat = dataIndex(0,"nRepeat", "nPeaks")

        lastConditionalPDF = parametericAvg(peaks[nDim], simplifyExpr(prodExpr(vcat(amplitude, peaks[1:nDim-1] ) ) ), nRepeat, x,p )
        pdf = parametericAvg(simplifyExpr(prodExpr( peaks[1:nDim] ) ),amplitude , nRepeat, x,p )

        if probabilityOutputSelection == :FullProbability
                probOutput = pdf
        end
        if probabilityOutputSelection == :ConditionalProbability
                probOutput = lastConditionalPDF
        end

#        return juliaLiteral(lastMarginalPDF)
        cdfs = integral.(peaks)
        marginalCDFs = Vector{derivedParametericAvg}(undef, nDim)

        for i in 1:nDim
                marginalCDFs[i] = derivedParametericAvg(cdfs[i],prodExpr(vcat(amplitude,peaks[1:i-1])),nRepeat,x,p,i)
        end

        #return marginalCDFs[nDim]

        return openclVectorFunction(marginalCDFs), openclVectorFunction( logPeakProb, x, p ), juliaLiteral(probOutput)#,juliaLiteral(lastConditionalPDF)


end

function exampleDiffParmAvgVector(len::Int, nRepeat::Int = 1)
        x = X("x");
        p = X("p");
        ex = Vector{derivedParametericAvg}(undef, len)
        xSeed = Vector{ScalarDataIndex}(undef, len)
        wSeed = Vector{ScalarDataIndex}(undef, len)

        ampl = newVar(p, "amplitude")
        for i in 1:len
                xSeed[i] = dotProdExpr( x[1:i], newVars(p,i) )
                wSeed[i] = dotProdExpr( x[1:i], newVars(p,i) )

        end

        nRep = dataIndex(0,"nRepeat")
        for i in 1:len
                ex[i] = derivedParametericAvg(xSeed[i], prodExpr(vcat(ampl,wSeed[1:i])), nRep, x,p)
        end

        ##
        xr = X("x");
        pr = X("p");
        exr = Vector{ScalarDataIndex}(undef, len )
        xrep= Array{ScalarDataIndex,2}(undef, (nRepeat,len))
        wrep = Array{ScalarDataIndex,2}(undef, (nRepeat,len))
        amplr = Vector{ScalarDataIndex}(undef, nRepeat)

        wm = ScalarDataIndex[]

        for j in 1:nRepeat
                amplr[j] = newVar(pr, "amplitude")
                for i in 1:len
                        xrep[j,i] = dotProdExpr( xr[1:i], newVars(pr,i) )
                        wrep[j,i] = dotProdExpr( xr[1:i], newVars(pr,i) )
                end
        end

        nRep = dataIndex(0,"nRepeat")
        for i in 1:len
                temp = Vector{ScalarDataIndex}(undef, nRepeat)
                for j in 1:nRepeat
                        temp[j] = prodExpr(vcat(amplr[j],wrep[j,1:i]))
                        # wt[j,i] = temp[j]
                        # append!(wm, xr[] )
                end
                a = dotProdExpr(xrep[:,i], temp)
                b = sumExpr(temp)
                exr[i] = divExpr(a , b)


                append!(wm,xrep[:,i])
                append!(wm,temp)
                wm2 = ScalarDataIndex[]

                for j in 1:nRepeat
                        c = derive(xrep[j,i], xr, pr)
                        append!(wm,  c[1])
                        append!(wm2, c[2])
                end
                for j in 1:nRepeat
                        c = derive(temp[j], xr, pr)
                        append!(wm, c[1] )
                        append!(wm2, c[2])
                end

                append!(wm,wm2)
                # da = derive(a, xr, pr)
                # db = derive(b, xr, pr)

        end


        dFdx = Array{ScalarDataIndex}(undef,(len,len))
        dFdp = Array{ScalarDataIndex}(undef,(len,pr.counter))
        for i in 1:len
                d = derive(exr[i], xr, pr)
                dFdx[i,:] .= d[1]
                dFdp[i,:] .= d[2]

        end

        return ex, exr, simplifyExpr.( vcat(dFdx[:], dFdp[:] )), simplifyExpr.(wm)
end

function writeOpenclFile(prg::String, filename::String)



        open(filename, "w") do f
                write(f, prg)
        end
        return nothing
end
