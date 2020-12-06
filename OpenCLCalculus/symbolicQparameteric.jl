using OpenCL
const cl = OpenCL
using Plots
using Distributions
using Calculus
using NLopt

function createQKernel(pmcdf::literal, plogProb::literal)

    #pmcdf,plogProb = vector_marginalCDF_opencl(2)

    localMem = false
    nPeakAssume = 32
    nDim   = maximum([i[3] for i in filter(x->x[2]=="x",pmcdf.input)])
    nParamSinglePeak = maximum([i[3] for i in filter(x->x[2]=="p",pmcdf.input)])


    pmcdf_function_declarations = filter(x->x[4]=="function declaration", pmcdf.input)
    pmcdf_declarations = string([a[2] for a in pmcdf_function_declarations]...)
    plogProb_function_declarations = filter(x->x[4]=="function declaration" && findfirst(y->y[1]==x[1], pmcdf_function_declarations) == nothing, plogProb.input)
    plogProb_declarations = string([a[2] for a in plogProb_function_declarations]...)


    pmcdf_functions = filter(x->x[4]=="function", pmcdf.input)
    pmcdf_prg = string([a[2] for a in pmcdf_functions]...)

    plogProb_functions = filter(x->x[4]=="function" && findfirst(y->y[1]==x[1], pmcdf_functions) == nothing, plogProb.input)
    plogProb_prg = string([a[2] for a in plogProb_functions]...)
    # pmcdf.expression
    # plogProb.expression

    mcdfObjName = filter(x->x[3]==-2, pmcdf.input)[1][2]
    logProbObjName = filter(x->x[3]==-2, plogProb.input)[1][2]

    extraInput = string(",
   __global float *workmemory,
   __global float *",mcdfObjName,",
   __global float *",logProbObjName)
   globalIndicator = "__global "
   objDeclarations = ""
   mcdfObjOffset = "+gid*mcdfObjSize"
   logProbObjOffset = "+gid*nRepeat"
    if localMem
        extraInput = ""
        globalIndicator = ""
        objDeclarations = string("float workmemory[",nPeakAssume*(2*nDim+2*nParamSinglePeak+2),"];
                                    float ",mcdfObjName,"[",nDim*(nDim+nParamSinglePeak*nPeakAssume),"];
                                    float ",logProbObjName,"[",nPeakAssume,"];\n")
        mcdfObjOffset = ""
        logProbObjOffset = ""
    end

    Q_opencl = string("

        __kernel void Q(__global float *data,
                        __global float *p,
                        __global float *Dp,
                        __global float *sumlogp,
                         int nRepeat,
                         bool ghost
                         ",extraInput,") {\n",
                        " const int nDim=",nDim,";\n",
                        " const int nParam=", nParamSinglePeak, "*nRepeat;\n",

                            "const int gid = get_global_id(0);
                            __global float* x = &data[gid*nDim];

                            ",objDeclarations,"
                            if(ghost)
                            {",pmcdf.definitions,"}\n
                            ",plogProb.definitions,"\n

                            const int mcdfObjSize = ",nDim,"*(",nDim," + nParam);
                            ",globalIndicator," float * dFdx = ",mcdfObjName,mcdfObjOffset,";
                            ",globalIndicator," float * dFdp = dFdx + nDim*nDim;
                            ",globalIndicator," float * logprob = ",logProbObjName,logProbObjOffset,";

                            float amplNorm = 0.0f;
                            for(int i=0; i<nRepeat;i++)
                                   {amplNorm+=p[i*",nParamSinglePeak,"] ;}
                            amplNorm = log((float)(amplNorm));
                            for(int i=0; i<nRepeat; i++)
                                {logprob[i] -= amplNorm;}

                            float y[",nDim,"] ;
                            float pert[",nDim,"] = {0.0f};
                            float q = 0.0;
                            if(ghost){
                                for(int i=0;i!=nParam;i++) {
                                        bool rankAlertNaN = false;
                                        bool rankAlertInf = false;
                                        float tempNaN = dFdp[i*nDim]/dFdx[0];
                                        y[0] = tempNaN;
                                        if(isnan(tempNaN))
                                            {
                                                rankAlertNaN = true;
                                            }
                                        if(isinf(tempNaN))
                                            {
                                                rankAlertInf = true;
                                            }
                                        for(int j=1;j<nDim; j++) {
                                               y[j] = dFdp[ i*nDim + j];
                                               for(int k=0; k<j;k++ ) {
                                                       y[j] -= y[k]* dFdx[k*nDim + j];
                                                   }
                                                float tempNaN = y[j]/dFdx[ j + j*nDim ];
                                                y[j] = tempNaN;
                                                if(isnan(tempNaN))
                                                    {
                                                    rankAlertNaN = true;
                                                }
                                                if(isinf(tempNaN))
                                                    {rankAlertInf = true;}
                                                //y[j] = y[j]/dFdx[ j + j*nDim ];
                                            }
                                        if(rankAlertInf || rankAlertNaN)
                                            {
                                                for(int j=0; j<nDim; j++)
                                                    {
                                                        if(isinf(y[j]) ||isnan(y[j]))
                                                            {
                                                                //y[j] = 1.0f;
                                                                y[j] = 1e-4f;
                                                            } else {
                                                                //y[j] = 0.0f;
                                                            }
                                                    }
                                            } else {
                                            if(rankAlertNaN) {
                                                /*for(int j=0; j<nDim; j++)
                                                    {
                                                        if(isnan(y[j]))
                                                            {
                                                                y[j] = 1.0f;
                                                            } else {
                                                                y[j] = 0.0f;
                                                            }
                                                    }*/

                                                }

                                            }
                                        for(int j=0; j< nDim; j++)
                                            {
                                              pert[j] += fabs((float)(y[j]))*exp((float)(Dp[i]));
    //                                        sumlogp[gid*nParam*nDim + i*nDim + j] = y[j];
                                            }

                                    }

    //                            for(int i=0; i!= nDim*nDim + nDim*nParam; i++) {
    //                                sumlogp[gid*(nParam*nDim+nDim*nDim)+i] = output[i];
    //                            }
                                float perturbation = 1.0f;
                                for(int i=0; i!= nDim; i++)
                                    {perturbation -= pert[i];}


                                if (perturbation <= 0) {
                                        q = 9.9e36;
                                    } else {
                                        q = -log((float)(perturbation) );
                                    }

                            }
                            float maxprob = logprob[0];
                            int maxindex = 0;
                            for(int i=1; i<nRepeat; i++) {
                                if(logprob[i] > maxprob) {
                                        maxindex = i;
                                        maxprob = logprob[i];
                                    }
                            }

                            if (maxprob < -90.0f)
                                {
                                    q-= maxprob;
                                    for(int i=0; i<nRepeat; i++) {
                                        if(i!=maxindex)
                                            {
                                                q -= exp(logprob[i] - maxprob);
                                            }
                                        }
                                } else {
                                    float prob = 0.0;
                                    for(int i=0; i!= nRepeat; i++)
                                        {
                                          prob += exp( (float)( logprob[i] ) );
                                        }
                                    q-=log((float)(prob));
                                }

//                                for(int i=0; i!= nParam; i++)
//                                    {q-=Dp[i];}

                            if (isinf(q) ){ q=9.9e36; }
                            sumlogp[gid] = q;
                        }



      ")

      #pragmaString = "#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable\n#pragma OPENCL EXTENSION cl_khr_fp64 : disable\n"

      # device, ctx, queue = OpenCL.cl.create_compute_context()
      # prg = OpenCL.cl.Program(ctx, source = string(pmcdf_declarations, plogProb_declarations, pmcdf_prg, plogProb_prg, Q_opencl) ) |> OpenCL.cl.build!

      return string(pmcdf_declarations, plogProb_declarations, pmcdf_prg, plogProb_prg, Q_opencl)
end

function createQRegressionKernel(pmcdf::literal, plogProb::literal, plogProbCond::literal)

    #pmcdf,plogProb = vector_marginalCDF_opencl(2)

    localMem = false
    nPeakAssume = 32
    nDim   = maximum([i[3] for i in filter(x->x[2]=="x",pmcdf.input)])
    nParamSinglePeak = maximum([i[3] for i in filter(x->x[2]=="p",pmcdf.input)])


    pmcdf_function_declarations = filter(x->x[4]=="function declaration", pmcdf.input)
    pmcdf_declarations = string([a[2] for a in pmcdf_function_declarations]...)
    plogProb_function_declarations = filter(x->x[4]=="function declaration" && findfirst(y->y[1]==x[1], pmcdf_function_declarations) == nothing, plogProb.input)
    plogProb_declarations = string([a[2] for a in plogProb_function_declarations]...)
    plogProbCond_function_declarations = filter(x->x[4]=="function declaration" && findfirst(y->y[1]==x[1], pmcdf_function_declarations) == nothing && findfirst(y->y[1]==x[1], plogProb_function_declarations) == nothing, plogProbCond.input)
    plogProbCond_declarations = string([a[2] for a in plogProbCond_function_declarations]...)


    pmcdf_functions = filter(x->x[4]=="function", pmcdf.input)
    pmcdf_prg = string([a[2] for a in pmcdf_functions]...)

    plogProb_functions = filter(x->x[4]=="function" && findfirst(y->y[1]==x[1], pmcdf_functions) == nothing, plogProb.input)
    plogProb_prg = string([a[2] for a in plogProb_functions]...)

    plogProbCond_functions = filter(x->x[4]=="function" && findfirst(y->y[1]==x[1], pmcdf_functions) == nothing && findfirst(y->y[1]==x[1], plogProb_functions) == nothing, plogProbCond.input)
    println(plogProbCond_functions)
    plogProbCond_prg = string([a[2] for a in plogProbCond_functions]...)
    # pmcdf.expression
    # plogProb.expression

    mcdfObjName = filter(x->x[3]==-2, pmcdf.input)[1][2]
    logProbObjName = filter(x->x[3]==-2, plogProb.input)[1][2]
    logProbCondObjName = filter(x->x[3]==-2, plogProbCond.input)[1][2]

    extraInput = string(",
   __global float *workmemory,
   __global float *",mcdfObjName,",
   __global float *",logProbObjName)
   globalIndicator = "__global "
   objDeclarations = ""
   mcdfObjOffset = "+gid*mcdfObjSize"
   logProbObjOffset = "+gid*nRepeat"
    if localMem
        extraInput = ""
        globalIndicator = ""
        objDeclarations = string("float workmemory[",nPeakAssume*(2*nDim+2*nParamSinglePeak+2),"];
                                    float ",mcdfObjName,"[",nDim*(nDim+nParamSinglePeak*nPeakAssume),"];
                                    float ",logProbObjName,"[",nPeakAssume,"];\n")
        mcdfObjOffset = ""
        logProbObjOffset = ""
    end

    Q_opencl = string("
        float logLikelihood(global float* logprob, int nRepeat);
        float logLikelihood(global float* logprob, int nRepeat){
            float maxprob = logprob[0];
            float loglikelihood = 0.0;
            int maxindex = 0;
            for(int i=1; i<nRepeat; i++) {
                if(logprob[i] > maxprob) {
                        maxindex = i;
                        maxprob = logprob[i];
                    }
            }

            if (maxprob < -90.0f)
                {
                    loglikelihood= maxprob;
                    for(int i=0; i<nRepeat; i++) {
                        if(i!=maxindex)
                            {
                                loglikelihood += exp(logprob[i] - maxprob);
                            }
                        }
                } else {
                    float prob = 0.0;
                    for(int i=0; i!= nRepeat; i++)
                        {
                          prob += exp( (float)( logprob[i] ) );
                        }
                    loglikelihood+=log((float)(prob));
                }

            return loglikelihood;
        }

        __kernel void Q(__global float *data,
                        __global float *p,
                        __global float *Dp,
                        __global float *sumlogp,
                         int nRepeat,
                         bool ghost
                         ",extraInput,") {\n",
                        " const int nDim=",nDim,";\n",
                        " const int nParam=", nParamSinglePeak, "*nRepeat;\n",

                            "const int gid = get_global_id(0);
                            __global float* x = &data[gid*nDim];

                            ",objDeclarations,"
                            global float *",logProbCondObjName," = ",logProbObjName,";
                            if(ghost)
                            {",pmcdf.definitions,"}\n
                            ",plogProb.definitions,"\n




                            const int mcdfObjSize = ",nDim,"*(",nDim," + nParam);
                            ",globalIndicator," float * dFdx = ",mcdfObjName,mcdfObjOffset,";
                            ",globalIndicator," float * dFdp = dFdx + nDim;
                            ",globalIndicator," float * logprob = ",logProbObjName,logProbObjOffset,";


                            float llLast = logLikelihood(logprob, nRepeat);
                            ",plogProbCond.definitions,"
                            float llCond = logLikelihood(logprob, nRepeat);
/*                            float amplNorm = 0.0f;
                            for(int i=0; i<nRepeat;i++)
                                   {amplNorm+=p[i*",nParamSinglePeak,"] ;}
                            amplNorm = log((float)(amplNorm));
                            for(int i=0; i<nRepeat; i++)
                                {logprob[i] -= amplNorm;}
*/

                            //float y[",nDim,"] = {0.0f};
                            float pert=0.0f;//  [",nDim,"] = {0.0f};
                            float q = 0.0;

                            float lh = llLast - llCond;

                            if(ghost){
                                for(int i=0;i!=nParam;i++) {
                                        float b = fabs((float)(dFdp[ i ]))*exp((float)(Dp[i]-lh));
                                        if(b<1e-40f){
                                            pert += 1e5f;
                                            }
                                        else {
                                                pert -= log(1.0f - b);
                                                //pert+=b;
                                            }

                                    }
                                }

                                /*q = exp((float)( lh )) - pert;

                                if (q <= 0) {
                                        q = 9.9e36;
                                    } else {
                                        q = -log((float)(q) );
                                    }*/

                                  q = -lh + pert;//*exp(-lh);
//                                    q = -lh - log((float)(1.0f - pert) );


//                                for(int i=0; i!= nParam; i++)
//                                    {q-=Dp[i];}

                            if (isinf(q) ){ q=9.9e36; }
                            sumlogp[gid] = q;
                        }



      ")

      #pragmaString = "#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable\n#pragma OPENCL EXTENSION cl_khr_fp64 : disable\n"

      # device, ctx, queue = OpenCL.cl.create_compute_context()
      # prg = OpenCL.cl.Program(ctx, source = string(pmcdf_declarations, plogProb_declarations, pmcdf_prg, plogProb_prg, Q_opencl) ) |> OpenCL.cl.build!

      return string(pmcdf_declarations, plogProb_declarations, plogProbCond_declarations, pmcdf_prg, plogProb_prg, plogProbCond_prg, Q_opencl)
end

function generateRangesQParameteric(input::Vector{Tuple{String,String,Int, String, Union{Nothing,String}}}, nRepeat::Int, variableName::String = "p")
        nDim = maximum([i[3] for i in filter(x->x[2]=="x",input)])
        nParam = maximum([i[3] for i in filter(x->x[2]==variableName,input)])
        nPeak = nRepeat

        params = Float64[]
        paramsMin = Float64[]
        paramsMax = Float64[]

        ranges = Dict{String, Tuple{Float64, Float64}}()
        ranges["mean"] = (-5.0, 5.0)
        ranges["affinePolynome"] = (-10.0, 10.0)
        ranges["sigma"] = (1e-7, 3.0)
        ranges["smear"] = (0.01, 1.0) #4
        ranges["width"]= (1e-5, 3.0) #3
        ranges["amplitude"] = (1e-5, 1.0)#(1e-5, 1.0)

        for peak in 1:nPeak
            par = zeros(maximum([i[3] for i in filter(x->x[2]==variableName,  input) ]))
            parMin = zero(par)
            parMax = zero(par)
            for p in filter(x->x[2]==variableName,  input)
                    mn = ranges[p[4]][1]
                    mx = ranges[p[4]][2]
                    par[p[3]] = rand()*(mx - mn) + mn
                    parMin[p[3]] = mn
                    parMax[p[3]] = mx
            end

            for p in filter(x->x[4] == "mean", input)
                 par[p[3]] = par[p[3]]/100.0
            end

            for p in filter(x->x[4] == "smear", input)
                 par[p[3]] = ranges["smear"][2]
            end

            for p in filter(x->x[4] == "sigma", input)
                 par[p[3]] = (1.0 + rand())*1e-2
            end

            for p in filter(x->x[4] == "width", input)
                 par[p[3]] = 1.0 + rand()
            end

            append!(params,par)
            append!(paramsMin, parMin)
            append!(paramsMax, parMax)

        end
        #nWorkMemory=(2*nDim+2*nParam+2)*nPeak*nDim # in case the execution of kernels is not sequential
        return params, paramsMin, paramsMax, ( nDim = nDim, nParam = nParam, nWorkMemory=(2*nDim+2*nParam+2)*nPeak, nRepeat = nPeak  )

end

function testOpenCL_paramAvgQ(dataF32::Array{Float32,2}, kernel::String, params::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1},NamedTuple{(:nDim, :nParam, :nWorkMemory, :nRepeat),NTuple{4,Int64}}};fullParam::Union{Nothing,Vector} = nothing,ghost::Bool = true)

    println("hello")
    println(typeof(params))

    function localQ(p::Vector{Float32})
        par[:] = p[1:nPar]
        dpar[:] = p[nPar+1:2*nPar]
        println(par)
        println(dpar)
        p_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = par[:])
        dp_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = dpar[:])

        #k[queue, (nData,1)]( d_buff, p_buff, dp_buff, o_buff )
        k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak, ghost, wm_buff, mcdf_buff, lp_buff)
        #k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak)#, wm_buff, mcdf_buff, lp_buff)
        OpenCL.cl.copy!(queue, out, o_buff )
        sl = Float32(sum(out)  - sum(dpar))
        println(sl)
        return sl
    end

    nDim, nData = size(dataF32)
    nPeak = params[4].nRepeat

    nPar = params[4].nParam*nPeak

    println(params[4])
    nWork = ((2*nDim+2)*nPeak + 2*nPar)#*nDim
    nMCDF = nDim*(nDim + nPar)
    println(nWork, "   ", params[4].nWorkMemory)
    nWork = params[4].nWorkMemory
    nOutPut = 1

    par = convert(Vector{Float32}, params[1])
    dpar = fill(-20.0f0,nPar)
    usepar =convert(Vector{Float32},vcat(par,fill(-20.0f0,nPar)))
    if fullParam!=nothing

            usepar = convert(Vector{Float32},fullParam)
            par = convert(Vector{Float32},fullParam[1:nPar])
            dpar = convert(Vector{Float32},fullParam[nPar+1:2*nPar])


    end

    device, ctx, queue = OpenCL.cl.create_compute_context()

    d_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = dataF32[:])
    p_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = par[:])
    dp_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = dpar[:])

    println(par)
    #out = Vector{Float64}(undef, nData)
    wm = Array{Float32,2}(undef, (nWork, nData ))
    wm_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData*nWork)
    mcdf_wm = Array{Float32,2}(undef, (nMCDF, nData ))
    mcdf_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData*nMCDF)
    lp_wm = Array{Float32,2}(undef, (nPeak, nData ))
    lp_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData*nPeak)
    out = Array{Float32,2}(undef, (nOutPut, nData ))
    o_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData*nOutPut)

    prg = OpenCL.cl.Program(ctx, source = kernel ) |> OpenCL.cl.build!

    k = OpenCL.cl.Kernel(prg,"Q")

    k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak, ghost, wm_buff, mcdf_buff, lp_buff)
#    k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak)#, wm_buff, mcdf_buff, lp_buff)

    #queue(k, (nData,1) , nothing, d_buff,p_buff, o_buff, Int32(nGauss) )
    OpenCL.cl.copy!(queue, wm, wm_buff)
    OpenCL.cl.copy!(queue, mcdf_wm, mcdf_buff)
    OpenCL.cl.copy!(queue, lp_wm, lp_buff)
    OpenCL.cl.copy!(queue, out, o_buff)

    # return  Calculus.gradient(p::Vector{Float32}->begin
    #     par[:] = p[1:nPar]
    #     dpar[:] = p[nPar+1:2*nPar]
    #     println(par)
    #     println(dpar)
    #     p_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = par[:])
    #     dp_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = dpar[:])
    #
    #     #k[queue, (nData,1)]( d_buff, p_buff, dp_buff, o_buff )
    #     k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak, wm_buff, mcdf_buff, lp_buff)
    #     #k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak)#, wm_buff, mcdf_buff, lp_buff)
    #     OpenCL.cl.copy!(queue, out, o_buff )
    #     sl = Float32(sum(out)  - sum(dpar))
    #     println(sl)
    #     return sl ;end
    #     ,


    #return  Calculus.gradient(localQ, usepar , :forward )
    return wm,mcdf_wm, lp_wm, out
end

# function convertParams( p1::Vector, p2::Vector, p3::Vector, pcopy::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1},NamedTuple{(:nDim, :nParam, :nWorkMemory, :nRepeat),NTuple{4,Int64}}} )
#     nPar = pcopy[4].nParam*pcopy[4].nRepeat
#     return p1[1:2*nPar], p2[1:2*nPar], p3[1:2*nPar], pcopy[4]
# end

function generateRangesQ(input::Vector{Tuple{String,String,Int, String}}, variableName::String = "p")
        nDim = maximum([i[3] for i in filter(x->x[2]=="x",input)])
        nParam = maximum([i[3] for i in filter(x->x[2]==variableName,input)])

        ffoutput = findfirst(i->i[4] == "output", input)
        if ffoutput!=nothing
            nOutPut = -input[ffoutput][3]
        else
            nOutPut = 0
        end
        nPeak = nOutPut - nParam*nDim-nDim*nDim

        ranges = Dict{String, Tuple{Float64, Float64}}()
        ranges["affinePolynome"] = (-1.0, 1.0)
        ranges["sigma"] = (1e-7, 3.0)
        ranges["smear"] = (0.01, 1.0) #4
        ranges["width"]= (1e-5, 3.0) #3
        ranges["amplitude"] = (1e-5, 1.0)

        par = zeros(maximum([i[3] for i in filter(x->x[2]==variableName,  input) ]))
        parMin = zero(par)
        parMax = zero(par)
        for p in filter(x->x[2]==variableName,  input)
                mn = ranges[p[4]][1]
                mx = ranges[p[4]][2]
                par[p[3]] = rand()*(mx - mn) + mn
                parMin[p[3]] = mn
                parMax[p[3]] = mx
        end

        for p in filter(x->x[4] == "smear", input)
             par[p[3]] = ranges["smear"][2]
        end

        for p in filter(x->x[4] == "sigma", input)
             par[p[3]] = 1.0 + rand()
        end

        for p in filter(x->x[4] == "width", input)
             par[p[3]] = 1.0 + rand()
        end

        return par, parMin, parMax, ( nDim = nDim, nParam = nParam, nOutPut=nOutPut, nPeak = nPeak  )

end

function generateRangesQData(input::Vector{Tuple{String,String,Int, String, Union{Nothing,String}}},
    data::Array{T,2}, nPeak::Int = 8, variableName::String = "p" ) where T <: Real


                nData = size(data)[2]

                nDim = maximum([i[3] for i in filter(x->x[2]=="x",input)])
                nParam = maximum([i[3] for i in filter(x->x[2]==variableName,input)])

                params = Float64[]
                paramsMin = Float64[]
                paramsMax = Float64[]
                ran = String[]

                ranges = Dict{String, Tuple{Float64, Float64}}()
                ranges["mean"] = (-5.0, 5.0)
                ranges["affinePolynome"] = (-10.0, 10.0)
                ranges["sigma"] = (1e-8, 10.0)
                ranges["smear"] = (1e-9, 20.0) #4
                ranges["width"]= (1e-8, 20.0) #3
                ranges["amplitude"] = (1e-10, 1.0)
                ranges["offset"] = (-10.0, 10.0)

                parSort = sort( filter(x->x[2] == variableName, input) , by=x->x[3])
                deviations = std(data,dims=2)
                sampleIndices = sample(1:nData,nPeak, replace = false)


                for i=1:nPeak
                    par = zeros(nParam)
                    parMin = zero(par)
                    parMax = zero(par)
                    j = 1
                    for p in parSort
                        pIndex = p[3]
                        mn = ranges[p[4]][1]
                        mx = ranges[p[4]][2]
                        if p[4] == "mean"
                            par[pIndex] = data[ j,sampleIndices[i]]#/3.0

                            #println(i," mean ",j, " ", par[pIndex])
                        elseif p[4] == "sigma"
                            #par[pIndex] = deviations[j] * log(nDim) #/nPeak^(1.0/nDim)  #/log(nPeak)
                            par[pIndex] = deviations[j]/nPeak^(1.0/nDim)*3.0  #/log(nPeak)
                            #println(i," sigma ",j, " ", par[pIndex])
                            j += 1
                        elseif p[4] == "width"
                            par[pIndex] = deviations[j]
                            #j += 1
                        elseif p[4] == "amplitude"
                            par[pIndex] = 1.0/(nPeak)
                        elseif p[4] == "affinePolynome"
                            par[pIndex] = (rand()*(mx - mn) + mn)*1e-4
                        elseif p[4] == "smear"
                            par[pIndex] = par[pIndex-1]/3.0
                        elseif p[4] == "offset"
                            par[pIndex] *= 1e-4
                        else
                            par[pIndex] = (rand()*(mx - mn) + mn)
                        end
                        parMin[p[3]] = mn
                        parMax[p[3]] = mx
                    end

                    append!(params, par)
                    append!(paramsMin, parMin)
                    append!(paramsMax, parMax)
                end
                dpar = fill(-10.0, length(params))
                dparmin = fill(-99.0, length(params))
                dparmax = fill(0.0, length(params))
                return vcat(params,dpar), vcat(paramsMin,dparmin), vcat(paramsMax, dparmax), (nDim = nDim, nParam = nParam, nWorkMemory=(2*nDim+2*nParam+2)*nPeak, nRepeat = nPeak)
end

function testBoundaries(startVector::Vector{Float64}, lowerBounds::Vector{Float64},upperBounds::Vector{Float64})
    sum([ (s<=u) && (s>=l) for (s,l,u) in zip(startVector, lowerBounds, upperBounds)])
end

function QOptF(inData::Array{Float32, 2},
                    params::Tuple{
                        Vector{Float64},
                        Vector{Float64},
                        Vector{Float64},
                        NamedTuple{(:nDim, :nParam, :nWorkMemory, :nRepeat),NTuple{4,Int64}}
                    },
                    kern,
                    callLimits = [1000,1000])

        fTolerance = 1e-20
        callCounter = 0

        nIData = size(inData)[2]
        nDim = params[4].nDim
        nPeak = params[4].nRepeat
        nParam = params[4].nParam*nPeak
        nWork = params[4].nWorkMemory
        nMCDF = nDim*(nDim + nParam)
        parStart = params[1]
        parMin = params[2]
        parMax = params[3]

        println(params[4])
        #mcdf_func_opencl = kern.mcdf_func_opencl
        Q_opencl = kern

        device, ctx, queue = OpenCL.cl.create_compute_context()
        prg = OpenCL.cl.Program(ctx, source = string(Q_opencl) ) |> OpenCL.cl.build!
        k = OpenCL.cl.Kernel(prg,"Q")

        par = Vector{Float32}(undef, nParam)
        par2 = Vector{Float32}(undef, 2*nParam)

        dPar = Vector{Float32}(undef, nParam)

        errorMsg = ""
        function OptimizeMe(algo::Symbol,
                    startVector::Vector{Float64},
                    lowerBounds::Vector{Float64},
                    upperBounds::Vector{Float64},
                    fTolerance::Float64, callLimit::Int64, verbose = false)

            localCallCounter = 0

            d_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = data[:])
            out = Array{Float32,1}(undef, ( nData ))
            o_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData)
            wm_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData*nWork)
            mcdf_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData*nMCDF)
            lp_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData*nPeak)


            function localQ(p::Vector{Float32})
                par[:] = p[1:nParam]
                dPar[:] = p[nParam+1:2*nParam]

                p_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = par[:])
                dp_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = dPar[:])

                #k[queue, (nData,1)]( d_buff, p_buff, dp_buff, o_buff )
                k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak, wm_buff, mcdf_buff, lp_buff)
                #k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak)#, wm_buff, mcdf_buff, lp_buff)
                OpenCL.cl.copy!(queue, out, o_buff )
                sl = Float32(sum(out)  - sum(dPar))
                return sl
            end
            prevParams = startVector
            function localQGrad(p::Vector{Float64}, g::Vector{Float64})
                #println(p)

                callCounter += 1
                localCallCounter += 1
                nGrad = length(g)
                par[:] = p[1:nParam]
                dPar[:] = p[nParam+1:2*nParam]

                if nGrad>0

                    grad = Calculus.gradient(localQ,convert(Vector{Float32}, p), :forward)
                    for i in 1:nGrad
                        g[i] = grad[i]
                    end
                    #g .= localQ(p)

                    #println(g)
                end
                #setParameters(dgTry, p[1:dgAParamSize])

                #sl = dgQ(p[1:dgAParamSize], p[dgAParamSize+1:end], probHandling, nDim, nGauss, nProbParams, dgAParamSize,  data)




                    # println("it's a me")
                    # println(par)
                    # println(dPar)
                p_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = par[:])
                dp_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = dPar[:])

                #k[queue, (nData,1)]( d_buff, p_buff, dp_buff, o_buff )
                k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak, wm_buff, mcdf_buff, lp_buff)
                #k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak)#, wm_buff, mcdf_buff, lp_buff)
                OpenCL.cl.copy!(queue, out, o_buff )
                sl = Float64(sum(out) - sum(dPar))
                #println(sl)
                if mod(callCounter, 500) == 0
                    println("Call  #",callCounter, " local Call #", localCallCounter , " sl = ", sl)
                end
                #print(sl, "   ")

                if isnan(sl) || sum(isnan.(p))>0 ||  sum(isnan.(g))>0
                    errorMsg = string(errorMsg, "***** sl == NaN *****\n",
                        "Call  #",callCounter, " sl = ", sl,"\n",
                        "prevParams:" , prevParams,"\n",
                        "p: ", p,"\n",
                        "g: ", g,"\n")
                    #error("NaN")

                    #error("** error **")
                end
                prevParams = p
                #println(sl)
                return sl
            end

            if testBoundaries(startVector,lowerBounds,upperBounds) != length(startVector)
                println("boundary test failed!")
                return (lowerBounds, upperBounds, startVector, 9e99, startVector, :FORCED_STOP )
            end

            # mygr = zero(startVector)
            # println(localQGrad(startVector,mygr) )
            # println(mygr)

            println("datasize = " , size(data))
            optLocal = Opt(algo, length(startVector))
            lower_bounds!(optLocal, lowerBounds)
            upper_bounds!(optLocal, upperBounds)
            min_objective!(optLocal, localQGrad)
            #ftol_abs!(optLocal,fTolerance)
            ftol_rel!(optLocal,fTolerance)
            #xtol_rel!(optLocal,fTolerance)
            maxeval!(optLocal, callLimit)
            (minf, minx, ret) = optimize(optLocal,startVector)

            if verbose
                println(algo)
                println("startValues")
                println(startVector)
                println("lowerBounds")
                println(lowerBounds)
                println("upperBounds")
                println(upperBounds)
                println("minf")
                println(minf)
                println("minx")
                println(minx)
                println("ret")
                println(ret)
            end
           return (lowerBounds, upperBounds, startVector, minf, minx, ret)
        end

        DlowerBounds = similar(parStart)
        DupperBounds = similar(parStart)
        startDxValues = similar(parStart)

        DlowerBounds[1:end] .= -99.0#smallDx
        DupperBounds[1:end] .= 0. #largeDx
        startDxValues[1:end] .= -20.

        lBoundsFull = vcat(parMin, DlowerBounds)
        uBoundsFull = vcat(parMax, DupperBounds)
        startValues = vcat(parStart, startDxValues)

        subSample = sample(1:nIData, div(nIData,10),replace = false)
        data = inData[:,subSample]
        nData = size(data)[2]

        println(startValues)
        (lv1, uv1, sv1, minf1, minx1, ret1) = OptimizeMe(:LN_SBPLX, startValues, lBoundsFull, uBoundsFull, fTolerance, callLimits[1])

        startValues2 = deepcopy(minx1)
        lBoundsFull2 = deepcopy(lBoundsFull)
        uBoundsFull2 = deepcopy(uBoundsFull)

        # setParameters(dgTry, startValues2[1:dgAParamSize])
        # #dgSumLogConservativeNoColons(dgTry, startValues2[dgAParamSize+1:end], data, true, true)
        # dgSumLogCDFTotalDxVary(dgTry, startValues2[dgAParamSize+1:end], data, true, false, varyDx)

        println(ret1)
        if length(callLimits) == 1
            return (lBoundsFull, uBoundsFull, startValues2, minf1, vcat(minx1,startDxValues), ret1)
        end


        data = inData
        nData = size(data)[2]

        # if dScale > 1.0
        #     reScale = (0.95-dScale)/nDim/nData
        #     #startValues2[dgAParamSize+1:dgAParamSize*2] .-= log(dScale)/nDim -0.05
        #     startValues2[dgAParamSize+1:dgAParamSize*2].+= reScale
        # end
        if callLimits[2]!=0
            (lv2, uv2, sv2, minf2, minx2, ret2) = OptimizeMe(:LN_SBPLX, startValues2, lBoundsFull2, uBoundsFull2, fTolerance, callLimits[2])
        else
            (lv2, uv2, sv2, minf2, minx2, ret2) = (lBoundsFull, uBoundsFull, startValues2, minf1, minx1, ret1)
        end

        println(ret2)
        if length(callLimits) == 2
        #        println(ret2)
            return (lBoundsFull, uBoundsFull, startValues2, minf2, minx2, ret2)
        end

        startValues3 = deepcopy(minx2)
        lBoundsFull3 = deepcopy(lBoundsFull)
        uBoundsFull3 = deepcopy(uBoundsFull)

        (lv3, uv3, sv3, minf3, minx3, ret3) = OptimizeMe(:LN_COBYLA, startValues3, lBoundsFull3, uBoundsFull3, fTolerance, callLimits[3])

        println(ret3)
        if length(callLimits) == 3
        #        println(ret3)
            return (lBoundsFull, uBoundsFull, minx2, minf3, minx3, ret3)
        end

        startValues4 = deepcopy(minx3)
        lBoundsFull4 = deepcopy(lBoundsFull)
        uBoundsFull4 = deepcopy(uBoundsFull)

        (lv4, uv4, sv4,minf4, minx4, ret4) = OptimizeMe(:LN_SBPLX, startValues4, lBoundsFull4, uBoundsFull4, fTolerance, callLimits[4] )

        println(ret4)
        if length(callLimits) == 4
        #        println(ret3)
            return (lBoundsFull, uBoundsFull, minx3, minf4, minx4, ret4)

        end

        startValues5 = deepcopy(minx4)
        lBoundsFull5 = deepcopy(lBoundsFull)
        uBoundsFull5 = deepcopy(uBoundsFull)
        #LD_MMA
        (lv5, uv5, sv5,minf5, minx5, ret5) = OptimizeMe(:LD_MMA, startValues5, lBoundsFull5, uBoundsFull5, fTolerance, callLimits[5] )

        println(ret5)
        if true #length(callLimits) == 5
        #        println(ret3)
            return (lBoundsFull, uBoundsFull, minx4, minf5, minx5, ret5)

        end


end

function QOpt(inData::Array{Float32, 2},
                    params::Tuple{
                        Vector{Float64},
                        Vector{Float64},
                        Vector{Float64},
                        NamedTuple{(:nDim, :nParam, :nWorkMemory, :nRepeat),NTuple{4,Int64}}
                    },
                    kern,
                    callLimits = [1000,1000])

        fTolerance = 1e-10
        callCounter = 0

        nIData = size(inData)[2]
        nDim = params[4].nDim
        nPeak = params[4].nRepeat
        nParam = params[4].nParam*nPeak
        nWork = params[4].nWorkMemory
        nMCDF = nDim*(nDim + nParam)
        parStart = params[1]
        parMin = params[2]
        parMax = params[3]

        println(params[4])
        #mcdf_func_opencl = kern.mcdf_func_opencl
        Q_opencl = kern

        # device = OpenCL.cl.devices()[1]
        # ctx = OpenCL.cl.Context(device)
        # queue = OpenCL.cl.CmdQueue(ctx, device)

        prg = OpenCL.cl.Program(ctx, source = string(Q_opencl) ) |> OpenCL.cl.build!
        k = OpenCL.cl.Kernel(prg,"Q")

        par = Vector{Float32}(undef, nParam)
        par2 = Vector{Float32}(undef, 2*nParam)

        dPar = Vector{Float32}(undef, nParam)

        errorMsg = ""
        function OptimizeMe(algo::Symbol,
                    startVector::Vector{Float64},
                    lowerBounds::Vector{Float64},
                    upperBounds::Vector{Float64},
                    fTolerance::Float64, callLimit::Int64, ghost = true, verbose = false)

            localCallCounter = 0

            d_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = data[:])
            out = Array{Float32,1}(undef, ( nData ))
            o_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData)
            wm_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData*nWork)
            mcdf_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData*nMCDF)
            lp_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData*nPeak)


            function localQ(p::Vector{Float32})
                par[:] = p[1:nParam]
                if ghost
                    dPar[:] = p[nParam+1:2*nParam]
                end

                p_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = par[:])
                dp_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = dPar[:])

                #k[queue, (nData,1)]( d_buff, p_buff, dp_buff, o_buff )
                k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak, ghost, wm_buff, mcdf_buff, lp_buff)
                #k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak)#, wm_buff, mcdf_buff, lp_buff)
                OpenCL.cl.copy!(queue, out, o_buff )
                sl = Float32(sum(out)  - sum(dPar))
                return sl
            end
            prevParams = startVector
            function localQGrad(p::Vector{Float64}, g::Vector{Float64})
                #println(p)

                callCounter += 1
                localCallCounter += 1
                nGrad = length(g)
                par[:] = p[1:nParam]
                if ghost
                    dPar[:] = p[nParam+1:2*nParam]
                end

                if nGrad>0

                    grad = Calculus.gradient(localQ,convert(Vector{Float32}, p), :forward)
                    for i in 1:nGrad
                        g[i] = grad[i]
                    end
                    #g .= localQ(p)

                    #println(g)
                end
                #setParameters(dgTry, p[1:dgAParamSize])

                #sl = dgQ(p[1:dgAParamSize], p[dgAParamSize+1:end], probHandling, nDim, nGauss, nProbParams, dgAParamSize,  data)




                    # println("it's a me")
                    # println(par)
                    # println(dPar)
                p_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = par[:])
                dp_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = dPar[:])

                #k[queue, (nData,1)]( d_buff, p_buff, dp_buff, o_buff )
                k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak, ghost, wm_buff, mcdf_buff, lp_buff)
                #k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak)#, wm_buff, mcdf_buff, lp_buff)
                OpenCL.cl.copy!(queue, out, o_buff )
                sl = Float64(sum(out) - sum(dPar))
                #println(sl)
                if mod(callCounter, 500) == 0
                    println("Call  #",callCounter, " local Call #", localCallCounter , " sl = ", sl)
                end
                #print(sl, "   ")

                if isnan(sl) || sum(isnan.(p))>0 ||  sum(isnan.(g))>0
                    errorMsg = string(errorMsg, "***** sl == NaN *****\n",
                        "Call  #",callCounter, " sl = ", sl,"\n",
                        "prevParams:" , prevParams,"\n",
                        "p: ", p,"\n",
                        "g: ", g,"\n")
                    #error("NaN")

                    #error("** error **")
                end
                prevParams = p
                #println(sl)
                return sl
            end

            if testBoundaries(startVector,lowerBounds,upperBounds) != length(startVector)
                println("boundary test failed!")
                return (lowerBounds, upperBounds, startVector, 9e99, startVector, :FORCED_STOP )
            end

            # mygr = zero(startVector)
            # println(localQGrad(startVector,mygr) )
            # println(mygr)

            println("datasize = " , size(data))
            optLocal = Opt(algo, length(startVector))
            lower_bounds!(optLocal, lowerBounds)
            upper_bounds!(optLocal, upperBounds)
            min_objective!(optLocal, localQGrad)
            #ftol_abs!(optLocal,fTolerance)
            ftol_rel!(optLocal,fTolerance)
            #xtol_rel!(optLocal,fTolerance)
            maxeval!(optLocal, callLimit)
            (minf, minx, ret) = optimize(optLocal,startVector)

            if verbose
                println(algo)
                println("startValues")
                println(startVector)
                println("lowerBounds")
                println(lowerBounds)
                println("upperBounds")
                println(upperBounds)
                println("minf")
                println(minf)
                println("minx")
                println(minx)
                println("ret")
                println(ret)
            end
           return (lowerBounds, upperBounds, startVector, minf, minx, ret)
        end


        lBoundsFull = parMin
        uBoundsFull = parMax
        startValues = parStart

        subSample = sample(1:nIData, div(nIData,10),replace = false)
        if length(callLimits) > 1
            data = inData[:,subSample]
        else
            data = inData
        end
        nData = size(data)[2]

        ghost1 = true
        ghostf = 1+ghost1

        dPar[:] = startValues[nParam+1:2*nParam]
        #@show startValues
        (lv1, uv1, sv1, minf1, minx1, ret1) = OptimizeMe(:LN_SBPLX, startValues[1:nParam*ghostf], lBoundsFull[1:nParam*ghostf], uBoundsFull[1:nParam*ghostf], 1e-5, callLimits[1], ghost1)

        if !ghost1
            startValues2 = deepcopy(vcat(minx1, dPar) )
        else
            startValues2 = deepcopy(minx1)
        end
        lBoundsFull2 = deepcopy(lBoundsFull)
        uBoundsFull2 = deepcopy(uBoundsFull)

        # setParameters(dgTry, startValues2[1:dgAParamSize])
        # #dgSumLogConservativeNoColons(dgTry, startValues2[dgAParamSize+1:end], data, true, true)
        # dgSumLogCDFTotalDxVary(dgTry, startValues2[dgAParamSize+1:end], data, true, false, varyDx)

        println(ret1)
        if length(callLimits) == 1
            return  (startValues2,lBoundsFull,uBoundsFull,params[4]),minf1,ret1  #(lBoundsFull, uBoundsFull, startValues, minf1, startValues2, ret1)
        end


        subSample = sample(1:nIData, div(nIData,3)*2,replace = false)
        data = inData[:,subSample]
        nData = size(data)[2]



        # if dScale > 1.0
        #     reScale = (0.95-dScale)/nDim/nData
        #     #startValues2[dgAParamSize+1:dgAParamSize*2] .-= log(dScale)/nDim -0.05
        #     startValues2[dgAParamSize+1:dgAParamSize*2].+= reScale

        # end
        ghost2 = true
        ghostf = 1+ghost2
        if callLimits[2]!=0
            (lv2, uv2, sv2, minf2, minx2, ret2) = OptimizeMe(:LN_SBPLX, startValues2[1:nParam*ghostf], lBoundsFull2[1:nParam*ghostf], uBoundsFull2[1:nParam*ghostf], 1e-7, callLimits[2],ghost2)
        else
            (lv2, uv2, sv2, minf2, minx2, ret2) = (lBoundsFull, uBoundsFull, startValues2, minf1, minx1, ret1)
        end

        println(ret2)
        if !ghost2
            startValues3 = deepcopy(vcat(minx2, dPar))
        else
            startValues3 = deepcopy(minx2)
        end
        lBoundsFull3 = deepcopy(lBoundsFull)
        uBoundsFull3 = deepcopy(uBoundsFull)


        if length(callLimits) == 2
        #        println(ret2)
            return (startValues3,lBoundsFull,uBoundsFull,params[4]),minf2,ret2 #(lBoundsFull, uBoundsFull, startValues2, minf2, startValues3, ret2)
        end


        data = inData
        nData = size(data)[2]

        (lv3, uv3, sv3, minf3, minx3, ret3) = OptimizeMe(:LN_SBPLX, startValues3, lBoundsFull3, uBoundsFull3, fTolerance, callLimits[3])

        println(ret3)
        if length(callLimits) == 3
        #        println(ret3)
            return   (minx3,lBoundsFull,uBoundsFull,params[4]),minf3,ret3 #(lBoundsFull, uBoundsFull, minx2, minf3, minx3, ret3)
        end

        startValues4 = deepcopy(minx3)
        lBoundsFull4 = deepcopy(lBoundsFull)
        uBoundsFull4 = deepcopy(uBoundsFull)

        (lv4, uv4, sv4,minf4, minx4, ret4) = OptimizeMe(:LN_SBPLX, startValues4, lBoundsFull4, uBoundsFull4, fTolerance, callLimits[4] )

        println(ret4)
        if length(callLimits) == 4
        #        println(ret3)
            return (minx4,lBoundsFull,uBoundsFull,params[4]),minf4,ret4 #(lBoundsFull, uBoundsFull, minx3, minf4, minx4, ret4)

        end

        startValues5 = deepcopy(minx4)
        lBoundsFull5 = deepcopy(lBoundsFull)
        uBoundsFull5 = deepcopy(uBoundsFull)
        #LD_MMA
        (lv5, uv5, sv5,minf5, minx5, ret5) = OptimizeMe(:LD_MMA, startValues5, lBoundsFull5, uBoundsFull5, fTolerance, callLimits[5] )

        println(ret5)
        if true #length(callLimits) == 5
        #        println(ret3)
            return (minx5,lBoundsFull,uBoundsFull,params[4]),minf5,ret5 #(lBoundsFull, uBoundsFull, minx4, minf5, minx5, ret5)

        end


end

function QOptBatch(inData::Array{Float32, 2},
                    params::Tuple{
                        Vector{Float64},
                        Vector{Float64},
                        Vector{Float64},
                        NamedTuple{(:nDim, :nParam, :nWorkMemory, :nRepeat),NTuple{4,Int64}}
                    },
                    kern,
                    callLimits = [1000,1000])

        fTolerance = 1e-10
        callCounter = 0

        nIData = size(inData)[2]
        nDim = params[4].nDim
        nPeak = params[4].nRepeat
        nParam = params[4].nParam*nPeak
        nWork = params[4].nWorkMemory
        nMCDF = nDim*(nDim + nParam)
        parStart = params[1]
        parMin = params[2]
        parMax = params[3]

        println(params[4])
        #mcdf_func_opencl = kern.mcdf_func_opencl
        Q_opencl = kern

        # device = OpenCL.cl.devices()[1]
        # ctx = OpenCL.cl.Context(device)
        # queue = OpenCL.cl.CmdQueue(ctx, device)

        prg = OpenCL.cl.Program(ctx, source = string(Q_opencl) ) |> OpenCL.cl.build!
        k = OpenCL.cl.Kernel(prg,"Q")

        par = Vector{Float32}(undef, nParam)
        par2 = Vector{Float32}(undef, 2*nParam)

        dPar = Vector{Float32}(undef, nParam)

        errorMsg = ""
        function OptimizeMe(algo::Symbol,
                    startVector::Vector{Float64},
                    lowerBounds::Vector{Float64},
                    upperBounds::Vector{Float64},
                    fTolerance::Float64, callLimit::Int64, ghost = true, verbose = false)

            localCallCounter = 0

            d_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = data[:])
            out = Array{Float32,1}(undef, ( nData ))
            o_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData)
            wm_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData*nWork)
            mcdf_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData*nMCDF)
            lp_buff = OpenCL.cl.Buffer(Float32,ctx,(:rw), nData*nPeak)


            function localQ(p::Vector{Float32})
                par[:] = p[1:nParam]
                if ghost
                    dPar[:] = p[nParam+1:2*nParam]
                end

                p_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = par[:])
                dp_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = dPar[:])

                #k[queue, (nData,1)]( d_buff, p_buff, dp_buff, o_buff )
                k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak, ghost, wm_buff, mcdf_buff, lp_buff)
                #k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak)#, wm_buff, mcdf_buff, lp_buff)
                OpenCL.cl.copy!(queue, out, o_buff )
                sl = Float32(sum(out)  - sum(dPar))
                return sl
            end
            prevParams = startVector
            function localQGrad(p::Vector{Float64}, g::Vector{Float64})
                #println(p)

                callCounter += 1
                localCallCounter += 1
                nGrad = length(g)
                par[:] = p[1:nParam]
                if ghost
                    dPar[:] = p[nParam+1:2*nParam]
                end

                if nGrad>0

                    grad = Calculus.gradient(localQ,convert(Vector{Float32}, p), :forward)
                    for i in 1:nGrad
                        g[i] = grad[i]
                    end
                    #g .= localQ(p)

                    #println(g)
                end
                #setParameters(dgTry, p[1:dgAParamSize])

                #sl = dgQ(p[1:dgAParamSize], p[dgAParamSize+1:end], probHandling, nDim, nGauss, nProbParams, dgAParamSize,  data)




                    # println("it's a me")
                    # println(par)
                    # println(dPar)
                p_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = par[:])
                dp_buff = OpenCL.cl.Buffer(Float32,ctx,(:r, :copy), hostbuf = dPar[:])

                #k[queue, (nData,1)]( d_buff, p_buff, dp_buff, o_buff )
                k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak, ghost, wm_buff, mcdf_buff, lp_buff)
                #k[queue, (nData,1) ](d_buff,p_buff, dp_buff, o_buff, nPeak)#, wm_buff, mcdf_buff, lp_buff)
                OpenCL.cl.copy!(queue, out, o_buff )
                sl = Float64(sum(out) - sum(dPar))
                #println(sl)
                if mod(callCounter, 500) == 0
                    println("Call  #",callCounter, " local Call #", localCallCounter , " sl = ", sl)
                end
                #print(sl, "   ")

                if isnan(sl) || sum(isnan.(p))>0 ||  sum(isnan.(g))>0
                    errorMsg = string(errorMsg, "***** sl == NaN *****\n",
                        "Call  #",callCounter, " sl = ", sl,"\n",
                        "prevParams:" , prevParams,"\n",
                        "p: ", p,"\n",
                        "g: ", g,"\n")
                    #error("NaN")

                    #error("** error **")
                end
                prevParams = p
                #println(sl)
                return sl
            end

            if testBoundaries(startVector,lowerBounds,upperBounds) != length(startVector)
                println("boundary test failed!")
                return (lowerBounds, upperBounds, startVector, 9e99, startVector, :FORCED_STOP )
            end

            # mygr = zero(startVector)
            # println(localQGrad(startVector,mygr) )
            # println(mygr)

            println("datasize = " , size(data))
            optLocal = Opt(algo, length(startVector))
            lower_bounds!(optLocal, lowerBounds)
            upper_bounds!(optLocal, upperBounds)
            min_objective!(optLocal, localQGrad)
            #ftol_abs!(optLocal,fTolerance)
            ftol_rel!(optLocal,fTolerance)
            #xtol_rel!(optLocal,fTolerance)
            maxeval!(optLocal, callLimit)
            (minf, minx, ret) = optimize(optLocal,startVector)

            if verbose
                println(algo)
                println("startValues")
                println(startVector)
                println("lowerBounds")
                println(lowerBounds)
                println("upperBounds")
                println(upperBounds)
                println("minf")
                println(minf)
                println("minx")
                println(minx)
                println("ret")
                println(ret)
            end
           return (lowerBounds, upperBounds, startVector, minf, minx, ret)
        end


        lBoundsFull = parMin
        uBoundsFull = parMax
        startValues = parStart

        ghost1 = true
        ghostf = 1+ghost1

        shuffle = sample(1:nIData, nIData,replace = false)
        shData = copy(inData[:,shuffle])
        #shData = copy(inData)
        minf1 = 1e99
        ret1 = :FORCED_STOP
        for segm in 500:500:nIData
            data = shData[:,1:segm]
            nData = size(data)[2]
            dPar[:] = startValues[nParam+1:2*nParam]
            (lv1, uv1, sv1, minf1, minx, ret1) = OptimizeMe(:LN_SBPLX, startValues[1:nParam*ghostf], lBoundsFull[1:nParam*ghostf], uBoundsFull[1:nParam*ghostf], 1e-5, callLimits[1], ghost1)
            if !ghost1
                startValues = deepcopy(vcat(minx, dPar) )
            else
                startValues = deepcopy(minx)
            end
        end


        startValues2 = deepcopy(startValues)
        lBoundsFull2 = deepcopy(lBoundsFull)
        uBoundsFull2 = deepcopy(uBoundsFull)

        # setParameters(dgTry, startValues2[1:dgAParamSize])
        # #dgSumLogConservativeNoColons(dgTry, startValues2[dgAParamSize+1:end], data, true, true)
        # dgSumLogCDFTotalDxVary(dgTry, startValues2[dgAParamSize+1:end], data, true, false, varyDx)

        println(ret1)
        if length(callLimits) == 1
            return  (startValues2,lBoundsFull,uBoundsFull,params[4]),minf1,ret1  #(lBoundsFull, uBoundsFull, startValues, minf1, startValues2, ret1)
        end


        #subSample = sample(1:nIData, div(nIData,3)*2,replace = false)
        data = inData#[:,subSample]
        nData = size(data)[2]



        # if dScale > 1.0
        #     reScale = (0.95-dScale)/nDim/nData
        #     #startValues2[dgAParamSize+1:dgAParamSize*2] .-= log(dScale)/nDim -0.05
        #     startValues2[dgAParamSize+1:dgAParamSize*2].+= reScale

        # end
        ghost2 = true
        ghostf = 1+ghost2
        if callLimits[2]!=0
            (lv2, uv2, sv2, minf2, minx2, ret2) = OptimizeMe(:LN_SBPLX, startValues2[1:nParam*ghostf], lBoundsFull2[1:nParam*ghostf], uBoundsFull2[1:nParam*ghostf], 1e-7, callLimits[2],ghost2)
        else
            (lv2, uv2, sv2, minf2, minx2, ret2) = (lBoundsFull, uBoundsFull, startValues2, minf1, minx1, ret1)
        end

        println(ret2)
        if !ghost2
            startValues3 = deepcopy(vcat(minx2, dPar))
        else
            startValues3 = deepcopy(minx2)
        end
        lBoundsFull3 = deepcopy(lBoundsFull)
        uBoundsFull3 = deepcopy(uBoundsFull)


        if length(callLimits) == 2
        #        println(ret2)
            return (startValues3,lBoundsFull,uBoundsFull,params[4]),minf2,ret2 #(lBoundsFull, uBoundsFull, startValues2, minf2, startValues3, ret2)
        end


        data = inData
        nData = size(data)[2]

        (lv3, uv3, sv3, minf3, minx3, ret3) = OptimizeMe(:LN_SBPLX, startValues3, lBoundsFull3, uBoundsFull3, fTolerance, callLimits[3])

        println(ret3)
        if length(callLimits) == 3
        #        println(ret3)
            return   (minx3,lBoundsFull,uBoundsFull,params[4]),minf3,ret3 #(lBoundsFull, uBoundsFull, minx2, minf3, minx3, ret3)
        end

        startValues4 = deepcopy(minx3)
        lBoundsFull4 = deepcopy(lBoundsFull)
        uBoundsFull4 = deepcopy(uBoundsFull)

        (lv4, uv4, sv4,minf4, minx4, ret4) = OptimizeMe(:LN_SBPLX, startValues4, lBoundsFull4, uBoundsFull4, fTolerance, callLimits[4] )

        println(ret4)
        if length(callLimits) == 4
        #        println(ret3)
            return (minx4,lBoundsFull,uBoundsFull,params[4]),minf4,ret4 #(lBoundsFull, uBoundsFull, minx3, minf4, minx4, ret4)

        end

        startValues5 = deepcopy(minx4)
        lBoundsFull5 = deepcopy(lBoundsFull)
        uBoundsFull5 = deepcopy(uBoundsFull)
        #LD_MMA
        (lv5, uv5, sv5,minf5, minx5, ret5) = OptimizeMe(:LD_MMA, startValues5, lBoundsFull5, uBoundsFull5, fTolerance, callLimits[5] )

        println(ret5)
        if true #length(callLimits) == 5
        #        println(ret3)
            return (minx5,lBoundsFull,uBoundsFull,params[4]),minf5,ret5 #(lBoundsFull, uBoundsFull, minx4, minf5, minx5, ret5)

        end


end

function pruneParams(params::Tuple{
    Vector{Float64},
    Vector{Float64},
    Vector{Float64},
    NamedTuple{(:nDim, :nParam, :nWorkMemory, :nRepeat),NTuple{4,Int64}}})

    #parNames = [s[4] for s=sort(filter(x->x[2] == "p",input), by = x->x[3])]
    nParam = params[4].nParam
    nRepeat = params[4].nRepeat
    nParamSize = nParam*nRepeat
    threshold = params[1][1:nParam:nParamSize] .> 1e-6
    # sigmaIndices = findall(x->x=="sigma", parNames)
    # for (i,p) in enumerate( Iterators.partition(params[1][1:nParamSize],nParam) )
    #     sigmaThresholds = p[sigmaIndices] .> 1e-6
    #     threshold[i] *= prod(sigmaThresholds)
    # end

    nPrunedRepeat = sum(threshold)
    prunedParams = Float64[]
    prunedParamsMin = Float64[]
    prunedParamsMax = Float64[]

    function pushParams( threshold::BitArray{1}, target::Vector{Float64}, source::Vector{Float64} )
        for i in 1:length(threshold)
            if threshold[i]
                append!(target, source[(i-1)*nParam+1:i*nParam])
            end
        end
    end

    pThreshold = vcat(threshold, threshold)
    pushParams( pThreshold, prunedParams, params[1] )
    pushParams( pThreshold, prunedParamsMin, params[2] )
    pushParams( pThreshold, prunedParamsMax, params[3] )

    return prunedParams, prunedParamsMin, prunedParamsMax, (nDim = params[4].nDim, nParam = nParam, nWorkMemory = div(params[4].nWorkMemory, params[4].nRepeat)*nPrunedRepeat, nRepeat = nPrunedRepeat)

end

function pruneParams(params::Tuple{
    Vector{Float64},
    Vector{Float64},
    Vector{Float64},
    NamedTuple{(:nDim, :nParam, :nWorkMemory, :nRepeat),NTuple{4,Int64}}},
    input::Vector{Tuple{String,String,Int, String, Union{Nothing,String}}})

    parNames = [s[4] for s=sort(filter(x->x[2] == "p",input), by = x->x[3])]
    nParam = params[4].nParam
    nRepeat = params[4].nRepeat
    nParamSize = nParam*nRepeat
    threshold = params[1][1:nParam:nParamSize] .> 1e-6
    sigmaIndices = findall(x->x=="sigma", parNames)
    for (i,p) in enumerate( Iterators.partition(params[1][1:nParamSize],nParam) )
        sigmaThresholds = p[sigmaIndices] .> 1e-3
        threshold[i] *= prod(sigmaThresholds)
    end

    nPrunedRepeat = sum(threshold)
    prunedParams = Float64[]
    prunedParamsMin = Float64[]
    prunedParamsMax = Float64[]

    function pushParams( threshold::BitArray{1}, target::Vector{Float64}, source::Vector{Float64} )
        for i in 1:length(threshold)
            if threshold[i]
                append!(target, source[(i-1)*nParam+1:i*nParam])
            end
        end
    end

    pThreshold = vcat(threshold, threshold)
    pushParams( pThreshold, prunedParams, params[1] )
    pushParams( pThreshold, prunedParamsMin, params[2] )
    pushParams( pThreshold, prunedParamsMax, params[3] )

    return prunedParams, prunedParamsMin, prunedParamsMax, (nDim = params[4].nDim, nParam = nParam, nWorkMemory = div(params[4].nWorkMemory, params[4].nRepeat)*nPrunedRepeat, nRepeat = nPrunedRepeat)

end


device, ctx, queue = OpenCL.cl.create_compute_context()
