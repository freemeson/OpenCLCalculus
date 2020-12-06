using Random
Random.seed!(201920261032)
#sample generation
d3 = MixtureModel(MvNormal[
    MultivariateNormal([1.1, 4.0], [1.0 0.0; 0.0 2.0]),
    MultivariateNormal([2.1, 2.0], [3.0 0.0; 0.0 0.4]),
    MultivariateNormal([7.5, 5.5], [1.0 0.0; 0.0 1.0])
    ], [0.1, 0.8, 0.1])

z3 = convert(Array{Float32,2},rand( d3, 10000))./10.0f0

#fit procedure
pmcdf, plogProb, plast = vector_marginalCDF_opencl(2)
writeFunctionFile(plast,"testGauss", "test.jl")
    include("test.jl")
kernel = createQKernel(pmcdf,plogProb)
params = generateRangesQData(pmcdf.input,z3[:,1:2000], 8)
gaussFit = QOpt(z3[:,1:2000], params , kernel, [15000])

#results
contour( -1.4:0.005:1.4, -0.2:0.01:1.0, (x,y)->begin b=log(testGauss([x,y],gaussFit[1][1],gaussFit[1][4].nRepeat))
                                                        if b<-30.0
                                                                b=-30.0
                                                        elseif b>5.0
                                                            b=5.0
                                                        end
                                                            return b;
                                                        end, fill = true)
scatter!(z3[1,1:2000], z3[2,1:2000], markersize = 2)
