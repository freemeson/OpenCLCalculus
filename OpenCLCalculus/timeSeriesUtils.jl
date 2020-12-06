using Mmap
using CodecZlib
using DataFrames
#using Pipe
using CSV
using Dates

function readCSV(filename = "krakenEUR.csv.gz")
    csvFile = CSV.File(transcode(GzipDecompressor, Mmap.mmap(filename))) |> DataFrame
    df = rename(csvFile,1=>:UnixTime, 2=>:Price, 3=>:Volume)
    df.year = year.(unix2datetime.(df[:,1]))
    df.day = dayofyear.( unix2datetime.(df[:,1]) )
    df.hour = hour.( unix2datetime.(df[:,1]) )
    df.divMinute = div.(second.(unix2datetime.(df[:,1])) + 60*minute.(unix2datetime.(df[:,1])),450)
    df.diffPrice = vcat(0,diff(df.Price))
    df.floorMinute = datetime2unix.(floor.(unix2datetime.(df.UnixTime), Dates.Second( 450 ))) #in correspondence with divMinute
    return df

    dfA = by(df, [:year,:day,:hour]) do df
        m = dot(df.Price,df.Volume)
        v = sum(df.Volume)
        md = dot(df.diffPrice,df.Volume)
        DataFrame(Price = m/v, Volume = v, diffPrice = md/v )
    end

    dfClose = filter(row->row.hour == 23,dfA )
    dfB = by(dfA,[:year,:day]) do df

        dayNorms = filter(row->row.hour == 0, df)
        if length(dayNorms.Price) !=0
            p = df.Price / dayNorms.Price[1]/10.0
            d = df.diffPrice / dayNorms.diffPrice[1]/10.0
            v = df.Volume /  dayNorms.Volume[1]/10.0
            dh = diff(p)
        #p = df.Price /  maximum(df.Price)
            return DataFrame(Price = p, Volume = v, hour = df.hour/24.0, diffPrice = d, dailyDiffPrice = vcat(0,dh))
        else
            p = df.Price/df.Price[1]
            dh = diff(p)

            return DataFrame(Price = p, Volume =df.Volume/df.Volume[1]  ,
                    hour = df.hour/24.0, diffPrice = df.diffPrice/df.diffPrice[1], dailyDiffPrice = vcat(0,dh) )
        end
    end

    return dfB
    aggregatedHours = collect(Matrix(filter(row->row.hour != 0, dfB[[:hour,:Price, :Volume, :diffPrice, :dailyDiffPrice]]))')
    aggregatedHours[1,:] .+= rand(size(aggregatedHours)[2])/24.0 .+ 1/48.0
    dayPrice = collect(hcat(dfClose.year.-2014.0 .+ dfClose.day/365.0,dfClose.Price./maximum(dfClose.Price))')
    scatter(df.day,df.Price)
end

function seriesProbability(params::Vector{Float64},
                    series::Vector{Float32},
                    x::Float64,y::Float64,nDim = 2, nRepeat = 2)
    xInt = Int64(floor(x))
    Coord = vcat(series[xInt-nDim+1:xInt-1],y)

    #println(dgMarginalProbability(dgA,dgCoord)[nDim])
    return testJLAvg(Coord, params, nRepeat)
end



function plotDiffPrediction(params::Vector{Float64},series::Vector{Float32}, nDim::Int64 = 2, nRepeat = 2, Range = -0.02:0.0001:0.02)
    dgAParamSize = div(length(params),2)

    nSize = length(series)
    #
    #temp = [log(seriesProbability(dgA,series,Float64(x),Float64(y),nDim)) for y=-30:0.1:30,x=nDim:nSize-1 ]
    #heatmap(temp)
    pl = heatmap(nDim:nSize-1, Range,(x,y)->(seriesProbability(params,series,Float64(x),Float64(y),nDim, nRepeat)))
    plot!(pl,nDim:nSize-1,series[nDim:nSize-1],linecolor = :blue, width = 3)
end

function seriesProbabilityOffset(params::Vector{Float64},
                    series::Vector{Float32},
                    x::Float64,y::Float64,nDim = 2, nRepeat = 2)
    xInt = Int64(floor(x))
    Coord = vcat(series[xInt-nDim+1:xInt-1] .- series[xInt-nDim] ,y- series[xInt-nDim])

    #println(Coord)
    #println(dgMarginalProbability(dgA,dgCoord)[nDim])
    return testJLAvg(Coord, params, nRepeat)
end


function plotDiffPredictionOffset(params::Vector{Float64},series::Vector{Float32}, nDim::Int64 = 2, nRepeat = 2, Range = -0.02:0.0001:0.02)
    dgAParamSize = div(length(params),2)

    nSize = length(series)
    #
    #temp = [log(seriesProbability(dgA,series,Float64(x),Float64(y),nDim)) for y=-30:0.1:30,x=nDim:nSize-1 ]
    #heatmap(temp)
    offSeries = series[nDim+1:nSize] #.- series[1:nSize-nDim]
    ex = extrema(offSeries)
    extra = (ex[2]-ex[1]) * 0.1
    Range = range(ex[1]-extra, ex[2]+extra, length = 40)
    pl = heatmap(nDim+1:nSize, Range,(x,y)->log(seriesProbabilityOffset(params,series,Float64(x),Float64(y),nDim, nRepeat)))

    #pl = plot()
    plot!(pl,nDim+1:nSize,offSeries,linecolor = :blue, width = 3)
    scatter!(pl,nDim+1:nSize,offSeries,linecolor = :red, markersize = 3)

    # norm,nDim-1,predictPoint
    # [1:nDim]

end

function plotM(params::Tuple{
    Vector{Float64},
    Vector{Float64},
    Vector{Float64},
    NamedTuple{(:nDim, :nParam, :nWorkMemory, :nRepeat),NTuple{4,Int64}}
},data::Array{Float32,2})

    nSize = size(data)
    ex = extrema(data[end,:])
    extra = (ex[2]-ex[1]) * 0.1
    Range = range(ex[1]-extra, ex[2]+extra, length = 40)

    pl = heatmap(1:nSize[2],  Range, (x,y)->testJLAvg( vcat(data[ 1:end-1, floor.(x)],y),  params[1],params[4].nRepeat ))
    plot!( pl, 1:nSize[2],data[end,:], linecolor = :blue, width = 3, label = false )
    scatter!( pl, 1:nSize[2],data[end,:],linecolor = :red, markersize = 3, label = false )
end

function plot1d( dataPoint::Vector, params::Tuple{
    Vector{Float64},
    Vector{Float64},
    Vector{Float64},
    NamedTuple{(:nDim, :nParam, :nWorkMemory, :nRepeat),NTuple{4,Int64}}
},r::StepRangeLen = -2:0.01:2 )
    function t0d(a,b)
        1.0/( (a-b)^2*800^2+1 ) *800.0
    end
    println(vcat(dataPoint[ 1:end-1 ]))
    plot(r, x->testJLAvg( vcat(dataPoint[ 1:end-1 ],  x ) ,
        params[1],params[4].nRepeat) )
end

function plotMV(params::Tuple{
    Vector{Float64},
    Vector{Float64},
    Vector{Float64},
    NamedTuple{(:nDim, :nParam, :nWorkMemory, :nRepeat),NTuple{4,Int64}}
},data::Array{Float32,2}, offset::Vector, vRange::UnitRange{Int64})
    function t1(a, b)
        (a/b -1.0)*100.0
    end
    function t1inv(a,b)
        (a/100.0+1.0)*b
    end
    cf = 100.0
    function t0(a,b)
        atan((a - b)*cf)
    end
    function t0d(a,b)
        1.0/( (a-b)^2*cf^2+1 ) *cf
    end
    function t0inv(a,b)
        tan(a)/cf + b
    end
    function t3(a,b)
        erf((a - b)*100.0)
    end
    function t3d(a,b)
        exp( -((a-b)*100.0)^2 )*2.0 /sqrt(π)*100.0
    end
    function t3inv(a,b)
        erfinv(a*0.999999)/100.0 + b
    end

    function t2(a,b)
        ((a - b)*100.0)
    end
    function t2d(a,b)
        100.0
    end
    function t2inv(a,b)
        (a)/100.0 + b
    end
    nSize = size(data[:,vRange])
    ex = extrema( t0inv.(data[end,vRange], offset[vRange]) )
    println(ex)
    extra = (ex[2]-ex[1]) *0.1#*0.1 #+0.02
    Range = range(ex[1]-extra, ex[2]+extra, length = 160)

    z = [  (testJLAvg( vcat(data[ 1:end-1, x],t0(y,offset[x]) ) ,
        params[1],params[4].nRepeat)*t0d(y,offset[x]) ) for y=Range,x=vRange ]

    pl = heatmap(vRange,Range, (z))
    avg = [  dot(c,Range)./sum(c) for c=eachcol(z)]
    # pl = heatmap(vRange,  Range,
    #     (x,y)->(testJLAvg( vcat(data[ 1:end-1, x],t0(y,offset[x]) ) ,
    #         params[1],params[4].nRepeat)*t0d(y,offset[x]) ) )

    plot!(pl, vRange, avg , linecolor = :black)
    plot!( pl, vRange,t0inv.(data[end,vRange], offset[vRange]), linecolor = :blue, width = 3, label = false )
    scatter!( pl, vRange,t0inv.(data[end,vRange], offset[vRange]),linecolor = :red, markersize = 3, label = false )
end

function prepareDF(data::DataFrame, symbols::Vector{Symbol} = [:year, :day, :hour, :divMinute])
    combine( groupby(data, symbols) ) do df
        m = dot(df.Price,df.Volume)
        v = sum(df.Volume)
        md = dot(df.diffPrice,df.Volume)
        sp = std(df.Price)
        if isnan(sp) || sp == 0.0
            sp = 1e-9
        end
        DataFrame(logPrice = log(m/v), Price = m/v, Volume = v, diffPrice = md/v, n= length(df.Volume), sp = sp )
    end
end

function trainMatrix(pSize::Int, dfI::DataFrame, nSize = 0, vSize = 0, lpSize = 0 )

    ns = size(dfI)[1]-pSize+1
    cM = Array{Float32,2}(undef,nSize+vSize+lpSize+pSize, ns-1)
    for i = 2:ns
        offset = 0
        c=Vector{Float32}(undef,pSize+nSize+vSize+lpSize)
        c[1:nSize] .=  log.(dfI.n[i+pSize-nSize-1:i-2+pSize]) .- log.(dfI.n[i+pSize-nSize-2:i-3+pSize])
        offset = nSize
        c[offset+1:offset+vSize] .= log.(dfI.Volume[i+pSize-vSize-1:i-2+pSize]) .-log.(dfI.Volume[i+pSize-vSize-2:i-3+pSize])# ./log(dfA.Volume[i-1]) .-1.0
        offset += vSize
        # c[offset+1:offset+lpSize] .= (dfA.logPrice[i+pSize-lpSize-1:i-2+pSize]./ dfA.logPrice[i-1] .- 1.0)*100.0
        # offset += lpSize
        #
        c[offset+1:offset+pSize] .= (dfI.logPrice[i:i-1+pSize].- dfI.logPrice[i-1:i-2+pSize])*100.0#or ./norm and -1.0
        #c[offset+1:offset+pSize] .= (dfA.logPrice[i:i-1+pSize]./ dfA.logPrice[i-1:i-2+pSize] .- 1.0)*100.0#or ./norm and -1.0
        #(dfA.logPrice[i:i-1+pSize] ./ dfA.logPrice[i-1].-1.0)*100.0#or ./norm and -1.0
        offset+=pSize
        cM[:,i-1] .= atan.(c)
    end
    return cM
end

function trainMatrixF(pSize::Int, dfI::DataFrame, nSize = 0, vSize = 0, lpSize = 0 )

    ns = size(dfI)[1]-2*nSize+1
    cM = Array{Float32,2}(undef,pSize, ns-1)

    function atanPow(x::Float64, p::Int )
        if p>0
            return erf(atanPow(x,p-1))
        else
            return x
        end
    end
    for i = 2:ns
        dv = ((dfI.logPrice[i:i-1+nSize].- dfI.logPrice[i-1:i-2+nSize])*200.0)
        c=Vector{Float32}(undef,pSize)

        c[1] = mean(dv[1:end-1])
        for j in 2:pSize-1
            c[j] = moment(dv[1:end-1],j)  #mean( x->x^j , dv[1:end-1] )#moment(dv[1:end-1],j) #atan(mean( x->x^j , dv[1:end-1] ))
        end
        c[2:pSize-1] ./= c[1:pSize-2]
#        c[end-1] = dv[end-1]
        c[end] = dv[end]
        #c[end] = (dfI.logPrice[i-1+2*nSize]- dfI.logPrice[i-2+nSize])*100.0
        #c .*= 100.0

        cM[:,i-1] .= atan.(c)
    end
    return cM
end

function trainMatrixFD(pSize::Int, dfI::DataFrame )

    ns = size(dfI)[1]-pSize+1
    cM = Array{Float32,2}(undef,pSize+pSize-1, ns-1)

    for i = 2:ns
        dv1 = (dfI.logPrice1[i:i-1+pSize].- dfI.logPrice1[i-1:i-2+pSize])*160.0
        dv2 = (dfI.logPrice2[i:i-2+pSize].- dfI.logPrice2[i-1:i-3+pSize])*160.0
        c=Vector{Float32}(undef, 2*pSize-1)

        for j in 1:pSize-1
            c[j] = atan(mean( x->x^j , dv1[1:end-1] ))
            c[j+pSize-1] = atan(mean( x->x^j , dv2[1:end] ))
        end

        c[end] = dv1[end]
        #c .*= 100.0

        cM[:,i-1] .= atan.(c)
    end
    return cM
end

function joinDF(df1::DataFrame, df2::DataFrame)
    names1 = names(df1)
    names2 = names(df2)
    renames1 = Pair{String,String}[]
    renames2 = Pair{String,String}[]
    gnames = ["year","day","hour"] #,"divMinute"]
    if "divMinute" in renames1
        push!(gnames,"divMinute")
    end

    for n in names1
        if  !(n in gnames)
            push!(renames1,Pair(n,string(n,"1")))
        end
    end
    for n in names2
        if !(n in gnames)
            push!(renames2,Pair(n,string(n,"2")))
        end
    end
    leftjoin( rename(df1,renames1), rename(df2,renames2), on = gnames )
end
