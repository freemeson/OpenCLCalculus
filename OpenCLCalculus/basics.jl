using LinearAlgebra
using SpecialFunctions
using StatsBase

mutable struct literal
        definitions::String
        expression::String
        input::Vector{Tuple{String, String, Int, String, Union{Nothing,String}}}
end

mutable struct referenceCount
        id::Int
        eval::String
        attribute::String
        count::Int
        headerType::Union{Nothing,String}
        extraAttributes::Union{Nothing,Vector{String}}
        function referenceCount( r::Tuple{Int, String, String}  )
                new(r[1], r[2], r[3], 1, nothing, nothing)
        end
        function referenceCount( r::Tuple{Int, String, String, String}  )
                new(r[1], r[2], r[3], 1, r[4], nothing)
        end
end

function Base.:setindex!(a::Dict{String,referenceCount}, b::Tuple{Int, String, String}, key::String)
        return a[key] = referenceCount(b)
end

function Base.:setindex!(a::Dict{String,referenceCount}, b::Tuple{Int, String, String, String}, key::String)
        return a[key] = referenceCount(b)
end

mutable struct hiddenCounter
        counter::Int
        #expressions::Dict{String, String}
        exprList::Dict{String, referenceCount}
        prodList::Dict{String, Vector{String}}
        function hiddenCounter()
                new(0, Dict{String, Tuple{Int, String, String}}(), Dict{String, Vector{String}}())
        end
end

function dataAttributes(hc::hiddenCounter)
        #return map(x->(x[1],x[2][2], -x[2][1], x[2][3]),  collect(filter(x->x[2][1] < 0,  hc.exprList)) )
        return map(x-> if x[2].extraAttributes == nothing
                                return (x[1],x[2].eval, -x[2].id, x[2].attribute, x[2].headerType)
                        elseif x[2].extraAttributes[1] == "extraInput"
                                return (x[1],x[2].eval, -1, x[2].attribute, x[2].headerType)
                        else return (x[1],x[2].eval, -2, x[2].attribute, x[2].headerType)
                        end
                        ,
                        collect(filter(x->x[2].id < 0 || x[2].attribute == "output" || x[2].attribute == "function" || x[2].attribute == "function declaration", hc.exprList))
                )
end

abstract type ScalarDataIndex end
abstract type VectorDataIndex end

mutable struct X <: VectorDataIndex
        counter::Int
        ref::Vector{ScalarDataIndex}
        name::String
        defaultAttribute::Union{Nothing, String}
        function X(name = "inputVector", defaultAttribute::Union{Nothing, String}="unspecified")
                new(0,ScalarDataIndex[], name, defaultAttribute)
        end
end

mutable struct dataIndex <: ScalarDataIndex
        id::Int
        name::String
        hash::Union{Nothing,String}
        attribute::Union{Nothing, String}
        function dataIndex(id::Int,name::String, attribute::Union{Nothing, String} = "unspecified")
                new(id,name,nothing, attribute)
        end
end

mutable struct constIndex <: ScalarDataIndex
        id::Int
        name::String
        value
        hash::Union{Nothing,String}
        function constIndex(id::Int,name::String,value)
                new(id,name,value,nothing)
        end
end

function newVar(x::X, attribute::Union{Nothing, String} = nothing)
        x.counter += 1
        if attribute == nothing
                attribute = x.defaultAttribute
        end
        push!(x.ref,dataIndex(x.counter,x.name, attribute))
        return last(x.ref)
end

function newVars(x::X, n::Int, attribute::Union{Nothing, String}=nothing)
        if attribute == nothing
                attribute = x.defaultAttribute
        end
        for i in 1:n
                push!(x.ref,dataIndex(x.counter+i,x.name, attribute))
        end
        x.counter += n

        return x.ref[end-n+1:end]
end

function Base.:getindex(x::X,r::Int)
        if r > x.counter
                for i in x.counter+1:r
                        push!(x.ref,dataIndex(i,x.name,x.defaultAttribute))
                end
                x.counter = r
        end
        return x.ref[r]
end

function Base.:getindex(x::X,r::UnitRange{Int})
        if length(r) > 0 && r.stop > x.counter
                for i in x.counter+1:r.stop
                        push!(x.ref,dataIndex(i,x.name,x.defaultAttribute))
                end
                x.counter = r.stop
        end
        return x.ref[r]
end

mutable struct dotProdExpr <: ScalarDataIndex
        id::Union{Int,Nothing}
        name::String
        pos::Vector{ScalarDataIndex}
        par::Vector{ScalarDataIndex}
        hash::Union{Nothing,String}
        function dotProdExpr(pos::Vector{T},par::Vector{S}, id::Union{Int,Nothing} = nothing) where T<:ScalarDataIndex where S<:ScalarDataIndex
                new(id, "dotProdExpr", pos, par, nothing)
        end
end

mutable struct sumExpr <: ScalarDataIndex
        id::Union{Int, Nothing}
        name::String
        px::Vector{ScalarDataIndex}
        hash::Union{Nothing,String}
        function sumExpr(px::Vector{T}, id::Union{Int,Nothing} = nothing) where T<:ScalarDataIndex
                sumPx = filter(x->typeof(x) == sumExpr,px)
                tempPx = filter(x->typeof(x) != sumExpr,px)
                #sum of dorProdExprs could be a dotProdExr
                if length(sumPx) > 0
                        tempPx = vcat(tempPx, [s.px for s in sumPx]...)
                end

                constPx = filter(x->typeof(x) == constIndex,tempPx)
                exprPx = filter(x->typeof(x) != constIndex,tempPx)
                if length(constPx)!= 0
                        sumCPX = sum(c->c.value, constPx)
                        if sumCPX == 0.0 && length(exprPx)!=0
                                return new(id,"sumExpr",exprPx, nothing)
                        else
                                return new(id,"sumExpr",vcat(constIndex(0,"sumExprConst",sumCPX), exprPx), nothing)
                        end
                else
                        return new(id,"sumExpr", exprPx, nothing)
                end
        end
end

mutable struct prodExpr <: ScalarDataIndex
        id::Union{Int, Nothing}
        name::String
        px::Vector{ScalarDataIndex}
        hash::Union{Nothing,String}
        keepAsBlock::Bool
        function prodExpr(px::Vector{T}, id::Union{Int,Nothing} = nothing, keepAsBlock::Bool = false) where T<:ScalarDataIndex

                nominator = ScalarDataIndex[]
                denominator = ScalarDataIndex[]
                collectDivisors(px, nominator, denominator)
                #simplDivExpr = vcat(simplifyExpr(divExpr(prodExpr([d.nominator for d in divs]), prodExpr([d.denominator for d in divs])), nonDivExpr))


                constPx = filter(x->typeof(x) == constIndex,nominator)
                exprPx = filter(x->typeof(x) != constIndex,nominator)
                if length(denominator) != 0
                        simpleProd = prodExpr(denominator)
                        denomConst = filter(x->typeof(x) == constIndex, simpleProd.px) # should be at most a sinlge element
                        noConstProd = filter(x->typeof(x) != constIndex, simpleProd.px)
                        if length(noConstProd) > 1
                                simpleProd.px = noConstProd
                                onlyDiv = divExpr(constIndex(0, "prodExprConst", 1.0), simpleProd)
                                push!(exprPx, onlyDiv)
                        elseif length(noConstProd) == 1
                                onlyDiv = divExpr(constIndex(0, "prodExprConst", 1.0), noConstProd[1])
                                push!(exprPx, onlyDiv)
                        end
                        if length(denomConst) == 1 #it really shouldn't be >1!!!
                                push!(constPx, constIndex(0, "prodExprConst", 1.0/denomConst[1].value) )
                        end
                end

                if length(constPx)!= 0
                        prodCPX = prod(c->c.value, constPx)
                        if prodCPX!=0
                                if prodCPX == 1.0 && length(exprPx)!=0
                                        return new(id,"prodExpr", exprPx, nothing, keepAsBlock )
                                else
                                        return new(id,"prodExpr",vcat(constIndex(0,"prodExprConst",prodCPX), exprPx), nothing, keepAsBlock )
                                end
                        else
                                return new(id,"prodExpr", [constIndex(0,"prodExprConst",0.0)], nothing, keepAsBlock )
                        end
                else
                        return new(id,"prodExpr", exprPx, nothing, keepAsBlock )
                end


                #new(id,"prodExpr", px)
        end
end

mutable struct divExpr <: ScalarDataIndex
        id::Union{Int, Nothing}
        name::String
        nominator::ScalarDataIndex
        denominator::ScalarDataIndex
        hash::Union{Nothing,String}
        function divExpr(nominator::ScalarDataIndex,denominator::ScalarDataIndex, id::Union{Int,Nothing} = nothing)


                nom = ScalarDataIndex[]
                den = ScalarDataIndex[]

                if typeof(nominator) in [divExpr, prodExpr]
                        collectDivisors(nominator, nom, den )
                else
                        push!(nom, nominator)
                end

                if typeof(denominator) in [divExpr, prodExpr]
                        collectDivisors(denominator, den, nom)
                else
                        push!(den, denominator)
                end

                noms = simplifyExpr(prodExpr(nom))
                dens = simplifyExpr(prodExpr(den))


                # if typeof(nom) == divExpr && typeof(den) != divExpr
                #         new(id,"divExpr", nom.nominator, prodExpr([den,nom.denominator]), nothing)
                # end
                #
                # if typeof(nom) != divExpr && typeof(den) == divExpr
                #         new(id, "divExpr", prodExpr( [nom, den.denominator] ), den.nominator , nothing)
                # end
                #
                # if typeof(nom) == divExpr && typeof(den) == divExpr
                #         new(id, "divExpr", prodExpr([nom.nominator, den.denominator]), prodExpr([nom.denominator, den.nominator]) , nothing)
                # end
                new(id,"divExpr", noms, dens, nothing)
        end
end



function simplifyExpr(d::dataIndex)
        return d
end

function simplifyExpr(c::constIndex)
        return c
end

function simplifyExpr(e::sumExpr)
        px = simplifyExpr.(e.px)

        #sum of const*products
        #sum of divisions
        prodPx =    filter(x->typeof(x) == prodExpr, px)
        dataPx = filter(x->typeof(x) == dataIndex, px)
        nonProdPx = filter(x->typeof(x) != prodExpr && typeof(x)!=dataIndex, px)

        if length(prodPx)+length(dataPx) > 1
                reorder.(prodPx) #convert to ScalarDataIndex[]
                prodPx = convert(Vector{ScalarDataIndex}, prodPx)
                # println(typeof(prodPx))
                # println(prodExpr.([[d] for d in dataPx]))
                append!(prodPx, prodExpr.([[d] for d in dataPx]))
                elements = Tuple{Real, Vector{ScalarDataIndex}, String}[] #Dict{String,Vector{Tuple{Real, ScalarDataIndex}}}()
                for p in prodPx
                        if typeof(p.px[1]) == constIndex
                                push!(elements, (p.px[1].value, p.px[2:end], string(p.px[2:end])) )
                        else
                                push!(elements, (1.0, p.px, string(p.px)))
                        end
                end
                uniqueHashes = unique([e[3] for e in elements])
                newPx = ScalarDataIndex[]
                if length(uniqueHashes)!=length(prodPx)
                        for uh in uniqueHashes
                                similarPXs = findall(e->e[3]==uh, elements)
                                overallConst = sum([ elements[s][1] for s in similarPXs])
                                #println(similarPXs)
                                if overallConst!=0.0
                                        #println(overallConst)
                                        if overallConst != 1.0
                                                push!(newPx, prodExpr( [constIndex(0, "simplSumConst", overallConst), elements[similarPXs[1]][2]...] ) )
                                        else
                                                if length(elements[similarPXs[1]][2])>1
                                                        push!(newPx, prodExpr(elements[similarPXs[1]][2]) )
                                                else
                                                        push!(newPx, elements[similarPXs[1]][2][1])
                                                end
                                        end
                                else
                                        push!(newPx, constIndex(0, "simplSumCancelConst", 0.0))
                                end
                        end

                        prodPx = simplifyExpr.(newPx)
                else
                        prodPx = simplifyExpr.(prodPx)
                end
        else
                prodPx = vcat(prodPx, dataPx)
        end

        px = vcat(prodPx, nonProdPx)
        divPx = filter(x->typeof(x)==divExpr, px)
        nonDivPx = filter(x->typeof(x)!=divExpr, px)

        if length(divPx) != 0
                reorder.(divPx)
                #println(length(divPx))
                elements = [ (d.nominator, d.denominator, d.denominator.hash ) for d in divPx  ]
                uniqueHashes = unique([e[3] for e in elements])
                newPx = ScalarDataIndex[]
                if length(uniqueHashes)!=length(divPx)
                        for uh in uniqueHashes
                                similarPXIndices = findall(e->e[3]==uh, elements)
                                similarPX = [e[1] for e in elements[ similarPXIndices ] ]
                                if length(similarPX) > 1
                                        push!(newPx, divExpr(sumExpr(similarPX), elements[similarPXIndices[1]][2]))
                                else
                                        push!(newPx, divExpr(similarPX[1], elements[similarPXIndices[1]][2] ))
                                end
                        end
                        divPx = simplifyExpr.(newPx)
                else

                        divPx = simplifyExpr.(divPx)
                end
        end

        px = vcat(divPx, nonDivPx)
        se = sumExpr(px,e.id)

        if length(se.px) == 1
                return se.px[1]
        else
                return se
        end
end

function simplifyExpr(e::prodExpr)
        px = simplifyExpr.(e.px)
#        println(px
        # if typeof(px) == Array{Any,1}
        #         println("dddd")
        #         println(length(px))
        #         println(length(e.px))
        #         println(e.px)
        #         for p in px
        #                 println("  ",typeof(p))
        #         end
        # end
        pe = prodExpr(px,e.id)

        if length(pe.px) == 1
                return pe.px[1]
        end

        divPx = filter(x->typeof(x) == divExpr,pe.px)
        exprPx = filter(x->typeof(x) != divExpr,pe.px)

        if length(divPx)!=0
                nom = [d.nominator for d in divPx]
                den = [d.denominator for d in divPx]
                return simplifyExpr(divExpr(prodExpr(vcat(nom,exprPx) ),
                                prodExpr(den),e.id ))
        else
                return pe
        end
end

function simplifyExpr(e::dotProdExpr)
        p = simplifyExpr.(e.pos)
        s = simplifyExpr.(e.par)
        if length(p) == 1
                return simplifyExpr(prodExpr([p[1], s[1]]))
        end
        if unique(typeof.(p)) == [constIndex] && unique(typeof.(s)) == [constIndex]
                pv = [c.value for c in p]
                sv = [c.value for c in s]
                return constExpr(e.id,string("simpl",e.id),dot(pv,sv) )
        else
                return dotProdExpr(p,s,e.id)
        end
end

function simplifyExpr(d::divExpr)
        # nom = ScalarDataIndex[]
        # den = ScalarDataIndex[]
        #
        # collectDivisors(d, nom, den)

        nom = simplifyExpr(d.nominator)
        den = simplifyExpr(d.denominator)

        if typeof(nom) == constIndex && typeof(den) == constIndex
                return constIndex(0,"simplDivConst",nom.value/den.value)
        end
        if typeof(nom) == constIndex && nom.value == 0.0
                return constIndex(0,"simplDivConst", 0.0) #assuming the denominator is never zero! :D
        end

        if typeof(den) == constIndex
                return prodExpr([nom, constIndex(0,"simplDivConst", 1.0/den.value)]   )
        end

        #if (typeof(nom) == dataIndex || typeof(nom) == prodExpr) && (typeof(den) == dataIndex || typeof(den) == prodExpr)
        if typeof(nom) == dataIndex && typeof(den) == dataIndex && nom.name == den.name && nom.id == den.id
                return constIndex(0,"simplDivConstData",1.0)
        end

        if typeof(nom) == typeof(den) && string(nom) == string(den) #constIndex was ecluded, but anyway, it works
                return constIndex(0,"simplDivConstExpr",1.0)
        end

        if typeof(nom) == constIndex && typeof(den) == prodExpr
                fIndex = findfirst(x->typeof(x) == constIndex, den.px)
                if fIndex != nothing
                        n = constIndex(0, "simplDivConstProd", nom.value/den.px[fIndex].value)
                        deleteat!(den.px, fIndex)
                        return divExpr(n, den)
                end
        end

        if typeof(den) == constIndex && typeof(nom) == prodExpr
                fIndex = findfirst(x->typeof(x) == constIndex, nom.px)
                if fIndex != nothing
                        nom.px[fIndex] = constIndex(0, "simplDivConstProd", nom.px[fIndex].value/den.value)
                        return nom
                else
                        return divExpr(nom,den,d.id)
                end
        end

        if typeof(den) !=  prodExpr && typeof(nom) == prodExpr
                fIndex = findfirst( x->x==den, nom.px)
                if fIndex != nothing
                        deleteat!(nom.px, fIndex)
                        if length(nom.px) == 1
                                return nom.px[1]
                        else
                                return nom
                        end
                else
                        return divExpr(nom,den,d.id)
                end
        end

        if typeof(nom) != prodExpr && typeof(den) == prodExpr
                fIndex = findfirst( x->x==nom, den.px)
                if fIndex != nothing
                        deleteat!(den.px, fIndex)
                        if length(den.px) == 1
                                return divExpr(  constIndex(0, "simplDivConst", 1.0) , den.px[1])
                        else
                                return divExpr(  constIndex(0, "simplDivConst", 1.0) , den )
                        end
                else
                        return divExpr(nom, den, d.id)
                end
        end

        if typeof(nom) == prodExpr && typeof(den) == prodExpr
                reorder(nom)
                reorder(den)
                if nom.hash == den.hash
                        return constIndex(0, "simplDivConst", 1.0)
                end
                nomConst = findfirst(x->typeof(x)==constIndex, nom.px)
                denConst = findfirst(x->typeof(x)==constIndex, den.px)
                if nomConst!=nothing && denConst!=nothing
                        nom.px[nomConst].value = nom.px[nomConst].value/den.px[denConst].value
                        deleteat!(den.px, denConst)
                end

                nomDeletes = Int[]
                denDeletes = Int[]

                denIndices = zeros(Int,length(den.px))
                for i in 1:length(nom.px)
                        found = findfirst(x->x.hash == nom.px[i].hash, den.px )
                        lastfind = found
                        while lastfind != nothing && denIndices[lastfind] != 0
                                found = findnext(x->x.hash == nom.px[i].hash, den.px, lastfind+1)
                                lastfind = found
                        end
#                        println(nom.px[i].hash)


                        if found != nothing
                                denIndices[found] = 1
                                push!(nomDeletes, i)
                                push!(denDeletes, found)
                        end
                end
                if length(nomDeletes) != 0
                        if length(denDeletes) != length(nomDeletes)
                                println("error in simplification!")
                        else
                                deleteat!(nom.px, nomDeletes)
                                #println(denDeletes)
                                deleteat!(den.px, sort(denDeletes))
                                if length(den.px) == 0 && length(nom.px) == 0
                                        return constIndex(0, "simplDivConst",1.0)
                                elseif length(den.px) == 0
                                        return simplifyExpr(prodExpr(nom.px))
                                elseif length(nom.px) == 0
                                        return simplifyExpr(divExpr(constIndex(0, "simplDivConst",1.0), prodExpr(den.px))) # as den.px could be a constIndex
                                else
                                        return simplifyExpr(divExpr(simplifyExpr(prodExpr(nom.px)), simplifyExpr(prodExpr(den.px))   ))
                                end

                        end
                end

        end

        return divExpr(nom,den,d.id)
end


function collectDivisors(px::Vector{T}, nom::Vector{ScalarDataIndex},den::Vector{ScalarDataIndex}) where T <: ScalarDataIndex
        for p in px
                if typeof(p) in [prodExpr, divExpr]
                        collectDivisors(p, nom, den)
                else
                        push!(nom, p)
                end
        end
end

function collectDivisors(pe::prodExpr, nom::Vector{ScalarDataIndex}, den::Vector{ScalarDataIndex})
        if pe.keepAsBlock
                push!(nom,pe)
        else
                collectDivisors(pe.px, nom, den)
        end
end

function collectDivisors(de::divExpr, nom::Vector{ScalarDataIndex}, den::Vector{ScalarDataIndex})
        if typeof(de.nominator) in [prodExpr, divExpr]
                collectDivisors(de.nominator, nom, den)
        else
                push!(nom, de.nominator)
        end

        if typeof(de.denominator) in [prodExpr, divExpr]
                collectDivisors(de.denominator, den, nom)
        else
                push!(den, de.denominator)
        end
end

function createConstIndices(x::X)
        xConst = Vector{ScalarDataIndex}(undef,x.counter)
        for i in 1:length(xConst)
                xConst[i] = constIndex(i,string("const",x.name),0.0)
        end
        return xConst
end

function createAbstractIndices(x::X)
        xConst = Vector{ScalarDataIndex}(undef,x.counter)
end

function sameTypeIsLess(a::constIndex,b::constIndex)
        if a.value!=b.value
                return a.value < b.value
        elseif a.name != b.name
                return a.name < b.name
        else
                return a.id < b.id
        end
end

function sameTypeIsLess(a::dataIndex,b::dataIndex)
        if a.name != b.name
                return a.name < b.name
        else
                return a.id < b.id
        end
end

function sameTypeIsLess(a::sumExpr, b::sumExpr)
        if a.id != nothing && b.id != nothing
                return a.id < b.id
        elseif length(a.px) != length(b.px)
                return length(a.px) < length(b.px)
        else
                return string(a.px) < string(b.px)
        end
end

function sameTypeIsLess(a::prodExpr, b::prodExpr)
        if a.id != nothing && b.id != nothing
                return a.id < b.id
        elseif length(a.px) != length(b.px)
                return length(a.px) < length(b.px)
        else
                return string(a.px) < string(b.px)
        end
end

function sameTypeIsLess(a::divExpr, b::divExpr)
        if a.id != nothing && b.id != nothing
                return a.id < b.id
        end
        a.denominator.hash = string(a.denominator)
        b.denominator.hash = string(b.denominator)

        if a.denominator.hash < b.denominator.hash
                return a.denominator.hash < b.denominator.hash
        else
                a.nominator.hash = string(a.nominator)
                b.nominator.hash = string(b.nominator)
                # a.hash = string(a)
                # b.hash = string(b)

                return a.nominator.hash < b.nominator.hash
        end
end

function sameTypeIsLess(a::dotProdExpr, b::dotProdExpr)
        if a.id != nothing && b.id != nothing
                return a.id < b.id
        elseif length(a.pos) != length(b.pos)
                return length(a.pos) < length(b.pos)
        else
                return string(a) < string(b)
        end
end

function Base.:isless(a::ScalarDataIndex, b::ScalarDataIndex)
        if typeof(a) == typeof(b)
                return sameTypeIsLess(a,b)
        else
                return exprOrder[typeof(a)] < exprOrder[typeof(b)]
        end
end

#the purpose of reorder is to generate fingerprint hash, that identifies expressions
function reorder(pe::prodExpr)
        pe.px = sort(pe.px)
        reorder.(pe.px)
        pe.hash = string(pe)
        return nothing
end

function reorder(de::divExpr)
        reorder(de.nominator)
        reorder(de.denominator)
        de.hash = string(de)
        return nothing
end

function reorder(se::sumExpr)
        se.px = sort(se.px)
        reorder.(se.px)
        se.hash = string(se)
        return nothing
end

function reorder(dpe::dotProdExpr)
        #sorting pos and par is allowed, but complicated
        reorder.(dpe.pos)
        reorder.(dpe.par)
        dpe.hash = string(dpe)
        return nothing
end

function reorder(d::dataIndex)
        d.hash = string(d)
        return nothing
end

function reorder(c::constIndex)
        c.hash = string(c)
        return nothing
end

function Base.:string(c::constIndex)
        string(c.value)
end

function Base.:string(d::dataIndex)
        string(d.name,"[", d.id,"]")
end

function Base.:string(dv::Vector{T}) where T <: ScalarDataIndex

        if length(dv) == 1
                return string(dv[1])
        else

                if unique(typeof.(dv)) == [dataIndex] &&
                        issorted([d.id for d in dv]) &&
                        #(dv[end].id - dv[1].id) == length(dv)-1 &&
                        [d.id for d in dv] == collect(dv[1].id:dv[end].id) &&
                        unique([d.name for d in dv])==[dv[1].name]

#                        println(dv[end].id, "-", dv[1].id, " l=", length(dv))
                        return string(dv[1].name, "[",dv[1].id,":",dv[end].id, "]")
                else

                        return string("[",[string(string(d),",") for d in dv]...,"]")
                end
        end
end

function Base.:string(se::sumExpr)
        if se.id==nothing
                return string("sum(",string(se.px), ")")
        else
                return string("sumExpr",se.id)
        end
end

function Base.:string(pe::prodExpr)
        if pe.id==nothing
                return string("prod(",string(pe.px), ")")
        else
                return string("prodExpr",pe.id)
        end
end

function Base.:string(de::divExpr)
        if de.id==nothing
                return string(string(de.nominator),"/(",string(de.denominator), ")")
        else
                return string("divExpr",de.id)
        end
end

function Base.:string(dpe::dotProdExpr)
        if dpe.id==nothing
                return string("dot(",string(dpe.pos),",",string(dpe.par), ")")
        else
                return string("dotProdExpr",dpe.id)
        end
end

function derive(expr::dataIndex, byX::X, byP::X) #could be a single X vector, but 2 is more useful
        #find exp.name in byX or byP, return nonzero for only that element
        xConst = createConstIndices(byX)
        pConst = createConstIndices(byP)

        if expr.name == byX.name
                xConst[expr.id].value = 1.0
        end
        if expr.name == byP.name
                pConst[expr.id].value = 1.0
        end
        return xConst,pConst
end

function derive(expr::constIndex, byX::X, byP::X) #should have a hiddenCounter too, so the names won't collide
        xConst = createConstIndices(byX)
        pConst = createConstIndices(byP)

        return xConst,pConst
end

function derive(expr::sumExpr, byX::X, byP::X)
        xExpr = createAbstractIndices(byX)
        pExpr = createAbstractIndices(byP)

        exprVectors = map(x->derive(x,byX,byP), expr.px)

        for i in 1:byX.counter
                xExpr[i] = sumExpr([expr[1][i] for expr in exprVectors ])
        end

        for i in 1:byP.counter
                pExpr[i] = sumExpr([expr[2][i] for expr in exprVectors ])
        end

        return simplifyExpr.(xExpr),simplifyExpr.(pExpr)
end

function derive(expr::prodExpr, byX::X, byP::X)
        xExpr = createAbstractIndices(byX)
        pExpr = createAbstractIndices(byP)

        exprVectors = map(x->derive(x,byX,byP), expr.px)

        for i in 1:byX.counter
                #sum(prod(replace(myCopyPX[i],i) )
                subExprVecX = Vector{ScalarDataIndex}(undef,length(expr.px))

                for j in 1:length(expr.px)
                        myCopyX = deepcopy(expr.px)
                        myCopyX[j] = exprVectors[j][1][i]

                        subExprVecX[j] = prodExpr(myCopyX)
                end
                xExpr[i] = sumExpr(subExprVecX)

        end

        for i in 1:byP.counter
                #sum(prod(replace(myCopyPX[i],i) )

                subExprVecP = Vector{ScalarDataIndex}(undef,length(expr.px))
                for j in 1:length(expr.px)

                        myCopyP = deepcopy(expr.px)
                        myCopyP[j] = exprVectors[j][2][i]

                        subExprVecP[j] = prodExpr(myCopyP)
                end

                pExpr[i] = sumExpr(subExprVecP)
        end

        return simplifyExpr.(xExpr),simplifyExpr.(pExpr)
end

function derive(expr::dotProdExpr, byX::X, byP::X)
        # typeS = unique(typeof.(expr.s))
        # typeP = unique(typeof.(expr.p))
#        if  (typeS == [dataIndex] || typeS == [constIndex]) && (typeR ==[dataIndex] || typeR == [constIndex])

        if unique([x.name for x in expr.pos]) == [byP.name] && unique([p.name for p in expr.par]) == [byX.name]
                #swap in case there is was a mix up
                pos = expr.par
                par = expr.pos
        else
                pos = expr.pos
                par = expr.par
        end

        #for most cases, this is the fast derivative of dotProdExpr
        if unique([x.name for x in pos]) == [byX.name] && unique(typeof.(pos)) == [dataIndex]
                if unique([p.name for p in par]) == [byP.name] && unique(typeof.(par)) == [dataIndex]
                        xExpr = createConstIndices(byX)
                        pExpr = createConstIndices(byP)

                        # xExpr = createAbstractIndices(byX)
                        # pExpr = createAbstractIndices(byP)
                        # for i in 1:length(xExpr)
                        #         xExpr[i] = constIndex(i,string("const",byX.name),0.0)
                        # end
                        for i in 1:length(par)
                                xExpr[pos[i].id] = par[i]
                        end

                        # for i in 1:length(pExpr)
                        #         pExpr[i] = constIndex(i,string("const",byP.name),0.0)
                        # end
                        for i in 1:length(par)
                                pExpr[par[i].id] = pos[i]
                        end
                        return xExpr,pExpr
                end
        end

        #this is the general expression to derive a dotProdExpr
        gExpr = Vector{ScalarDataIndex}(undef, length(pos))
        for i in 1:length(gExpr)
                gExpr[i] = prodExpr([pos[i],par[i]])
        end

        return derive(simplifyExpr(sumExpr(deepcopy(gExpr))) , byX, byP)

end

function derive(expr::divExpr, byX, byP)

        nom = deepcopy(simplifyExpr(expr.nominator))
        den = deepcopy(simplifyExpr(expr.denominator))

        #special caseses (should write one for const*div(a,b)
        if typeof(nom) == dataIndex && typeof(den) == dataIndex
                 if nom.name == byX.name && den.name == byX.name
                        xExpr = createConstIndices(byX)
                        pExpr = createConstIndices(byP)

                        if den.id != nom.id
                                xExpr[nom.id] = divExpr(constIndex(0,"deriveDivConst",1.0), den)
                                xExpr[den.id] = divExpr(nom,prodExpr([constIndex(1,"deriveDivConst",-1.0),den,den]) )
                                return xExpr, pExpr
                        else
                                return xExpr, pExpr
                        end
                end
                if nom.name == byP.name && den.name == byP.name
                       xExpr = createConstIndices(byX)
                       pExpr = createConstIndices(byP)

                       if den.id != nom.id
                               pExpr[nom.id] = divExpr(constIndex(0,"deriveDivConst",1.0), den)
                               pExpr[den.id] = divExpr(nom,prodExpr([constIndex(1,"deriveDivConst",-1.0),den,den]) )
                               return xExpr, pExpr
                       else
                               return xExpr, pExpr
                       end
               end

                if nom.name == byX.name && den.name == byP.name
                        xExpr = createConstIndices(byX)
                        pExpr = createConstIndices(byP)
                        xExpr[nom.id] = divExpr(constIndex(0,"deriveDivConst",1.0), den)
                        pExpr[den.id] = divExpr(nom,prodExpr([constIndex(1,"deriveDivConst",-1.0),den,den]) )
                        return xExpr, pExpr
                end
                if nom.name == byP.name && den.name == byX.name
                        xExpr = createConstIndices(byX)
                        pExpr = createConstIndices(byP)
                        pExpr[nom.id] = divExpr(constIndex(0,"deriveDivConst",1.0), den)
                        xExpr[den.id] = divExpr(nom,prodExpr([constIndex(1,"deriveDivConst",-1.0),den,den]) )
                        return xExpr, pExpr
                end
        end
        dNom = derive(nom, byX, byP)
        dDen = derive(den, byX, byP)

        xExpr = createConstIndices(byX)
        pExpr = createConstIndices(byP)


        for i in 1:byX.counter
                if !(typeof(dNom[1][i]) == constIndex && dNom[1][i].value==0.0)
                        xExpr[i] = divExpr(dNom[1][i],den)
                end
        end

        for i in 1:byP.counter
                if !(typeof(dNom[2][i]) == constIndex && dNom[2][i].value==0.0)
                        pExpr[i] = divExpr(dNom[2][i],den)
                end
        end

        for i in 1:byX.counter
                if !(typeof(dDen[1][i]) == constIndex && dDen[1][i].value==0.0)
                        tempExpr = divExpr(prodExpr([constIndex(0,"derDivConst",-1.0  ),nom,dDen[1][i]]), prodExpr([ den,den ] ) )
                        xExpr[i] = sumExpr([ xExpr[i], tempExpr ] )
                end
        end

        for i in 1:byP.counter
                if !(typeof(dDen[2][i]) == constIndex && dDen[2][i].value==0.0)
                        tempExpr = divExpr(prodExpr([constIndex(0,"derDivConst",-1.0  ),nom,dDen[2][i]]), prodExpr([ den,den ] ) )
                        pExpr[i] = sumExpr([ pExpr[i], tempExpr ] )
                end
        end

        return xExpr,pExpr
end

exprOrder = Dict([ constIndex => 1, dataIndex => 2, sumExpr => 3, prodExpr => 4, divExpr => 5, dotProdExpr => 6])


function juliaLiteral(c::constIndex, hc::hiddenCounter = hiddenCounter())
        c.hash = string(c)
        if !haskey(hc.exprList, c.hash)
                #const are not so important, maybe they don't need to be tracked
                hc.exprList[c.hash] = (0,c.name, string(c.value) )
        end
        literal("", string(c.value), dataAttributes(hc)  )
end

function juliaLiteral(d::dataIndex, hc::hiddenCounter = hiddenCounter())
        # if hc.counter == 0
        #         firstPass(d, hc)
        # end
        d.hash = string(d)
        if !haskey(hc.exprList, d.hash)
                hc.exprList[d.hash] = (-d.id,d.name, string(d.attribute), "Float") #I want to track if a dataIndex was used or is obsolete
        end
        literal("",string(d),dataAttributes(hc))
end

function juliaLiteral(dv::Vector{T}, hc::hiddenCounter = hiddenCounter()) where T<:ScalarDataIndex

        if length(dv) == 1
                return juliaLiteral(dv[1], hc)
        else
                if unique(typeof.(dv)) == [dataIndex]
                        for d in dv
                                d.hash = string(d)
                                if !haskey(hc.exprList, d.hash)
                                        hc.exprList[d.hash] = (-d.id,d.name, string(d.attribute), "Float") #I want to track if a dataIndex was used and isn't obsolete
                                end
                        end
                        if issorted([d.id for d in dv]) && length(unique([d.name for d in dv]))==1 && [d.id for d in dv] == collect(dv[1].id:dv[end].id)

                                return literal("",string(dv[1].name, "[",dv[1].id,":",dv[end].id, "]"), dataAttributes(hc) )
                        else
                                return literal("",string("[",[string(string(d),",") for d in dv]...,"]"), dataAttributes(hc) )
                        end
                else
                        lv = [juliaLiteral(d, hc) for d in dv]
                        #hc tracks if a variable was defined before, and removes it from l.definitions, sets correct name in l.expressions
                        return literal(string([l.definitions for l in lv]...), string("[",[string(l.expression,",") for l in lv]...,"]") , dataAttributes(hc) )
                end
        end
end

function juliaLiteral(e::sumExpr,hc::hiddenCounter = hiddenCounter())

        if e.id == nothing
                hc.counter+=1
                baseName = string(e.name,"Anon",hc.counter)
        else
                baseName = string(e.name,e.id)
        end


        px = juliaLiteral(e.px, hc)

        e.hash = string(e)
        if !haskey(hc.exprList, e.hash)
                eval = baseName
                hc.exprList[e.hash] = (hc.counter, eval, baseName)
                if length(e.px) == 1
                        defs = string(baseName, "=", px.expression)
                else
                        defs = string(baseName,"=","sum(",px.expression,")\n")
                end
                # defs = ""
                # if length(e.px) == 1
                #         eval = px.expression
                # else
                #         eval = string("sum(",px.expression,")")
                # end

        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end

        literal(string(px.definitions,defs), eval, dataAttributes(hc))
end

function juliaLiteral(e::prodExpr, hc::hiddenCounter = hiddenCounter())
        if e.id == nothing
                hc.counter+=1
                baseName = string(e.name,"Anon",hc.counter)
        else
                baseName = string(e.name,e.id)
        end
        px = juliaLiteral(e.px,hc)

        e.hash = string(e)
        if !haskey(hc.exprList, e.hash)
                eval = baseName
                hc.exprList[e.hash] = (hc.counter, eval, baseName)

                if length(e.px) == 1
                        defs = string(baseName,"=", px.expression)

                else
                        defs = string(baseName, "=","prod(",px.expression,")\n")
                end
                # defs = ""
                # if length(e.px) == 1
                #         eval = px.expression
                # else
                #         eval = string("prod(",px.expression,")")
                # end
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end

        literal(string(px.definitions,defs), eval, dataAttributes(hc))
end

function juliaLiteral(e::divExpr,hc::hiddenCounter = hiddenCounter())
        nominator = juliaLiteral(e.nominator,hc)
        denominator = juliaLiteral(e.denominator,hc)

        if e.id == nothing
                hc.counter+=1
                baseName = string(e.name,"Anon",hc.counter)
        else
                baseName = string(e.name,e.id)
        end

        # if e.id == nothing
        #         hc.counter+=1
        #         baseName = string(e.name,"Anon",hc.counter)
        #         defs = ""
        #         eval = string(nominator.expression,"/(",denominator.expression,")")
        # else
        #         baseName = string(e.name,e.id)
        #         defs = string(baseName, "=",nominator.expression,"/(",denominator.expression,")\n" )
        #         eval = baseName
        # end

        e.hash = string(e)
        if !haskey(hc.exprList, e.hash)
                eval = baseName
                hc.exprList[e.hash] = (hc.counter, eval, baseName)
                defs = string(baseName, "=", nominator.expression,"/(",denominator.expression, ")\n")
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end

        literal(string(nominator.definitions, denominator.definitions, defs), eval, dataAttributes(hc))
end

function juliaLiteral(e::dotProdExpr,hc::hiddenCounter = hiddenCounter())
        # baseName = string(e.name,e.id)
        # pos = string(baseName,"pos")
        # par = string(baseName,"par")
        # defs = string(pos, "=", string(e.pos), "\n",
        #                 par,"=", string(e.par), "\n")
        # eval = string("dot(",pos,",",par,")")
        # return literal(defs, eval)
        if e.id == nothing
                hc.counter+=1
                baseName = string(e.name,"Anon",hc.counter)
        else
                baseName = string(e.name,e.id)
        end
        pos = juliaLiteral(e.pos,hc)
        par = juliaLiteral(e.par,hc)

        e.hash = string(e)
        if !haskey(hc.exprList, e.hash)
                eval = baseName
                hc.exprList[e.hash] = (hc.counter, eval, baseName)
                defs = string(baseName,"=","dot(",pos.expression,",",par.expression,")\n")
                # defs = ""
                # if length(e.pos) == 1
                #         eval = string(pos.expression,"*",par.expression)
                # else
                #         eval = string("dot(",pos.expression,",",par.expression,")")
                # end
        else
                defs = ""
                eval = hc.exprList[e.hash].eval
        end


        literal(string(pos.definitions, par.definitions,defs), eval, dataAttributes(hc))
end

function neg(e::T) where T <: ScalarDataIndex
        return simplifyExpr(prodExpr([constIndex(0, "negation", -1.0), e]))
end

function members(e::dataIndex)
        return []
end

function members(e::constIndex)
        return []
end

function members(e::sumExpr)
        return e.px
end

function members(e::prodExpr)
        return e.px
end

function members(e::divExpr)
        return [e.nominator, e.denominator]
end

function members(e::dotProdExpr)
        return vcat(e.pos, e.par)
end

function deRef(e::dataIndex)
        return deepcopy(e)
end

function deRef(e::constIndex)
        return deepcopy(e)
end

function deRef(e::prodExpr)
        pr = prodExpr(ScalarDataIndex[])
        pr.px = deRef.(e.px)
        pr.id = deepcopy(e.id)
        pr.hash = deepcopy(e.hash)
        pr.name = deepcopy(e.name)
        pr.keepAsBlock = e.keepAsBlock
        return pr
end

function deRef(e::sumExpr)
        pr = sumExpr(ScalarDataIndex[])
        pr.px = deRef.(e.px)
        pr.id = deepcopy(e.id)
        pr.hash = deepcopy(e.hash)
        pr.name = deepcopy(e.name)
#        pr.keepAsBlock = e.keepAsBlock
        return pr
end

function deRef(e::divExpr)
        pr = divExpr(constIndex(0, "dummy", 1.0), constIndex(0, "dummy", 1.0))
        #pr = divExpr( deRef(e.nominator), deRef(e.denominator),deepcopy(e.id) )
        pr.nominator = deRef(e.nominator)
        pr.denominator = deRef(e.denominator)
        pr.id = deepcopy(e.id)
        pr.name = deepcopy(e.name)
        pr.hash = deepcopy(e.hash)
        return pr
end

function deRef(e::dotProdExpr)
        pr = dotProdExpr(ScalarDataIndex[], ScalarDataIndex[])
        pr.pos = deRef.(e.pos)
        pr.par = deRef.(e.par)
        pr.id = deepcopy(e.id)
        pr.name = deepcopy(e.name)
        pr.hash = deepcopy(e.hash)
        return pr
end
# mutable struct absExpr <: ScalarDataIndex
#         id::Union{Int, Nothing}
#         name::String
#         pos::ScalarDataIndex
#         hash::Union{Nothing,String}
#
#         function absExpr(pos::ScalarDataIndex, id::Union{Int, Nothing} = nothing)
#                 new(id, "absExpr", pos, nothing)
#         end
# end
#
# #abs(const) can be computed
# #abs(prod) = prod(abs) if needed
# #abs(gauss) = gauss
# #abs(pdf) = pdf
# #abs(cdf) = cdf, usually, but since there's a const/nonconst offset, this is not true generally
# function simplifyExpr(ae::absExpr)
#         pos = simplifyExpr(ae.pos)
#         if typeof(pos) == constIndex
#                 return constIndex(0, "absExprConst", abs(pos.value))
#         end
#         return absExpr(simplifyExpr(ae.pos))
# end
#
# function sameTypeIsLess(a::absExpr, b::absExpr)
#         if a.id!=b.id
#                 return a.id < b.id
#         else
#                 return string(a) < string(b)
#         end
# end
#
# exprOrder[absExpr] = 13
#
#
# function reorder(ae::absExpr)
#         reorder(ae.pos)
#         ae.hash = string(ae)
#         return nothing
# end
#
# function deRef(e::absExpr)
#         pr = absExpr(  constIndex(0, "dummy", 1.0) )
#         pr.pos = deRef(e.pos)
#
#         pr.id = deepcopy(e.id)
#         pr.name = deepcopy(e.name)
#         pr.hash = deepcopy(e.hash)
#         return pr
# end
