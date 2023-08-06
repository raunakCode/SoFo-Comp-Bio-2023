using CSV, DataFrames, NextGenSeqUtils, Distances, LSHFunctions, Graphs, Plots

cd("C:\\Users\\Admin\\Desktop\\programming\\julia stuff\\SoFo 2023\\Human data")

df = DataFrame(CSV.File("sample_2-4-5-7_filtered_assignments_100k_members.tsv", limit = 100000, delim = "\t"))
dropmissing!(df)
df[!,"cdr3_len"] = [length(cdr3) for cdr3 in df[!,"cdr3"]]
df[!,"clone_id_sl"] = df[!,"clone_id"]
df[!,"clone_id"] .= -1

function oneHotEncode(s::String) 
    bitVec = zeros(Bool, 4*length(s))
    pos = 1
    for nuc in s
        if nuc == 'A'
            bitVec[pos] = 1
        elseif nuc == 'C'
            bitVec[pos+1] = 1
        elseif nuc == 'G'
            bitVec[pos+2] = 1
        elseif nuc == 'T'
            bitVec[pos+3] = 1
        end
        pos += 4
    end
    return bitVec
end

b = 5
r = 5

function localitySH(onehots)
    bands = b
    nrows = r
    n_hash = bands * nrows
    clusterGraph = path_graph(0)
    add_vertices!(clusterGraph, length(onehots))
    ids = [-1 for i in 1:length(onehots)]
    hashfn1 = LSHFunction(ℓ1, n_hash)
    initialHashvals = hashfn1.([el[2] for el in onehots])
    finalHashvals = [Dict{Int32, Vector{Int64}}() for i in 1:bands]
    curBand = 0
    for i in 1:nrows:n_hash
        curBand += 1
        curHashfn = LSHFunction(ℓ1, 1)
        for seq in 1:length(onehots)
            seqBand = Vector{Int32}()
            for j in i:i+nrows-1
                push!(seqBand, initialHashvals[seq][j]) 
            end
            curHash = curHashfn(seqBand)
            push!(get!(finalHashvals[curBand], curHash[1], Vector{Int64}()), seq)
        end
    end

    for map in finalHashvals
        for band in map 
            if length(band.second) <= 1
                continue
            end
            for i in 2:length(band.second)
                add_edge!(clusterGraph, band.second[i], band.second[1])
            end
        end
    end  

    components = connected_components(clusterGraph)
    for i in 1:length(components)
        for seq in components[i]
            ids[seq] = i
        end
    end
    return ids, length(components)
end

function assign_clone_ids(sub, NNFunc::Function)
    oneHots = NamedTuple{(:id, :oneHot), Tuple{Int64, Vector{Bool}}}[]
    for i in 1:length(sub.cdr3)
        push!(oneHots, (id = i, oneHot = oneHotEncode(sub.cdr3[i])))
    end
    clone_ids, clusters = NNFunc(oneHots)
    return clone_ids, clusters
end

biggroup = df[(df.v_call .== "IGHV2-5*02") .& (df.j_call .== "IGHJ5*02") .& (df.cdr3_len .== 51),:]
using Clustering

scaleToCluster = Vector{Vector{Int64}}(undef, 0)
x = Vector{Int64}()

@time ids, clusters = assign_clone_ids(biggroup, localitySH)
RI = randindex(biggroup.clone_id_sl, ids)

minBiggroup = minimum(biggroup.clone_id_sl)
biggroupClusters = [Set{Int64}() for i in 1:(maximum(biggroup.clone_id_sl)-minimum(biggroup.clone_id_sl)+1)]
lshClusters = [Set{Int64}() for i in 1:maximum(ids)]

for i in 1:length(biggroup.clone_id_sl)
    push!(biggroupClusters[biggroup.clone_id_sl[i] - minBiggroup + 1], i)
end

for i in 1:length(ids)
    push!(lshClusters[ids[i]], i)
end

sort(biggroupClusters, by=length, rev = true)
sort(lshClusters, by=length, rev = true)

comparison = zeros(Float64, length(biggroupClusters), length(lshClusters))

for i in 1:length(biggroupClusters) 
    for j in 1:length(lshClusters)
        comparison[i, j] = LSHFunctions.jaccard(biggroupClusters[i], lshClusters[j])
        counts = length(intersect(biggroupClusters[i], lshClusters[j]))
    end
end

println("#Clusters: ", clusters)
println("randIndex: ", RI)

plot(heatmap(comparison), xlabel = "LSH Clusters", ylabel = "Single linkage clusters", title = "Rows: $r | Bands: $b", 
    xticks = 1:1:length(lshClusters), yticks = 1:1:length(biggroupClusters)
)

cnts = zeros(Float64, length(biggroupClusters), length(lshClusters))

for i in 1:length(biggroupClusters) 
    for j in 1:length(lshClusters)
        cnts[i, j] = length(intersect(biggroupClusters[i], lshClusters[j]))
    end
end

plot(x, y, title = "Onehot: hash scale to #clusters", ylabel = "#Clusters", xlabel= "Hash scale")
