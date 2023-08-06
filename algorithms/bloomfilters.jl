using FASTX
using GZip
using CodecZlib
using BloomFilters
using ProfileView
using Plots

cd("C:\\Users\\Admin\\Desktop\\programming\\julia stuff\\SoFo 2023\\Macaques")

BF = []
BloomFilterCount = 18
ORcontaminants = []
finalContaminants = []
numberOfSeq = zeros(Int64, BloomFilterCount) 
seqCount = 0
bitCount = 0
libNames = []

function createFilters(number::Int64) 
    for i in 1:number
        push!(BF, BloomFilter(6964500, 0.01))
    end
end

function addSeq(seqs, names, BFnumber, excludedOR; f_kwargs...)
    # seqCount += length(seqs)
    for i in 1:length(seqs)
        add!(BF[BFnumber], seqs[i])
    end
end

function getOR(excludedFilter::Int64) 
    startBitVector = BloomFilter(6964500, 0.01)
    for i in 1:BloomFilterCount
        if i != excludedFilter
            startBitVector.array .|= BF[i].array
        end
    end
    return startBitVector.array
end

function checkClash(sequences, names, BFnumber::Int64, excludedOR; f_kwargs...) 
    tempBF = BloomFilter(6964500, 0.01)
    tempBF.array .= excludedOR
    for seq in sequences
        if (in(seq, tempBF)) 
            push!(ORcontaminants, (seq, BFnumber))
        end
    end
end

function chunked_fastq_apply(fpath, func::Function, BFnumber::Int64; chunk_size=10000, f_kwargs = [], verbose = false)
    if length(libNames) < BloomFilterCount
        push!(libNames, fpath[5:9])
    end
    if endswith(fpath, ".gz")
        reader = FASTQ.Reader(GzipDecompressorStream(open(fpath)))
    else
        reader = FASTQ.Reader(open(fpath))
    end
    cnt = 0
    seqs, names = [], []
    chunk = 0				
    i = 0
    record = FASTQ.Record()
    excludedOR = getOR(BFnumber)
    # read_counts = [0, 0, 0]
    while !eof(reader)
        read!(reader, record)
        push!(seqs, FASTQ.sequence(String, record))
        push!(names, FASTQ.identifier(record))
        i += 1
        if i == chunk_size
            #apply func...
        chunk += 1						
        cnt += i
            # seqCount += length(seqs)
            # read_counts += func(seqs, names; f_kwargs...)
            func(seqs, names, BFnumber, excludedOR; f_kwargs...)
            seqs, names = [], []
            i = 0
            if verbose print(".") end
        end
    end
    close(reader)
    if i > 0
        #apply func...
    chunk += 1
    cnt += i
        func(seqs, names, BFnumber, excludedOR; f_kwargs...)
    end
    if verbose print(".") end
    numberOfSeq[BFnumber] = cnt
    # return read_counts
end

# create BloomFilters
createFilters(BloomFilterCount)

bitCount = length(BF[1].array)

# add sequences to corresponding Bloom Filters | Stored in BF
chunked_fastq_apply("Mac_Cy1_L_merged.fastq.gz", addSeq, 1) 
chunked_fastq_apply("Mac_Cy2_L_merged.fastq.gz", addSeq, 2) 
chunked_fastq_apply("Mac_Cy3_L_merged.fastq.gz", addSeq, 3)
chunked_fastq_apply("Mac_Cy4_L_merged.fastq.gz", addSeq, 4)
chunked_fastq_apply("Mac_Cy5_L_merged.fastq.gz", addSeq, 5)
chunked_fastq_apply("Mac_Cy6_L_merged.fastq.gz", addSeq, 6)
chunked_fastq_apply("Mac_Cy7_L_merged.fastq.gz", addSeq, 7)
chunked_fastq_apply("Mac_Cy1_U_merged.fastq.gz", addSeq, 8) 
chunked_fastq_apply("Mac_Cy2_U_merged.fastq.gz", addSeq, 9) 

chunked_fastq_apply("Mac_Cy3_U_merged.fastq.gz", addSeq, 10)
chunked_fastq_apply("Mac_Cy4_U_merged.fastq.gz", addSeq, 11)
chunked_fastq_apply("Mac_Cy5_U_merged.fastq.gz", addSeq, 12)
chunked_fastq_apply("Mac_Cy6_U_merged.fastq.gz", addSeq, 13)
chunked_fastq_apply("Mac_Cy7_U_merged.fastq.gz", addSeq, 14)
chunked_fastq_apply("Mac_D11_U_merged.fastq.gz", addSeq, 15)
chunked_fastq_apply("Mac_D16_U_merged.fastq.gz", addSeq, 16)
chunked_fastq_apply("Mac_D19_U_merged.fastq.gz", addSeq, 17)
chunked_fastq_apply("Mac_Rh5_U_merged.fastq.gz", addSeq, 18)

# find sequences that don't pass the OR test | Stored in ORcontaminants
chunked_fastq_apply("Mac_Cy1_L_merged.fastq.gz", checkClash, 1) 
chunked_fastq_apply("Mac_Cy2_L_merged.fastq.gz", checkClash, 2) 
chunked_fastq_apply("Mac_Cy3_L_merged.fastq.gz", checkClash, 3)
chunked_fastq_apply("Mac_Cy4_L_merged.fastq.gz", checkClash, 4)

chunked_fastq_apply("Mac_Cy5_L_merged.fastq.gz", checkClash, 5)
chunked_fastq_apply("Mac_Cy6_L_merged.fastq.gz", checkClash, 6)
chunked_fastq_apply("Mac_Cy7_L_merged.fastq.gz", checkClash, 7)
chunked_fastq_apply("Mac_Cy1_U_merged.fastq.gz", checkClash, 8) 
chunked_fastq_apply("Mac_Cy2_U_merged.fastq.gz", checkClash, 9) 

chunked_fastq_apply("Mac_Cy3_U_merged.fastq.gz", checkClash, 10)
chunked_fastq_apply("Mac_Cy4_U_merged.fastq.gz", checkClash, 11)
chunked_fastq_apply("Mac_Cy5_U_merged.fastq.gz", checkClash, 12)
chunked_fastq_apply("Mac_Cy6_U_merged.fastq.gz", checkClash, 13)
chunked_fastq_apply("Mac_Cy7_U_merged.fastq.gz", checkClash, 14)
chunked_fastq_apply("Mac_D11_U_merged.fastq.gz", checkClash, 15)
chunked_fastq_apply("Mac_D16_U_merged.fastq.gz", checkClash, 16)
chunked_fastq_apply("Mac_D19_U_merged.fastq.gz", checkClash, 17)
chunked_fastq_apply("Mac_Rh5_U_merged.fastq.gz", checkClash, 18)

println("#Sequences failing the OR test: ", sizeof(ORcontaminants))

# check OR-failing sequences against all libraries
for x in ORcontaminants
    seq = x[1]
    possibleOrigins = []
    for i in 1:length(BF)
        if i != x[2] && in(seq, BF[i])
            push!(possibleOrigins, i)
        end
    end
    if length(possibleOrigins) > 0
        push!(finalContaminants, [seq, x[2], possibleOrigins])
    end
end

println("#Final-contaminants: ", sizeof(finalContaminants))
println("#sequences: ", seqCount)

leakage = zeros(Int64, BloomFilterCount, BloomFilterCount)

for i in 1:length(finalContaminants)
    data = finalContaminants[i][2]
    for origin in finalContaminants[i][3]
        leakage[data, origin] += 1
    end
end

hm = Array{Float64}(undef, BloomFilterCount, BloomFilterCount)
for i in 1:BloomFilterCount
    for j in 1:BloomFilterCount
        hm[i, j] = leakage[i, j] / numberOfSeq[i]
    end
end


plot(
    heatmap(hm), 
    title = "Leakage extent", ylabel = "Sequence Destination", xlabel = "Sequence Origin", 
    xrotation = 45, xticks = (1:1:BloomFilterCount, libNames), yticks = (1:1:BloomFilterCount, libNames)
)
