using CSV
using GZip
using JSON
using FASTX
using IterTools
using CodecZlib
using DataFrames
using Combinatorics

if length(ARGS) != 2
    println("Usage: julia readfastq.jl Spatial.csv rowindex")
    @assert false
end

################################################################################
##### READ THE SHEET ###########################################################
################################################################################

csv = ARGS[1]                  ; println("csv: " * csv) 
rowindex = parse(Int, ARGS[2]) ; println("rowindex: " * string(rowindex))

row = CSV.read(csv, DataFrame, header=false)[rowindex, :]

RNAINDEX = row[1]       ; println("RNAINDEX: " * RNAINDEX)
SBINDEX = row[2]        ; println("SBINDEX: " * SBINDEX)
GSFASTQS = row[3]       ; println("GSFASTQS: " * GSFASTQS)
GSPUCK = row[4]         ; println("GSPUCK: " * GSPUCK)
GSCB = row[5]           ; println("GSCB: " * GSCB)
GSCBFULL = row[6]       ; println("GSCBFULL: " * GSCBFULL)
GSINTRON = row[7]       ; println("GSINTRON: " * GSINTRON)
GSNOINTRON = row[8]     ; println("GSNOINTRON: " * GSNOINTRON)

################################################################################
##### CONVERT GS PATHS TO LOCAL PATHS ##########################################
################################################################################

run(`mkdir -p FASTQS`)
run(`mkdir -p PUCK`)
run(`mkdir -p CB`)
run(`mkdir -p CBFULL`)
run(`mkdir -p INTRON`)
run(`mkdir -p NOINTRON`)

run(`gcloud storage cp "$GSFASTQS/**/$SBINDEX*_R*.fastq.gz" FASTQS`)
run(`gcloud storage cp "$GSPUCK" PUCK`)
run(`gcloud storage cp "$GSCB" CB`)
run(`gcloud storage cp "$GSCBFULL" CBFULL`)
run(`gcloud storage cp "$GSINTRON" INTRON`)
run(`gcloud storage cp "$GSNOINTRON" NOINTRON`)

################################################################################
##### LOAD LOCAL FILES #########################################################
################################################################################

# Define the paths
fastqs = readdir("FASTQS",join=true)       ; println(fastqs)
sbspath = "PUCK/"*basename(GSPUCK)         ; println(sbspath)
cbspath = "CB/"*basename(GSCB)             ; println(cbspath)
fullcbspath = "CBFULL/"*basename(GSCBFULL) ; println(fullcbspath)
cb_dict_path = "3M-february-2018.txt"      ; println(cb_dict_path)
output_path = RNAINDEX                     ; println(output_path)

# Process and Load
R1s = filter(s -> occursin("_R1_", s), fastqs) ; println("R1s: ", R1s)
R2s = filter(s -> occursin("_R2_", s), fastqs) ; println("R2s: ", R2s, "\n")
@assert length(R1s) == length(R2s)

sbdf = CSV.read(sbspath, DataFrame, header=false)
rename!(sbdf, [:sb, :x, :y])
sbs = Set(sbdf.sb)
println("Loaded $sbspath: $(length(sbs)) sbs found")
println("First few spatial barcodes: ", collect(sbs)[1:3], "\n")
@assert all(length(s) == 14 for s in sbs)

function loadCBs(cbspath)
    # Read file (only .tsv and .tsv.gz are recognized, can add more below)
    if length(cbspath)>=7 && cbspath[end-6:end]==".tsv.gz"
        cbs = readlines(GZip.open(cbspath))
    elseif length(cbspath)>=4 && cbspath[end-3:end]==".tsv"
        cbs = readlines(cbspath)
    else
        error("Unrecognized input file type: $(cbspath)")
    end
    println("loaded $(cbspath): $(length(cbs)) cbs found")
    # Trim off trailing -1
    if length(cbs[1])==18 && cbs[1][end-1:end]=="-1"
        println("trimming -1 from cell barcodes")
        cbs = [s[1:end-2] for s in cbs]
    end
    println("first few cell barcodes: ", cbs[1:5], "\n")
    @assert all(length(s) == 16 for s in cbs)
    cbs = Set(cbs)
    return(cbs)
end
cbs = loadCBs(cbspath)
fullcbs = loadCBs(fullcbspath)

dict = CSV.File(cb_dict_path, delim='\t', header=false) |> DataFrame
rename!(dict, ["GEX", "FEATURE"])
CB10X = Set(dict.GEX)
gex2feature = Dict(zip(dict.GEX, dict.FEATURE))
feature2gex = Dict(zip(dict.FEATURE, dict.GEX))

GGseq = "GGGGGGGGGGGGGGGGGG"
UPseq = "TCTTCAGCGTTCCCGAGA"
const GG_HD = 3
const UP_HD = 3

################################################################################
##### LEARN REMAP AND R1/R2 ####################################################
################################################################################

switch = false
remap = false

# Switch R1 and R2 if needed
R1len = R1s[1] |> open |> GzipDecompressorStream |> FASTQ.Reader |> first |> FASTQ.sequence |> length
R2len = R2s[1] |> open |> GzipDecompressorStream |> FASTQ.Reader |> first |> FASTQ.sequence |> length
if (R1len < 32) & (R2len < 32)
    error("Unexpected read structure: one fastq should have at least 32 characters")
elseif (R1len < 32) & (R2len >= 32)
    switch = false
elseif (R1len >= 32) & (R2len < 32)
    switch = true
else
    println("Learning the correct R1 and R2 assignment")
    s1 = 0 ; s2 = 0
    i1 = R1s[1] |> open |> GzipDecompressorStream |> FASTQ.Reader
    i2 = R2s[2] |> open |> GzipDecompressorStream |> FASTQ.Reader
    for (i,record) in enumerate(zip(i1, i2))
        i > 100000 ? break : nothing
        global s1 += FASTQ.sequence(record[1])[9:26] == UPseq
        global s2 += FASTQ.sequence(record[2])[9:26] == UPseq
    end
    println("R1: ", s1, " R2: ", s2)
    switch = s2 > s1 ? false : true
end
if switch == true
    println("Switching R1 and R2")
    temp = R1s
    R1s = R2s
    R2s = temp
end

# Remap CB if needed
if all(x -> x in CB10X, cbs)
    println("Learning the remap status")
    cbs1 = Set(cbs) ; cbs2 = Set((k->gex2feature[k]).(cbs))
    c1 = 0 ; c2 = 0
    for (i,record) in R1s[1] |> open |> GzipDecompressorStream |> FASTQ.Reader |> enumerate
        i > 1000000 ? break : nothing
        global c1 += in(FASTQ.sequence(record)[1:16], cbs1)
        global c2 += in(FASTQ.sequence(record)[1:16], cbs2)
    end
    println("No remap: ", c1, " Remap: ", c2)
    remap = c1 > c2 ? false : true
else
    println("Some CB barcodes are not in the 10X dictionary - skipping remap")
    remap = false
end
if remap == true
    println("Remapping CB")
    cbs = Set((k->gex2feature[k]).(cbs))
    fullcbs = Set((k->gex2feature[k]).(fullcbs))
else
    println("Skipping CB remapping")
end

################################################################################
##### PREPARE DICTIONARIES FOR FAST PROCESSING #################################
################################################################################

println("Creating Dictionaries")

@assert isa(cbs, Set)
@assert isa(fullcbs, Set)
@assert isa(sbs, Set)
CBlength = cbs |> first |> length
SBlength = sbs |> first |> length
@assert all(x -> length(x) == CBlength, cbs)
@assert all(x -> length(x) == CBlength, fullcbs)
@assert all(x -> length(x) == SBlength, sbs)

charlist=['A','C','G','T','N']
function listHDneighbors(str, hd)
    @assert occursin(r"^[AGCTN]*$", str)
    @assert all(char -> isa(char, Char), charlist)
    @assert 0 <= hd <= length(str)
    res = Set()
    for inds in combinations(1:length(str), hd)
        chars = [str[i] for i in inds]
        pools = [setdiff(charlist, [char]) for char in chars]
        prods = product(pools...)
        @assert length(prods) == (length(charlist)-1)^hd
        for prod in prods
            s = str
            for (i, c) in zip(inds, prod)
                s = s[1:i-1]*string(c)*s[i+1:end]
            end
            push!(res,s)
        end
    end
    @assert length(res) == length(combinations(1:length(str), hd))*((length(charlist)-1)^hd)
    return(res)
end

function listLD1neighbors(str)
    # Substitutions
    res = listHDneighbors(str,1)
    # Insertions    
    [[push!(res, string(str[1:i], char, str[i+1:end])) for char in charlist] for i in 0:length(str)]
    # Deletions
    [push!(res, string(str[1:i-1], str[i+1:end])) for i in 1:length(str)]
    return(res)
end

# Create sets for fast lookup
GGset = reduce(union, [listHDneighbors(GGseq,i) for i in 0:GG_HD])
UPset = reduce(union, [listHDneighbors(UPseq,i) for i in 0:UP_HD])
UPset = union(UPset, listLD1neighbors(UPseq))

CB1dict = Dict()
[[push!(get!(CB1dict,cb_fuzzy,Set()), cb) for cb_fuzzy in listHDneighbors(cb,1)] for cb in cbs]
CB1set = Set(keys(CB1dict))

SB1dict = Dict()
[[push!(get!(SB1dict,sb_fuzzy,Set()), sb) for sb_fuzzy in listLD1neighbors(sb)] for sb in sbs]
SB1set = Set(keys(SB1dict))

################################################################################
##### PROCESS THE FASTQS #######################################################
################################################################################

println("Reading FASTQs")

function CBprocess(cb)
    if in(cb,cbs)
        cbmatchtype = "exact"
    elseif in(cb,CB1set)
        matches = CB1dict[cb]
        cbmatchtype = (length(matches) == 1) ? "fuzzy" : "fuzzy.many"
        cb = first(matches)
    elseif in(cb,fullcbs)
        cbmatchtype = "uncalled"
    else
        cbmatchtype = "none"
        # CBnone[cb] = get(CBnone,cb,0) + 1
    end
    return cb, cbmatchtype
end

function UPprocess(r2)
    up = r2[9:26]; sb = r2[[1:8; 27:32]]
    if in(up, GGset) # all G
        upmatchtype = "allG"
    elseif up == UPseq # exact match
        upmatchtype = "exact"
    elseif in(up,UPset)
        upmatchtype = "fuzzy"
    elseif in(r2[8:25],UPset) # deletion before
        upmatchtype = "fuzzy"
        sb = r2[[1:7; 26:31]]
    elseif in(r2[9:25],UPset) # deletion inside
        upmatchtype = "fuzzy"
        sb = r2[[1:8; 26:31]]
    else
        upmatchtype = "none"
        # UPnone[up] = get(UPnone,up,0) + 1
    end
    return sb, upmatchtype
end

function SBprocess(sb)
    if in(sb,sbs)
        sbmatchtype = "exact"
    elseif in(sb,SB1set)
        matches = SB1dict[sb]
        sbmatchtype = (length(matches) == 1) ? "fuzzy" : "fuzzy.many"
        sb = first(matches)
    else
        sbmatchtype = "none"
        # SBnone[sb] = get(SBnone,sb,0) + 1
    end
    return sb, sbmatchtype
end

mat = Dict()
metrics = Dict(cbmatchtype=>Dict(upmatchtype=>Dict("exact"=>0,"fuzzy"=>0,"fuzzy.many"=>0,"none"=>0,"noUP"=>0) for upmatchtype in ["exact","fuzzy","allG","none"]) for cbmatchtype in ["exact","fuzzy","fuzzy.many","uncalled","none"])
# CBnone = Dict() ; UPnone = Dict() ; SBnone = Dict()
reads = 0

for fastqpair in zip(R1s,R2s)
    it1 = fastqpair[1] |> open |> GzipDecompressorStream |> FASTQ.Reader;
    it2 = fastqpair[2] |> open |> GzipDecompressorStream |> FASTQ.Reader;

    for record in zip(it1, it2)
        r1 = record[1] |> FASTQ.sequence
        r2 = record[2] |> FASTQ.sequence
        global reads += 1

        cb, cbmatchtype = CBprocess(r1[1:16])
        sb, upmatchtype = UPprocess(r2)
        sb, sbmatchtype = in(upmatchtype,["exact","fuzzy"]) ? SBprocess(sb) : (sb,"noUP")
        
        metrics[cbmatchtype][upmatchtype][sbmatchtype] += 1

        umi = r1[17:28]
        if in(cbmatchtype,["exact","fuzzy"]) & in(sbmatchtype,["exact","fuzzy"])
            if haskey(mat, (cb,sb))
                mat[(cb,sb)][umi] = get(mat[(cb,sb)],umi,0) + 1
            else
                mat[(cb,sb)] = Dict(umi => 1)
            end
        end
    end
end

println(reads)
println(metrics)

vals = []
[[[[push!(vals,val)] for val in values(d2)] for d2 in values(d1)] for d1 in values(metrics)];
@assert sum(vals) == reads

################################################################################
##### SAVE OUTPUT ##############################################################
################################################################################

# Create a dataframe
cb_vec = String[]
sb_vec = String[]
umi_vec = Int[]
reads_vec = Int[]
for (key, dict) in mat
    push!(cb_vec, key[1])
    push!(sb_vec, key[2])
    push!(umi_vec, dict|>keys|>length)
    push!(reads_vec, dict|>values|>sum)
end
df = DataFrame(cb = cb_vec, sb = sb_vec, umi = umi_vec, reads = reads_vec)
df = leftjoin(df, sbdf, on=:sb)
sort!(df, :umi, rev=true)
remap ? transform!(df, :cb => (x -> (k->feature2gex[k]).(x)) => :cb) : nothing

# Save the outout
run(`mkdir -p $output_path`)
CSV.write(joinpath(output_path,"coords.csv"), df)
open(joinpath(output_path,"metrics.json"), "w") do f write(f, JSON.json(metrics)) end
params = [("RNA index",RNAINDEX),
          ("SB index",SBINDEX),
          ("CB path",GSCB),
          ("CBfull path",GSCBFULL),
          ("SB path",GSPUCK),
          ("CB num",length(cbs)),
          ("full CB num",length(fullcbs)),
          ("SB num",length(sbs)),
          ("switch",switch),
          ("remap",remap),
          ("reads",reads),
          ("fastqs",fastqs)]
params = DataFrame(col1 = [x[1] for x in params], col2 = [x[2] for x in params])
CSV.write(joinpath(output_path,"params.csv"), params, header = false)

# @assert isfile("coords.csv")
# @assert isfile("metrics.json")
# @assert isfile("params.csv")
# @assert isdir("INTRON")
# @assert isdir("NOINTRON")

#sort(DataFrame(CB=collect(keys(CBnone)),count=collect(values(CBnone))), :count, rev=true)
#sort(DataFrame(UP=collect(keys(UPnone)),count=collect(values(UPnone))), :count, rev=true)
#sort(DataFrame(SB=collect(keys(SBnone)),count=collect(values(SBnone))), :count, rev=true)