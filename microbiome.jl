using DataFrames
using DataStructures

type microbiome_data
    otu::DataFrame
    tax::DataFrame
    meta::DataFrame
    seq::DataFrame
end

source_dir = dirname(Base.source_path())

#Input/output
joinpath(source_dir, "mothur_import.jl") |> include
joinpath(source_dir, "mothur_export.jl") |> include
joinpath(source_dir, "qiime_export.jl") |> include
joinpath(source_dir, "dada2_import.jl") |> include


#Plumbing
joinpath(source_dir, "OTU_filtering.jl") |> include
joinpath(source_dir, "sample_filtering.jl") |> include
joinpath(source_dir, "OTU_combining.jl") |> include
joinpath(source_dir, "OTUs_in_sample.jl") |> include

#Statistics
joinpath(source_dir, "alpha_diversity_indices.jl") |> include
joinpath(source_dir, "rarefaction.jl") |> include
joinpath(source_dir, "beta_diversity.jl") |> include

#Plotting
joinpath(source_dir, "plotting.jl") |> include
#joinpath(source_dir, "deseq2_processing.jl") |> include

