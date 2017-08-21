using DataFrames
using DataStructures

type microbiome_data
    otu::DataFrame
    tax::DataFrame
end

source_dir = dirname(Base.source_path())

#Input/output
joinpath(source_dir, "mothur_import.jl") |> include
joinpath(source_dir, "mothur_export.jl") |> include
joinpath(source_dir, "qiime_export.jl") |> include

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
#joinpath(source_dir, "plotting.jl") |> include
#joinpath(source_dir, "deseq2_processing.jl") |> include

function export_ancom_wide_format(OTU_table::DataFrame,metadata::DataFrame,attribute::Symbol,filename::String)
    rotated_OTU_table = convert(Array,OTU_table[filter(name->name!=:OTU,names(OTU_table))])' |> DataFrame
    renaming_dict = zip(names(rotated_OTU_table),map(Symbol,OTU_table[:OTU]))

    rename!(rotated_OTU_table,renaming_dict)

    rotated_OTU_table[:Sample] = map(String,filter(name->name!=:OTU,names(OTU_table)))
    minimal_metadata = metadata[[:Sample,attribute]]
    ancom_wide_format = join(rotated_OTU_table,minimal_metadata,on=:Sample,kind=:left)
    delete!(ancom_wide_format,:Sample)
    rename!(ancom_wide_format,attribute,:Group)
    writetable_noquotes(filename,ancom_wide_format)
    return ancom_wide_format
end

