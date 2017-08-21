function combine_OTUs(OTU_table::DataFrame,OTU_dict::Dict{String,Array{String,1}})
    for new_OTU_name in keys(OTU_dict)
        summed_OTU_table_subsection = OTU_table[convert(Array{Bool,1},map(OTU -> in(OTU,OTU_dict[new_OTU_name]),OTU_table[:OTU])),:] |> sample_sums
        newline = DataFrame()
        newline[:OTU] = new_OTU_name
        for (sample_name,sample_sum) in zip(filter(name-> name!=:OTU,names(OTU_table)),summed_OTU_table_subsection)
            newline[sample_name] = sample_sum
        end
        OTU_table = vcat(OTU_table,newline)
        OTU_table = remove_OTUs(OTU_table,OTU_dict[new_OTU_name])
    end
    return OTU_table
end

function collapse_taxon(data::microbiome_data, taxonomic_level::Symbol; other_threshold::Float64=0.0)
    taxa = unique(data.tax[taxonomic_level])
    OTU_dict = Dict{String,Array{String,1}}()
    for taxon in taxa
        OTU_dict[taxon] = convert(Array{String,1},data.tax[data.tax[taxonomic_level] .== taxon,:OTU])
    end
    taxon_df = combine_OTUs(data.otu,OTU_dict) |> combined_df -> rename(combined_df,:OTU,taxonomic_level)
    sum_of_samples = map(sample_name->sum(taxon_df[sample_name]),filter(name-> name!=taxonomic_level,names(taxon_df)))
    other_row = fill(0,length(sum_of_samples))
    output_df = DataFrame()
    for row in eachrow(taxon_df)
        sample_array = filter(x->typeof(x)<:Real,map(x->x[2],row))
        if maximum(sample_array./sum_of_samples) < other_threshold
            other_row = other_row .+ sample_array
        else
            newline = DataFrame(reshape(map(x->x[2],row),(1,length(names(taxon_df)))))
            output_df = vcat(output_df,newline)
        end
    end
    rename!(output_df,Dict(zip(names(output_df),names(taxon_df))))
    if sum(other_row) > 0
        otherline = DataFrame()
        otherline[taxonomic_level] = "Others <"*string(round(other_threshold,2))
        for (name,number) in zip(filter(name-> name!=taxonomic_level,names(taxon_df)),other_row)
            otherline[name] = number
        end
        output_df = vcat(output_df,otherline)
    end
    return output_df
end