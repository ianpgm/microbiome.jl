function sample_subset(OTU_table::DataFrame,samples_to_keep::Array{Symbol,1},taxonomic_level::Symbol)
    columns_to_keep = push!(samples_to_keep,taxonomic_level)
    new_otu_table = OTU_table[columns_to_keep]
    return new_otu_table
end

function remove_samples(OTU_table::DataFrame,samples_to_remove::Array{Symbol,1},taxonomic_level::Symbol)
    samples_to_keep = setdiff(filter(name->name!=taxonomic_level,names(OTU_table)),samples_to_remove)
    return sample_subset(OTU_table,samples_to_keep,taxonomic_level)
end

function sample_subset(data::microbiome_data,samples_to_keep::Array{Symbol,1})
    new_otu_table = sample_subset(data.otu,samples_to_keep,:OTU)
    data.tax[:Size] = sum_each_OTU(new_otu_table)
    if size(data.meta)[1] > 0
        new_metadata = data.meta[convert(Array{Bool,1},map(sample_name -> in(sample_name,convert(Array{String,1},samples_to_keep)),data.meta[:Sample])),:]
    else
        new_metadata = data.meta
    end
    return microbiome_data(new_otu_table,data.tax,new_metadata,data.seq)
end

function remove_samples(data::microbiome_data,samples_to_remove::Array{Symbol,1})
    samples_to_keep = setdiff(filter(name->name!=:OTU,names(data.otu)),samples_to_remove)
    return sample_subset(data,samples_to_keep)
end

function combine_samples(data::microbiome_data,sample_dict::Dict{Symbol,Array{Symbol,1}})
    for combined_sample in keys(sample_dict)
        sub_table = data.otu[sample_dict[combined_sample]]
        data.otu[combined_sample] = sum_each_OTU(sub_table)
        data.tax[:Size] = data.tax[:Size] .+ data.otu[combined_sample]
        data = remove_samples(data,sample_dict[combined_sample])
    end
    data.meta = DataFrame(Sample = convert(Array{String,1},collect(keys(sample_dict))))
    return data
end

function combine_samples(OTU_table::DataFrame,sample_dict::Dict{Symbol,Array{Symbol,1}},taxonomic_level::Symbol)
    for combined_sample in keys(sample_dict)
        sub_table = OTU_table[sample_dict[combined_sample]]
        OTU_table[combined_sample] = sum_each_OTU(sub_table)
        OTU_table = remove_samples(OTU_table,sample_dict[combined_sample],taxonomic_level)
    end
    return OTU_table
end

function group_samples(data::microbiome_data,metadata_field::Symbol;normalize_sample_sizes=true)
    output_data = deepcopy(data)

    if normalize_sample_sizes
        output_data = microbiome_data(fractionate_table(data.otu,:OTU),data.tax,data.meta)
    end

    combine_sample_dict = Dict{Symbol,Array{Symbol,1}}()
    for sample_type in unique(data.meta[metadata_field])
        combine_sample_dict[Symbol(sample_type)] = data.meta[data.meta[metadata_field] .== sample_type,:Sample] |> x -> convert(Array{Symbol,1},x)
    end

    output_data = combine_samples(output_data,combine_sample_dict)

    return output_data
end



