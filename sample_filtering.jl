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
    return microbiome_data(new_otu_table,data.tax)
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
