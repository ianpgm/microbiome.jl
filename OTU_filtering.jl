function sample_sums(OTU_table::DataFrame)
    taxonomic_levels = [:OTU,:Genus,:Family,:Order,:Class,:Phylum,:Domain,:Kingdom]
    return map(sample_name->sum(OTU_table[sample_name]),filter(name-> in(name,taxonomic_levels)==false,names(OTU_table)))
end

function check_negative_controls(OTU_table::DataFrame,negative_control_sample_names::Array{Symbol,1},multiplication_factor::Float64)
    sample_names_no_negatives = setdiff(filter(name->name!=:OTU,names(OTU_table)),negative_control_sample_names)
    suspect_OTUs = String[]
    for line in eachrow(OTU_table)
        sample_read_counts = map(x->x[2],line[sample_names_no_negatives])
        negative_read_counts = map(x->x[2],line[negative_control_sample_names])
        if maximum(negative_read_counts) > multiplication_factor*maximum(sample_read_counts)
            push!(suspect_OTUs,line[:OTU])
        end
    end
    return suspect_OTUs
end

function OTU_subset(input_data::DataFrame,OTUs_to_keep::Array)
    new_table = DataFrame(OTU = OTUs_to_keep)
    return join(new_table,input_data, on=:OTU, kind=:inner)
end

function OTU_subset(input_data::microbiome_data,OTUs_to_keep::Array)
    output_data = microbiome_data(OTU_subset(input_data.otu,OTUs_to_keep), OTU_subset(input_data.tax,OTUs_to_keep))
    println("Retaining "*string(length(OTUs_to_keep))*" out of "*string(size(input_data.otu)[1])*" OTUs")
    input_sample_sizes = sample_sums(input_data.otu)
    output_sample_sizes = sample_sums(output_data.otu)
    fraction_remaining = output_sample_sizes ./ input_sample_sizes
    println("Retaining the following fraction of reads from samples:")
    println(summarystats(fraction_remaining))
    return output_data
end

function remove_OTUs(input_data::microbiome_data,OTUs_to_remove::Array)
    OTUs_to_keep = setdiff(input_data.otu[:OTU],OTUs_to_remove)
    return OTU_subset(input_data,OTUs_to_keep)
end

function remove_OTUs(input_data::DataFrame,OTUs_to_remove::Array)
    OTUs_to_keep = setdiff(input_data[:OTU],OTUs_to_remove)
    return OTU_subset(input_data,OTUs_to_keep)
end

function sum_matrix_rows(matrix::Array{Real,2})
    sums = Real[]
    for i in 1:size(matrix)[1]
        push!(sums,sum(matrix[i,1:end]))
    end
    return sums
end

function sum_each_OTU(OTU_table::DataFrame)
    OTU_array = convert(Array{Real,2},OTU_table[filter(name->name!=:OTU,names(OTU_table))])
    return sum_matrix_rows(OTU_array)
end

function remove_rare_OTUs(data::microbiome_data, threshold::Int64)
    OTUs_to_keep = data.otu[map(readcount -> readcount > threshold,sum_each_OTU(data.otu)),:OTU]
    return OTU_subset(data,convert(Array,OTUs_to_keep))
end

function identify_OTUs_by_taxon(taxon_table::DataFrame,taxon_dict::Dict{Symbol,Array{String,1}})
    OTUs_identified = Set(String[])
    for taxonomic_level in keys(taxon_dict)
        for taxon_name in taxon_dict[taxonomic_level]
            OTUs_identified = union(OTUs_identified,taxon_table[taxon_table[taxonomic_level] .== taxon_name,:OTU])
        end
    end
    return OTUs_identified
end

function remove_taxon(data::microbiome_data,taxon_dict::Dict{Symbol,Array{String,1}})
    OTUs_to_remove = identify_OTUs_by_taxon(data.tax,taxon_dict)
    return remove_OTUs(data,OTUs_to_remove)
end

function taxon_subset(data::microbiome_data,taxon_dict::Dict{Symbol,Array{String,1}})
    OTUs_to_keep = identify_OTUs_by_taxon(data.tax,taxon_dict)
    return OTU_subset(data,OTUs_to_keep)
end