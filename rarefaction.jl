function expand_sample(OTU_table::DataFrame,sample::Symbol)
    expanded_OTU_array = String[]
    for (OTU_name, OTU_count) in zip(OTU_table[:OTU],OTU_table[sample])
        append!(expanded_OTU_array,fill(OTU_name,OTU_count))
    end
    return expanded_OTU_array
end

function compress_sample(expanded_sample::Array{String,1})
    OTU_counter = counter(expanded_sample)
    return (map(string,OTU_counter|>keys),map(Int,OTU_counter|>values))
end

function NA_to_zero(input_array::DataArray)
    NA_rows = map(isna,input_array) |> x-> convert(Array{Bool,1},x)
    input_array[NA_rows] = 0
    return input_array
end

function rarefy(OTU_table::DataFrame, reads::Int64, sample::Symbol; taxonomic_level::Symbol=:OTU, remove_small_samples=false)
    samplename = string(sample)
    total_reads = sum(OTU_table[sample]) |> string
    necessary_reads = string(reads)
    output_OTU_table = deepcopy(OTU_table)
    if reads > sum(OTU_table[sample])
        if remove_small_samples
            delete!(output_OTU_table,sample)
            println("$samplename has only $total_reads, needs $necessary_reads, deleting sample as remove_small_samples is set to true")
            return output_OTU_table
        else
            error("$samplename has only $total_reads, needs $necessary_reads")
        end
    end
    expanded_sample = expand_sample(OTU_table,sample)
    rarefied_OTUs,rarefied_OTU_counts = expanded_sample[randperm(length(expanded_sample))[1:reads]] |> compress_sample
    rarefied_sample = DataFrame()
    rarefied_sample[taxonomic_level] = rarefied_OTUs
    rarefied_sample[sample] = rarefied_OTU_counts
    delete!(output_OTU_table,sample)
    output_OTU_table = join(output_OTU_table,rarefied_sample,on=taxonomic_level,kind=:outer)
    output_OTU_table[sample] = output_OTU_table[sample] |> NA_to_zero
    return output_OTU_table
end

function rarefy(OTU_table::DataFrame, reads::Int64; taxonomic_level::Symbol=:OTU, remove_small_samples=false)
    rarefied_OTU_table = deepcopy(OTU_table)
    sample_names = filter(name->name!=taxonomic_level,names(OTU_table))
    for sample in sample_names
        rarefied_OTU_table = rarefy(rarefied_OTU_table, reads, sample; taxonomic_level=taxonomic_level,remove_small_samples=remove_small_samples)
    end
    return rarefied_OTU_table
end

function rarefy(data::microbiome_data, reads::Int64; taxonomic_level::Symbol=:OTU, remove_small_samples=false)
    new_OTU_table = rarefy(data.otu,reads,taxonomic_level=taxonomic_level,remove_small_samples=remove_small_samples)
    new_tax_table = deepcopy(data.tax)
    new_tax_table[:Size] = sum_each_OTU(new_OTU_table)
    return microbiome_data(new_OTU_table,new_tax_table)
end
