function simpson_index(sample_array::Array{Int64,1})
    numerator = map(n->n*(n-1),sample_array) |> sum
    denominator = sum(sample_array)*(sum(sample_array)-1)
    return numerator/denominator
end

function simpson_index(OTU_table::DataFrame,taxonomic_level::Symbol=:OTU)
    sample_names = filter(name->name!=taxonomic_level,names(OTU_table))
    simpson_indices = map(sample_name -> simpson_index(convert(Array{Int64,1},OTU_table[sample_name])),sample_names)
    summarystats(simpson_indices) |> println
    return simpson_indices
end

function shannon_index(sample_array::Array{Int64,1})
    N = sum(sample_array)
    sample_array_nozeroes = filter(x->x!=0,sample_array)
    return map(n->(n/N)*log(n/N),sample_array_nozeroes) |> sum |> -
end

function shannon_index(OTU_table::DataFrame,taxonomic_level::Symbol=:OTU)
    sample_names = filter(name->name!=taxonomic_level,names(OTU_table))
    shannon_indices = map(sample_name -> shannon_index(convert(Array{Int64,1},OTU_table[sample_name])),sample_names)
    summarystats(shannon_indices) |> println
    return shannon_indices
end

function richness_summary(OTU_table::DataFrame,taxonomic_level::Symbol=:OTU)
    sample_names = filter(name->name!=taxonomic_level,names(OTU_table))
    richnesses = map(sample_name -> length(OTUs_in_sample(OTU_table,sample_name)),sample_names)
    summarystats(richnesses) |> println
    return richnesses
end
