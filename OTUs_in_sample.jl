function OTUs_in_sample(OTU_table::DataFrame, sample::Symbol)
    return convert(Array{String,1},OTU_table[convert(Array{Bool,1},OTU_table[sample] .> 0),:OTU])
end

function shared_OTUs(OTU_table::DataFrame,samples::Array{Symbol,1})
    OTUs_in_each_sample = Array{String,1}[]
    for sample in samples
        push!(OTUs_in_each_sample,OTUs_in_sample(OTU_table,sample))
    end
    return intersect(OTUs_in_each_sample...)
end

function unique_OTUs(OTU_table::DataFrame,sample::Symbol)
    OTUs_in_target_sample = OTUs_in_sample(OTU_table,sample)
    OTUs_in_other_samples = Set(String[])
    for other_sample in names(OTU_table)
        if other_sample == :OTU || other_sample == sample
            continue
        else
            union!(OTUs_in_other_samples,Set(OTUs_in_sample(OTU_table,other_sample)))
        end
    end
    return setdiff(OTUs_in_target_sample,OTUs_in_other_samples)
end
