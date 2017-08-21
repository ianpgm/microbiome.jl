function bray_curtis(OTU_table::DataFrame,samples::Array{Symbol,1})
    shared_OTU_array = shared_OTUs(OTU_table,samples)
    shared_rows = map(x->in(x,shared_OTU_array),OTU_table[:OTU]) |> x -> convert(Array{Bool,1},x)
    shared_OTU_table = OTU_table[shared_rows,:]
    C = map(minimum,zip(shared_OTU_table[samples[1]],shared_OTU_table[samples[2]])) |> sum
    Si = sum(OTU_table[samples[1]])
    Sj = sum(OTU_table[samples[2]])
    BC = 1 - ((2*C)/(Si +Sj))
    return BC
end