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

function bray_curtis(sample_i::Array{Int64,1},sample_j::Array{Int64,1})
    Si = sum(sample_i)
    Sj = sum(sample_j)
    Cij = 0
    for (i,j) in zip(sample_i,sample_j)
        if i > 0 && j > 0
            if i <= j
                Cij = Cij + i
            elseif j < i
                Cij = Cij + j
            end
        end
    end

    BC = 1-(2*Cij)/(Si+Sj)

    return BC
end
