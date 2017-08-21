function make_mothur_shared_df(OTU_table::DataFrame; label::AbstractString = "0.03")
    shared_df = DataFrame()
    shared_df[:label] = fill(label,size(OTU_table)[2]-1)
    shared_df[:Group] = filter(name->name!=:OTU,names(OTU_table))
    shared_df[:numOtus] = fill(size(OTU_table)[1],size(OTU_table)[2]-1)
    rotated_OTU_table = convert(Array,OTU_table[filter(name->name!=:OTU,names(OTU_table))])' |> DataFrame
    shared_df = hcat(shared_df,rotated_OTU_table)
    for (old_name,new_name) in zip(names(shared_df)[4:end],map(Symbol,OTU_table[:OTU]))
        rename!(shared_df, old_name, new_name)
    end
    return shared_df
end

function make_mothur_taxonomy_df(tax_table::DataFrame)
    taxonomy_column = []
    for row_number in 1:size(tax_table)[1]
        push!(taxonomy_column,join(convert(Array,tax_table[row_number,3:end]),";")*";")
    end
    output_df = tax_table[[:OTU,:Size]]
    output_df[:Taxonomy] = taxonomy_column
    return output_df
end

function export_mothur_data(microbiome::microbiome_data,file_prefix::AbstractString)
    shared_table = make_mothur_shared_df(microbiome.otu)
    taxonomy_table = make_mothur_taxonomy_df(microbiome.tax)
    writetable_noquotes(file_prefix*".shared",shared_table)
    writetable_noquotes(file_prefix*".cons.taxonomy",taxonomy_table)
end

function writetable_noquotes(filename::AbstractString,input_df::DataFrame)
    writetable("temp_table",input_df, separator='\t')
    write(filename,map(line->replace(line,'"',"")*"\n",readlines("temp_table")))
    rm("temp_table")
end
