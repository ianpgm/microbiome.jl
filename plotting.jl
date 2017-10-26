using Plots
using StatPlots

#gr()

function stacked_barplot(taxon_df::DataFrame,taxonomic_level::Symbol)
    sample_names = filter(name->name!=taxonomic_level,names(taxon_df))
    groupedbar(map(string,sample_names),convert(Array{Real,2},taxon_df[sample_names])',bar_position=:stack,label=convert(Array,taxon_df[taxonomic_level]))
end

function fractionate_table(taxon_df::DataFrame,taxonomic_level::Symbol)
    taxon_matrix = convert(Array{Real,2},taxon_df[filter(name->name!=taxonomic_level,names(taxon_df))])
    sum_array = reshape(convert(Array{Real,1},map(sample_name->sum(taxon_df[sample_name]),filter(name-> name!=taxonomic_level,names(taxon_df)))),(1,length(names(taxon_df))-1))
    output_df = DataFrame(taxon_matrix./sum_array)
    for (old_name,new_name) in zip(names(output_df),filter(name->name!=taxonomic_level,names(taxon_df)))
        rename!(output_df,old_name,new_name)
    end
    output_df[taxonomic_level] = taxon_df[taxonomic_level]

    return output_df
end