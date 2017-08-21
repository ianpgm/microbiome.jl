function make_deseq2_significant_taxa_table(taxon_table::DataFrame,deseq2_results_table::DataFrame,taxonomic_level::Symbol;numerator_condition="Numerator Condition",denominator_condition="Denominator Condition",significance_level=0.05)
    minimal_taxon_table = taxon_table[[:OTU,taxonomic_level]]
    significant_results = deseq2_results_table[deseq2_results_table[:padj] .< significance_level,:]
    significant_results = join(significant_results,minimal_taxon_table,on=:OTU,kind=:left)
    
    significant_positive_results = significant_results[significant_results[:log2FoldChange] .> 0.0,:]
    significant_positive_counter = counter(significant_positive_results[taxonomic_level])
    positive_taxon_table = DataFrame()
    positive_taxon_table[taxonomic_level] = significant_positive_counter |> keys |> collect
    positive_taxon_table[:num_positive_OTUs] = significant_positive_counter |> values |> collect

    significant_negative_results = significant_results[significant_results[:log2FoldChange] .< 0.0,:]
    significant_negative_counter = counter(significant_negative_results[taxonomic_level])
    negative_taxon_table = DataFrame()
    negative_taxon_table[taxonomic_level] = significant_negative_counter |> keys |> collect
    negative_taxon_table[:num_negative_OTUs] = significant_negative_counter |> values |> collect

    complete_significant_taxa_table = join(positive_taxon_table,negative_taxon_table,on=taxonomic_level,kind=:outer)
    complete_significant_taxa_table[:num_negative_OTUs] = NA_to_zero(complete_significant_taxa_table[:num_negative_OTUs])
    complete_significant_taxa_table[:num_positive_OTUs] = NA_to_zero(complete_significant_taxa_table[:num_positive_OTUs])

    rename!(complete_significant_taxa_table,:num_positive_OTUs,Symbol(numerator_condition))
    rename!(complete_significant_taxa_table,:num_negative_OTUs,Symbol(denominator_condition))

    return complete_significant_taxa_table
end
    
function make_deseq2_significant_taxa_figure(taxon_table::DataFrame,deseq2_results_table::DataFrame,taxonomic_level::Symbol;numerator_condition="Numerator Condition",denominator_condition="Denominator Condition",significance_level=0.05)
    significant_taxa_table = make_deseq2_significant_taxa_table(taxon_table,deseq2_results_table,taxonomic_level;numerator_condition=numerator_condition,denominator_condition=denominator_condition,significance_level=significance_level)
    groupedbar(convert(Array,significant_taxa_table[taxonomic_level]),convert(Array{Real,2},significant_taxa_table[[Symbol(numerator_condition),Symbol(denominator_condition)]]),orientation=:vertical,bar_position=:dodge, label=["Greater in "*numerator_condition,"Greater in "*denominator_condition])
    yaxis!("Number of OTUs")
    xaxis!(String(taxonomic_level))
end
