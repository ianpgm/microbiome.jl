function export_qiime_otu_table(OTU_table::DataFrame,output_filename::String)
    qiime_otu_table = deepcopy(OTU_table)
    rename!(qiime_otu_table,:OTU,Symbol("#OTU ID"))
    new_column_order = [Symbol("#OTU ID")]
    append!(new_column_order,filter(name->name!=Symbol("#OTU ID"),names(qiime_otu_table)))
    qiime_otu_table = qiime_otu_table[new_column_order]
    writetable_noquotes(output_filename,qiime_otu_table)
    new_table = ["#\n"]
    append!(new_table,map(line->line*"\n",readlines(output_filename)))
    write(output_filename,new_table)
end
