function dada2_seqtab_import(seqtab_file::AbstractString)
    seqtab = readtable(seqtab_file,header=true,separator='\t');
    SV_numbers = 1:size(seqtab)[2]-1 |> collect .|> string
    max_SV_ID_length = map(length,SV_numbers) |> maximum
    SV_names = "SV".*map(SV_number -> lpad(SV_number,max_SV_ID_length,"0"),SV_numbers)
    seq_df = DataFrame(OTU = SV_names, Sequence = map(string,names(seqtab)[2:end]))
    #rename!(seqtab,Dict(zip(names(seqtab)[2:end],SV_names)));
    
    OTU_table = DataFrame(transpose(convert(Array,seqtab[2:end])))
    rename!(OTU_table,Dict(zip(names(OTU_table),map(Symbol,seqtab[1]))))
    
    OTU_table[:OTU] = SV_names

    return(OTU_table,seq_df)
end

function dada2_taxa_import(taxa_file::AbstractString,seq_df::DataFrame,OTU_table::DataFrame)
    taxon_table = readtable(taxa_file,separator='\t',header=true)
    rename!(taxon_table,:x,:Sequence)
    rename!(taxon_table,:Kingdom,:Domain)
    taxon_table = join(taxon_table,seq_df,on=:Sequence,kind=:left)
    delete!(taxon_table,:Sequence)
    OTU_sizes = DataFrame(OTU = OTU_table[:OTU],Size = sum_each_OTU(OTU_table))
    taxon_table = join(taxon_table,OTU_sizes,on=:OTU,kind=:left)
    for taxonomic_level in [:Domain,:Phylum,:Class,:Order,:Family,:Genus]
        taxon_table[isna.(taxon_table[taxonomic_level]),taxonomic_level] = string(taxonomic_level)*"_unclassified"
    end
    return taxon_table
end
    
function dada2_import(seqtab_file::AbstractString,taxa_file::AbstractString)
    OTU_table,seq_df = dada2_seqtab_import(seqtab_file)
    taxon_table = dada2_taxa_import(taxa_file,seq_df,OTU_table)
    data = microbiome_data(OTU_table,taxon_table,DataFrame(),seq_df)
    return data
end
