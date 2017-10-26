function import_mothur_shared_file(shared_file::AbstractString)
    OTU_table = readtable(shared_file, header=true, separator='\t')

    OTU_table_rot = DataFrame(transpose(convert(Array,OTU_table[4:end])))
    for (old_name,sample_name) in zip(names(OTU_table_rot),map(Symbol,OTU_table[:Group]))
        rename!(OTU_table_rot,old_name,sample_name)
    end

    OTU_table_rot[:OTU] = map(string,names(OTU_table)[4:end])
    return OTU_table_rot
end

function remove_number_from_taxon_string(taxon_string::AbstractString)
    return split(taxon_string,"(")[1]
end

function split_mothur_taxonomy_string(taxonomy_string::AbstractString)
    return split(taxonomy_string,';')[1:end-1] |> taxa_with_numbers -> map(remove_number_from_taxon_string,taxa_with_numbers)
end

function import_mothur_taxonomy_file(taxonomy_file::AbstractString, taxonomic_levels::Int64)
    taxonomy_table = readtable(taxonomy_file::AbstractString, separator='\t')
    split_taxonomy_df = DataFrame()
    for line in eachrow(taxonomy_table)
        newline = split_mothur_taxonomy_string(line[:Taxonomy]) |> x -> reshape(x, (1,taxonomic_levels)) |> DataFrame
        split_taxonomy_df = vcat(split_taxonomy_df,newline)
    end
    taxonomic_level_names = Dict(:x1=>:Domain,:x2=>:Phylum,:x3=>:Class,:x4=>:Order,:x5=>:Family,:x6=>:Genus)
    rename!(split_taxonomy_df,taxonomic_level_names)
    delete!(taxonomy_table, :Taxonomy)
    taxonomy_table_with_split_taxa = hcat(taxonomy_table,split_taxonomy_df)
    return taxonomy_table_with_split_taxa
end

function import_mothur_data(shared_file::AbstractString,taxonomy_file::AbstractString,metadata_file::AbstractString;taxonomic_levels::Int64=6)
    data = microbiome_data(import_mothur_shared_file(shared_file::AbstractString),import_mothur_taxonomy_file(taxonomy_file::AbstractString, taxonomic_levels),readtable(metadata_file,separator='\t',header=true))
    otu_length = string(length(data.otu[:OTU]))
    if Set(data.otu[:OTU]) == Set(data.tax[:OTU])
        println("$otu_length OTUs found in both shared and taxonomy files.")
    else
        error("OTU identifiers do not match in shared and taxonomy files. $otu_length OTUs found in shared file and "*string(length(data.tax[:OTU]))*" OTUs found in taxonomy file.")
    end
    println(string(length(names(data.otu))-1)*" samples found in shared file.")
    return data
end

function import_mothur_data(shared_file::AbstractString,taxonomy_file::AbstractString;taxonomic_levels::Int64=6)
    data = microbiome_data(import_mothur_shared_file(shared_file::AbstractString),import_mothur_taxonomy_file(taxonomy_file::AbstractString, taxonomic_levels),DataFrame(),DataFrame())
    otu_length = string(length(data.otu[:OTU]))
    if Set(data.otu[:OTU]) == Set(data.tax[:OTU])
        println("$otu_length OTUs found in both shared and taxonomy files.")
    else
        error("OTU identifiers do not match in shared and taxonomy files. $otu_length OTUs found in shared file and "*string(length(data.tax[:OTU]))*" OTUs found in taxonomy file.")
    end
    println(string(length(names(data.otu))-1)*" samples found in shared file.")
    return data
end