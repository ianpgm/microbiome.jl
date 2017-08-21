#!/usr/bin/env julia

include("/home/mepe-db/programs/microbiome.jl/microbiome.jl")

function show_usage()
    println("usage: taxon_fraction_table.jl -s mothur_shared_file.shared -t mothur_const_taxonomy_file.cons.taxonomy -o output_filename.tsv -l taxonomic_level -h other_threshold")
    println("shared file and cons.taxonomy file are outputs from mothur")
    println("taxonomic_level is one of the following: Domain Phylum Class Order Family Genus - the default is Phylum")
    println("other_threshold is the threshold for a fraction of reads that a taxon must exceed in at least one sample to not be added to \"other\" - the default is 0.0 (0%)")
end


function parse_parameters(arguments)
    println("Taxon table script")
    println("Written by Ian Marshall for Microbial Element Cycling and Population Ecology, 2017")

    parameter_dict = Dict()

    parameter_dict["shared_file"] = ""
    parameter_dict["taxonomy_file"] = ""
    parameter_dict["output_filename"] = "taxon_fraction_table_output.tsv"
    parameter_dict["taxonomic_level"] = "Phylum"
    parameter_dict["other_threshold"] = 0.0



    if length(arguments) == 0
        println("No arguments given")
        show_usage()
        exit()
    end

    if first(arguments) == "-h"||first(arguments) == "--help"
        show_usage()
        exit()
    end

    for (marker,value) in zip(arguments[1:2:length(arguments)],arguments[2:2:length(arguments)])
        println(marker*" = "*value)
        if marker == "-s"||marker == "--shared_file"
            parameter_dict["shared_file"] = value
        elseif marker == "-t"||marker == "--taxonomy_file"
            parameter_dict["taxonomy_file"] = value
        elseif marker == "-l"||marker == "--taxonomic_level"
            parameter_dict["taxonomic_level"] = value
        elseif marker == "-o"||marker == "--output_filename"
            parameter_dict["output_filename"] = value
        elseif marker == "-h"||marker == "--other_threshold"
            parameter_dict["other_threshold"] = parse(Float64,value)
        end

    end
    
    if parameter_dict["shared_file"] == ""
        println("Error: no shared file given")
        show_usage()
        exit()
    elseif parameter_dict["taxonomy_file"] == ""
        println("Error: no taxonomy file given")
        show_usage()
        exit()
    end

    return parameter_dict
end

parameter_dict = parse_parameters(ARGS)
data = import_mothur_data(parameter_dict["shared_file"],parameter_dict["taxonomy_file"])
collapsed_data = collapse_taxon(data, Symbol(parameter_dict["taxonomic_level"]), other_threshold=parameter_dict["other_threshold"]) |> x -> fractionate_table(x,Symbol(parameter_dict["taxonomic_level"]))
writetable(parameter_dict["output_filename"],collapsed_data,separator='\t')