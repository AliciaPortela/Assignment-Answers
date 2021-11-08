require 'rest-client'
require 'json'
require './InteractionNetwork.rb'


file = "./"+ARGV[0]+""
report = "./"+ARGV[1]+""

InteractionNetwork.getting_agi(file)
InteractionNetwork.search_protein_int()
InteractionNetwork.create_networks()


##### WRITING A REPORT OF THE RESULTS OBTAINED #####
# creating a new .txt file, writing it line by line
File.open(report, 'w') do |line|
  line.puts "---------Analysis of direct and indirect interactions between a set of co-expressed genes----------"
  line.puts "\n"
  line.puts "For the #{File.open(file){|f| f.read.count("\n")}} co-expressed genes, direct and indirect interactions between them have been analysed through PSICQUIC REST API"
  line.puts "\n"
  line.puts "SCHEMATIC REPRESENTATION OF INTERACTIONS:" 
  line.puts "-Direct interactions: A-B, in which A and B are genes that belong to the co-expressed set of genes"
  line.puts "-Indirect interactions: A-B-C, in which A and C are genes that belong to the co-expressed set of genes, but not B, which is the intermediary"
  line.puts "\n"
  line.puts "CRITERIA FOR SEARCHING FOR INTERACTIONS:"
  line.puts "1. Arabidopsis thaliana species, filtered by a regex"
  line.puts "2. MIscore > 0.485, as mentioned in 10.1093/database/bau131"
  line.puts "\n"
  line.puts "RESULTS:"
  line.puts "In total, #{InteractionNetwork.get_number_networks()} networks have been found"
  line.puts "Below, it is shown one by one the networks, with the following information:" 
  line.puts "1. Number of network"
  line.puts "2. Type (direct/indirect) of interaction"
  line.puts "3. Members of the network"
  line.puts "4. KEGG annotations: KEGG ID and KEGG pathway"
  line.puts "5. GO annotations: GO ID and GO terms"
  InteractionNetwork.get_networks.each do |network| # for each interaction...
    line.puts "\n"
    line.puts "----------------------------------------------------------------------------------------------------------"
    line.puts "Network number #{network.network_num}"
    if network.members.count == 2
      line.puts "Direct interaction"
      line.puts "Members: #{network.members[0]} and #{network.members[1]}"
    elsif network.members.count == 3
      line.puts "Indirect interaction"
      line.puts "Members: #{network.members[0]} and #{network.members[2]}, with #{network.members[1]} as intermediary"
    end
    # kegg annotations
    line.puts "\n"
    line.puts "KEGG annotations:"
    # only if there is information about the KEGG in the hash
    if !network.kegg_pathway.empty?
      # if there is only one key-values pair...
      if network.kegg_pathway.count == 1
        line.puts "The #{network.kegg_pathway.keys[0]} member has KEGG ID: #{network.kegg_pathway.values[0][0]} and KEGG pathway: #{network.kegg_pathway.values[0][1]}"
      # if there is more than one key-values pair...
      else 
        network.kegg_pathway.each do |key,values|
          line.puts "The #{key} member has KEGG ID: #{values[0]} and KEGG pathway: #{values[1]}"
        end 
      end 
    else 
      line.puts "There are not any KEGG ID/pathway for any member of this network"
    end
    # go annotations
    line.puts "\n"
    line.puts "GO annotations:"
    # only if there is information about the GO in the hash
    if !network.go_terms.empty?
      if network.go_terms.count == 1
        line.puts "The #{network.go_terms.keys[0]} member has GO ID: #{network.go_terms.values[0][0]} and GO terms: #{network.go_terms.values[0][1]}"
      else 
        network.go_terms.each do |key,values|
          line.puts "The #{key} member has GO ID: #{values[0]} and GO terms: #{values[1]}"
        end
      end
    else 
      line.puts "There are not any GO ID/terms for any member of this network"
    end 
    line.puts "----------------------------------------------------------------------------------------------------------"
  end
end