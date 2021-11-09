######## this scritps does the following ########
######## 1.gets the AGIs from the ArabidopsisSubNetwork_GeneList.txt file ########
######## 2.searches for interactions between genes from that file ########
######## 3.creates networks (instances) in the class InteractionNetwork with network, KEGG and GO annotations ########
######## 4.writes a report with all the information about the performed analysis ########


# in order to make a request to a server: PSICQUIC/TOGO
require 'rest-client'
# in order to handle '.json' format
require 'json'
# required to get the definition of the InteractionNetwork class
require './InteractionNetwork.rb'

##### This script is run in the command line as follows #####
##### $ruby Assignment_2.rb ArabidopsisSubNetwork_GeneList.txt report.txt #####
# using ARGV[]
# passing the command line arguments to variables 
file = "./"+ARGV[0]+""
report = "./"+ARGV[1]+""

##### getting the AGIs #####
InteractionNetwork.getting_agi(file)
##### searching for interactions #####
InteractionNetwork.search_protein_int()
##### creating networks (instances) #####
InteractionNetwork.create_networks()


##### WRITING A REPORT OF THE RESULTS OBTAINED #####
# creating a new .txt file, writing it line by line
File.open(report, 'w') do |line|
  line.puts "\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s---------ANALYSIS OF INTERACTIONS BETWEEN A SET OF PREDICTIVELY CO-EXPRESSED GENES----------"
  line.puts "\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s---------------------PERFORMED BY: Alicia Portela Estévez----------------------"
  line.puts "\n"
  line.puts "INTRODUCTION"
  line.puts "A recent paper (DOI: 10.1371/journal.pone.0108567) executes a meta-analysis of a few thousand published co-expressed gene sets from Arabidopsis thaliana. They break these co-expression sets into ~20 sub-networks of <200 genes each, that they find consistently co-expressed with one another. One of this sub-networks has #{File.open('ArabidopsisSubNetwork_GeneList.txt'){|f| f.read.count("\n")}} predictively co-expressed genes"
  line.puts "\n"
  line.puts "OBJECTIVE"
  line.puts "The objective of this further analysis is testing if there is already information that links those #{File.open(file){|f| f.read.count("\n")}} genes into known regulatory networks, either in a direct or indirect way, and annotate members of each network with any possible KEGG Pathway or GO Terms."
  line.puts "\n"
  line.puts "MATERIAL & METHODS"
  line.puts "Materials used in this analysis are:"
  line.puts "1. The sub-network with #{File.open(file){|f| f.read.count("\n")}} predicted co-expressed genes."
  line.puts "2. BAR database from University of Toronto, in order to search for interactions, either direct or indirect, between the #{File.open(file){|f| f.read.count("\n")}} co-expressed genes, through PSICQUIC REST API."
  line.puts "3. UniProt database, in order to search for KEGG and GO annotations, through TOGO REST API."
  line.puts "\n"
  line.puts "Methods followed in this analysis are:"
  line.puts "1. Type of interactions"
  line.puts "In this analysis, 2 kind of interactions have been analysed:"
  line.puts "-Direct interactions: A-B, in which A and B are genes that belong to the co-expressed sub-network."
  line.puts "-Indirect interactions: A-B-C, in which A and C are genes that belong to the co-expressed sub-network, but not B, which is the intermediary."
  line.puts "2. Criteria for searching for interactions"
  line.puts "In this analysis, 2 criteria have been followed:"
  line.puts "-Interactors must belong to Arabidopsis thaliana species. It has been achieved by filtering by a regex, corresponding to the AGI locus code."
  line.puts "-MIscore > 0.485, as mentioned in the article with DOI: 10.1093/database/bau131."
  line.puts "\n"
  line.puts "RESULTS"
  line.puts "In total, #{InteractionNetwork.get_number_networks()} networks have been found"
  line.puts "Below, it is shown one by one the networks, with the following information:" 
  line.puts "1. Number of network"
  line.puts "2. Type of interaction: direct/indirect"
  line.puts "3. Members of the network"
  line.puts "4. KEGG annotations: KEGG ID and KEGG pathway"
  line.puts "5. GO annotations: GO ID and GO terms"
  InteractionNetwork.get_networks.each do |network| # for each interaction...
    line.puts "\n"
    line.puts "----------------------------------------------------------------------------------------------------------"
    # write type and members
    line.puts "Network number #{network.network_num}"
    # if network contains only 2 members is a direct interaction (A-B)
    if network.members.count == 2
      line.puts "Direct interaction"
      line.puts "Members: #{network.members[0]} and #{network.members[1]}"
    # if network contains 3 members is an indirect interaction (A-B-C)
    elsif network.members.count == 3
      line.puts "Indirect interaction"
      line.puts "Members: #{network.members[0]} and #{network.members[2]}, with #{network.members[1]} as intermediary"
    end
    # write kegg annotations
    line.puts "\n"
    line.puts "KEGG annotations:"
    # only if the KEGG hash is not empty
    if !network.kegg_pathway.empty?
        # iterating over the hash (key-value)
        network.kegg_pathway.each do |key,values|
          # values are always array datatype
          # if values are not empty
          if values.is_a?(Array) && !values.empty?
            # if array has only one element...
            if values.count == 1 
              line.puts "The #{key} member has KEGG ID: #{values[0][0]} and KEGG pathway: #{values[0][1]}"
            # if array has more than one element...
            elsif values.count > 1
              line.puts "The #{key} member has the following KEGG IDs/pathways:"
              # iterating over the elements 
              values.each do |value|
                line.puts "\s#{value[0]},#{value[1]}\s"
              end
            end
          # if array values is empty
          else
            line.puts "The #{key} member has not got any KEGG ID/pathway"
          end 
        end
    # if KEGG hash is empty
    else 
      line.puts "There are not any KEGG ID/pathway for any member of this network"
    end
    # write go annotations
    # same code as keggs 
    line.puts "\n"
    line.puts "GO annotations:"
    if !network.go_terms.empty?
      network.go_terms.each do |key,values|
        if values.is_a?(Array) && !values.empty?
          if values.count == 1
            line.puts "The #{key} member has GO ID: #{values[0][0]} and GO terms: #{values[0][1]}"
          elsif values.count > 1
            line.puts "The #{key} member has the following GO ID/terms:"
            values.each do |value|
              line.puts "\s#{value[0]},#{value[1]}\s"
            end
          end
        else
          line.puts "The #{key} member has not got any GO ID/term"
        end
      end
    else 
      line.puts "There are not any GO ID/terms for any member of this network"
    end 
    line.puts "----------------------------------------------------------------------------------------------------------"
  end
  line.puts "\n"
  line.puts "DISCUSSION"
  line.puts "Given a sub-network of predictively co-expressed #{File.open(file){|f| f.read.count("\n")}} number of genes by the article DOI: 10.1371/journal.pone.0108567, the goal of this further analysis has been checking if there is already available information in databases about direct or indirect interactions between those genes, and thus, creating sub-networks of interaction. Also, each network has been annotated with possible KEGG pathways and GO terms of its members."
  line.puts "This analysis has been able to find #{InteractionNetwork.get_number_networks()} networks, most of them as a result of indirect interactions between members of the sub-network list, with only one intermediary. Thus, maximum interaction depth developed in this analysis is 2. That is to say, direct interactions include 2 members of the sub-network list (A-B), whereas indirect interactions take into account 2 members of the sub-network list (A-C), which interact through only one intermediary (B)." 
  line.puts "Although this analysis has performed well, further studies may do deeper interaction analyses in the future, with more intermediaries, that will lead to create greater networks."
end