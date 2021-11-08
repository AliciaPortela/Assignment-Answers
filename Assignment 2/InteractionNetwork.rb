###### DEFINING InteractionNetwork class ######
# this class contains all the interaction networks of'ArabidopsisSubNetwork_GeneList.txt' genes, as instances
# each interaction network (instance) will have members, KEGG Pathways and GO Terms

class InteractionNetwork
  ## defining and initializing class variables ##
  # 1st. counts number of instances in the class. Begins at 0
  @@number_networks = 0
  # 2nd. an empty array to store all the instances of the class
  @@all_networks = Array.new
  # 3rd. an empty array to store the AGIs from the 'ArabidopsisSubNetwork_GeneList.txt' file
  @@agis = Array.new
  
  ## defining attribute accesors (instance variables) ##
  # 1. assigning a number to each network
  attr_accessor :network_num
  # 2. members of the interacting network
  # it will be an array containing all the AGI codes of the members of the interaction network
  attr_accessor :members
  # 3. KEGG pathways
  attr_accessor :kegg_pathway
  # 4. GO Terms
  attr_accessor :go_terms
  
  ## initializing instance variables ##
  def initialize(params = {})
    # each time an instance is created with .new() method...
    # 1. assigning a number to the network
    @network_num = params.fetch(:network_num, 0)
    # 2. an array containing all the AGIs of members of the interaction network is created
    @members = params.fetch(:members, 'NA')
    # 3. anotation of KEGG ID and pathway name of all memebers in the interaction network
    @kegg_pathway = annotate_kegg(members = @members)
    # 4. anotation of GO ID and GO terms of all members of the interaction network
    @go_terms = annotate_go(members = @members)
    # 5. the total number of instances increases by 1
    @@number_networks +=1
    # 3. all instances stored in the networks array
    @@all_networks << self 
  
  end

  ## DEFINING CLASS METHODS WITH self ##
  ## 1. creating an array with all AGIs in lowercase codes from the file ('ArabidopsisSubNetwork_GeneList.txt')
  # lowercase is needed to compare with the AGIs from the web
  def self.getting_agi(file)
    File.open(file).each do |gene|
      @@agis << gene.strip.downcase
    end
    return @@agis
  end
  
  ## 2. attempts to contact the website and handles possible errors ##
  def self.fetch(url, headers = {accept: "*/*"}, user = "", pass="") 
  response = RestClient::Request.execute({
    method: :get,
    url: url.to_s,
    user: user,
    password: pass,
    headers: headers})
  return response
  # handling-errors
  rescue RestClient::ExceptionWithResponse => e # si pasa este error rescatalo de este codigo
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue RestClient::Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  rescue Exception => e
    $stderr.puts e.inspect
    response = false
    return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  end
  
  ## 3. searches for protein-protein interaction in the Bar database ##
  def self.search_protein_int() 
    # creating an empty class variable hash to store the interactions
    @@interactions = Hash.new
    # making a list equal to the @@agis class variable, where AGIs are stored
    agis = @@agis
      # for each AGI code inside that list...
      agis.each do |agi|
      # search interactions in BAR database 
      res = InteractionNetwork.fetch("http://bar.utoronto.ca:9090/psicquic/webservices/current/search/query/#{agi}")
      # res is either the response object (RestClient::Response), or false, test it with 'if'
      if res
        # res.body returns data, in which each interaction is separated by new line
        # creating an array in which each element is an interaction
        body = res.body.split("\n") 
        # for each interaction for this AGI...
        body.each do |int| 
          # ...give me an array in which each element is a column 
          array = int.split("\t") 
          # indices 2 and 3 array (columns 3 and 4 of data) --> interactors
          ##### filtering by only Arabidopsis thaliana with regex #####
          # go to the next iteration if nil values
          next if array[2].match(/(A[Tt]\w[Gg]\d\d\d\d\d)/) == nil
          next if array[3].match(/(A[Tt]\w[Gg]\d\d\d\d\d)/) == nil 
          # saving data in variables and converting to string
          int_1 = array[2].match(/(A[Tt]\w[Gg]\d\d\d\d\d)/)[0].to_s
          int_2 = array[3].match(/(A[Tt]\w[Gg]\d\d\d\d\d)/)[0].to_s
          # index 14 array (column 15 of data) --> MIscore
          # saving matched score with regex in variable and converting to float 
          miscore = array[14].match(/(0\.\d+)/)[0].to_f
          ##### filtering by MIscore #####
          # cutoff of MIscore is 0.485 (source: 10.1093/database/bau131)
          # go to the next iteration if MIscore < 0.485
          next if miscore < 0.485
          
          
          ##### searching fot DIRECT interactions (between genes in the list) #####
          # 1. if interactor 1 is the AGI used for search and interactor 2 is not the AGI used for search and interactor 2 is in the AGIs list...
          if int_1.downcase == agi && int_2.downcase != agi && agis.include?(int_2.downcase)
            # ...save interactor 2 in this hash with the AGI used for search as key
            @@interactions[agi.upcase] = int_2.upcase
          # 2. if interactor 2 is the AGI used for search and interactor 1 is not the AGI used for search and interactor 1 is in the AGIs list... 
          elsif int_1.downcase != agi && int_2.downcase == agi && agis.include?(int_1.downcase)
            # ...save interactor 1 in this hash with the AGI used for search as key
            @@interactions[agi.upcase] = int_1.upcase
            
          ##### searching for INDIRECT interactions (between genes in the list with an intermediary)
          # 3. if interactor 1 is the AGI used for search and interactor 2 is not the AGI used for search and interactor 2 is not in the AGIs list...
          elsif int_1.downcase == agi && int_2.downcase != agi && !agis.include?(int_2.downcase)
            # search interactors for int_2 that are in the list
            # same as before...
            res = InteractionNetwork.fetch("http://bar.utoronto.ca:9090/psicquic/webservices/current/search/query/#{int_2}")
            if res 
              body2 = res.body.split("\n") 
              body2.each do |int2| 
                array = int2.split("\t") 
                next if array[2].match(/(A[Tt]\w[Gg]\d\d\d\d\d)/) == nil
                next if array[3].match(/(A[Tt]\w[Gg]\d\d\d\d\d)/) == nil
                int_12 = array[2].match(/(A[Tt]\w[Gg]\d\d\d\d\d)/)[0].to_s
                int_22 = array[3].match(/(A[Tt]\w[Gg]\d\d\d\d\d)/)[0].to_s
                if int_12.downcase == int_2.downcase && int_22.downcase != int_2.downcase && agis.include?(int_22.downcase) && int_22.downcase != agi
                  # save the intermediary AGI that is not in the list and indirect interactor that is in the list with AGI used for search as key
                  @@interactions[agi.upcase] = [int_2.upcase,int_22.upcase]
                elsif int_12.downcase != int_2.downcase && int_22.downcase == int_2.downcase && agis.include?(int_12) && int_22.downcase != agi
                  # save the intermediary AGI that is not in the list and other interactor that is in the list with AGI used for search as key
                  @@interactions[agi.upcase] = [int_2.upcase,int_12.upcase]
                end
              end
            else
              puts "the Web call failed - see STDERR for details..."
            end  
          # 3. if interactor 2 is the AGI used for search and interactor 1 is not the AGI used for search and interactor 1 is not in the AGIs list...
          elsif int_1.downcase != agi && int_2.downcase == agi && !agis.include?(int_1.downcase)
            # search interactors for int_1 that are in the list
            # same as before...
            res = InteractionNetwork.fetch("http://bar.utoronto.ca:9090/psicquic/webservices/current/search/query/#{int_1}")
            if res 
              body3 = res.body.split("\n") 
              body3.each do |int3| 
                array = int3.split("\t") 
                next if array[2].match(/(A[Tt]\w[Gg]\d\d\d\d\d)/) == nil
                next if array[3].match(/(A[Tt]\w[Gg]\d\d\d\d\d)/) == nil
                int_123 = array[2].match(/(A[Tt]\w[Gg]\d\d\d\d\d)/)[0].to_s
                int_223 = array[3].match(/(A[Tt]\w[Gg]\d\d\d\d\d)/)[0].to_s
                if int_123.downcase == int_1.downcase && int_223.downcase != int_1.downcase && agis.include?(int_223.downcase) && int_223.downcase != agi
                  @@interactions[agi.upcase] = [int_1.upcase,int_223.upcase]
                elsif int_123.downcase != int_1.downcase && int_223.downcase == int_1.downcase && agis.include?(int_123.downcase) && int_123.downcase != agi
                  @@interactions[agi.upcase] = [int_1.upcase,int_123.upcase]
                end
              end 
            else
              puts "the Web call failed - see STDERR for details..."
            end 
          end  
        end
      else 
        puts "the Web call failed - see STDERR for details..."
      end
      end
    return @@interactions
  end
  
  ## 4. creates instances (networks) in the InteractionNetwork class ##
  def self.create_networks()
    # setting the number of the network at the beginning
    number = 0
    # for each entry on the hash where are all interactions...
    @@interactions.each do |key,value|
      # create an empty array to store members of interaction
      members = Array.new
      # increase by 1 the number of the network
      number += 1
      # store in members array members of indirect interactions
      if value.is_a?(Array) && key.is_a?(String)
        value.each do |val|
          members << val
        end
        members << key
        InteractionNetwork.new(:network_num => number, :members => members)
      # store in members array members of direct interactions
      elsif value.is_a?(String) && key.is_a?(String)
        members = [key,value]
        InteractionNetwork.new(:network_num => number, :members => members)
      end
    end
  end
  
  ## 5. return the class variable array containing all instances (networks) created in InteractionNetwork class ##
  def self.get_networks()
    return @@all_networks
  end
  
  ## 6. return the number of networks that InteractionNetwork class has ##
  def self.get_number_networks()
    return @@number_networks
  end
  
  
  ## DEFINING INSTANCE METHODS ##
  ## 1. annotate KEGG ID and KEGG pathways for members of a network ##
  def annotate_kegg(members)
    # creating an empty array to store KEGG annotations
    kegg_annotations = Hash.new
    # for each member of the network...
    members.each do |member|
      res = InteractionNetwork.fetch("http://togows.org/entry/genes/ath:#{member}/pathways.json")
      if res 
        data = JSON.parse(res.body)
        # data[0] is a hash 
        # iterating for each entry of this hash
        # key and value
        data[0].each do |kegg_id,kegg_p|
          # save key --> KEGG ID and value --> KEGG pathway in the array
          kegg_annotations[member] = [kegg_id,kegg_p]
        end 
      else 
        puts "the Web call failed - see STDERR for details...kegg"
      end
    end
    return kegg_annotations
  end
  ## 2. annotate GO ID and GO terms for members of a network ##
  def annotate_go(members)
    # creating an empty array to store GO annotations
    go_annotations = Hash.new
    # for each member of the network...
    members.each do |member|
      res = InteractionNetwork.fetch("http://togows.org/entry/ebi-uniprot/#{member}/dr.json")
      if res 
        data = JSON.parse(res.body)
        data[0]["GO"].each do |go|
          # go to the next iteration if the first entry does not match a GO ID
          next unless go[0].match(/GO:\d+/)
          # go to the next iteration if the second entry does not match a p = biological process, followed by string
          next unless go[1].match(/P:\w+/)
          go_id = go[0]
          go_term = go[1].match(/:(.+)/)[1]
          go_annotations[member] = [go_id,go_term]
        end
      else 
        puts "the Web call failed - see STDERR for details...go"
      end
    end
    return go_annotations
  end      
end