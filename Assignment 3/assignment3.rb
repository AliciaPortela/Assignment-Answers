##### This script has been created by Alicia Portela Estévez #####
##### It does many things: #####
##### 1. searches for CTTCTT repeats in exons of ArabidopsisSubNetwork_GeneList.txt gene list and adds repeats with positions as gene features #####
##### 2. writes a gff3 file with repeat features (AGI locus code as seqname and positions inside the gene) #####
##### 3 writes a gff3 file with repeat features (Chr number as seqname and positions inside the chromosome) #####
##### 4 writes a report with those genes (AGI locus code) that do NOT have exons with the CTTCTT repeat #####

# upload all BioRuby classes
require 'bio'
# in order to make a request to a server
require 'rest-client'

### DEFINING METHODS ###
# 1. contact the web and handle possible errors
def fetch(url, headers = {accept: "*/*"}, user = "", pass="") 
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

# 2. searches for 'CTTCTT' repeat in exons from each ArabidopsisSubNetwork_GeneList.txt genes and
# adds repeats as Bio::Feature object of each Bio::Sequence object
def add_repeat_features(file)
  # new hash to store each Bio::Sequence object
  @all_biosequence = Hash.new
  # for each gene in file...
  File.open(file).each do |gene|
    # AGI locus code without '\n'
    gene.strip! 
    # search for gene info in ensemblgenomesgene database 
    response = fetch("http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=ensemblgenomesgene&format=embl&id=#{gene}")
    if response
      record = response.body
      # store database record as a Bio::EMBL object
      entry = Bio::EMBL.new(record)
      # creating a new array to store positions of 'CTTCTT' repeats
      positions = Array.new 
      # convert EMBL record to Bio::Sequence object
      bioseq = entry.to_biosequence
      # features method of Bio::Sequence Class creates Bio::Feature array objects
      # iterating over array objects
      bioseq.features.each do |feature|
        # feature method of Bio::Feature Class is the type of feature (i.e: exon, gene...)
        # next iteration unless feature = exon
        next unless feature.feature == 'exon'
        # locations method of Bio::Feature Class creates Bio::Locations object
        locations = feature.locations
        # locations method of Bio::Locations Class creates an array of Bio::Location objects
        location = locations.locations
        # iterating over each Bio::Location object of the array
        location.each do |loc|
          # subtract exon sequence from original sequence with 'subseq' method --> always forward strand!
          # with positions from Bio::Location object obtained with 'from' and 'to' mehods 
          # these positions correspond always to forward strand 
          exon = bioseq.subseq(loc.from, loc.to)
          # next iteration if there is not exon sequence (nil)
          next if exon == nil 
          # if exon is in the forward strand...
            if loc.strand == 1
              # search for 'CTTCTT' repeat 
              # defining sequence repeat
              rep_for =(Bio::Sequence::NA.new("CTTCTT")).to_re
              # if exon contains the 'CTTCTT' repeat...
              if exon.match(rep_for)
                # get start position of match 'CTTCTT' in exon
                position_match_ex = exon.index(rep_for)
                # positions (start and end) of match 'CTTCTT' in whole sequence 
                position_match_seq = [position_match_ex + loc.from, position_match_ex + loc.from + 5]
                # adding repeat 'CTTCTT' feature (with create_features method defined below) unless same 'CTTCTT' positions have just been added --> avoid repeating the same 'CTTCTT' feature 
                bioseq.features << create_features("#{position_match_seq[0]}..#{position_match_seq[1]}",loc.strand) unless positions.include?(position_match_seq)
                # adding the Bio::Sequence object to the hash with AGI locus code as key
                @all_biosequence[gene] = bioseq
                # adding positions to array
                positions << position_match_seq
              else 
                @all_biosequence[gene] = bioseq 
              end
            # if exon is in the forward strand...
            elsif loc.strand == -1
              # search for 'AAGAAG' repeat (reverse complement of 'CTTCTT')
              # defining sequence repeat
              rep_rev =(Bio::Sequence::NA.new("AAGAAG")).to_re
               # if exon contains the 'AAGAAG' repeat...
              if exon.match(rep_rev)
                # start position of match 'AAGAAG' in exon
                position_match_ex = exon.index(rep_rev)
                # positions (start and end) of match 'AAGAAG' in whole sequence 
                position_match_seq = [position_match_ex + loc.from, position_match_ex + loc.from + 5]
                bioseq.features << create_features("complement(#{position_match_seq[0]}..#{position_match_seq[1]})",loc.strand) unless positions.include?(position_match_seq)
                @all_biosequence[gene] = bioseq
                positions << position_match_seq
              else
                @all_biosequence[gene] = bioseq
              end
            end  
        end
      end
    end
  end
  return @all_biosequence 
end

# 3. creating a new Bio::Feature object (used in method above)
def create_features(position,strand) 
  ft=Bio::Feature.new('myrepeat',position) 
  ft.append(Bio::Feature::Qualifier.new('repeat_region','CTTCTT'))
  ft.append(Bio::Feature::Qualifier.new('note', 'found by repeatfinder 2.0'))
  if strand == 1
    ft.append(Bio::Feature::Qualifier.new('strand', '+'))
  elsif strand == -1
    ft.append(Bio::Feature::Qualifier.new('strand', '-'))
  end
end

# 4. creating a gff3 file for all repeat regions: genes as seqname and positions inside the gene
def gff3_file_genes(file)
  # open a file and write it line by line
  File.open(file, 'w') do |line|
  # The first line of a GFF3 file must be a comment that identifies the version
  line.puts "##gff-version 3\n"
  # iterating over the hash for each key (AGI locus code) and value (Bio::Sequence object)
  @all_biosequence.each do |gene, biosequence|
    # iterating over each Bio::Feature object
    biosequence.features.each do |feature|
      # next iteration if the name of the feature is not "myrepeat"
      next unless feature.feature == "myrepeat"
      ## defining variables (each gff3 column)
      # source will be always BioRuby
      source = "BioRuby"
      # Hash constructed from qualifier objects
      qualifiers = feature.assoc
      # type will be the key of CTTCTT value
      type = qualifiers.key('CTTCTT')
      # don't care about the score
      score = '.'
      # strand will be the value of 'strand' key
      strand = qualifiers['strand']
      # don't care about the phase
      phase = '.'
      # don't care about the attributes
      attributes = '.'
      # creates Bio::Locations object
      locations = feature.locations
      # creates an array of Bio::Location objects
      location = locations.locations
      # iterating over the array
      location.each do |loc|
        # start position 
        start = loc.from
        # end position
        finish = loc.to
        # writting the 9 col., tab-separated format for each repeat with predefined variables 
        line.puts "#{gene}\t#{source}\t#{type}\t#{start}\t#{finish}\t#{score}\t#{strand}\t#{phase}\t#{attributes}\n"
      end 
    end
  end
end
end

# 5. creating a gff3 file for all repeat regions: chromosomes as seqname and positions inside the chromosome
def gff3_file_chr(file) 
  File.open(file, 'w') do |line|
  # The first line of a GFF3 file must be a comment that identifies the version
  line.puts "##gff-version 3\n"
  @all_biosequence.each do |gene, biosequence|
    # primary_accession method gives the primary accession number (i.e. chromosome:TAIR10:3:5506931:5508414:1)
    # separating them by ':' and select 3 column --> start position of gene inside the chromosome 
    start_gene = biosequence.primary_accession.split(":")[3]
    # number of chromosome 
    chr = biosequence.entry_id
    # same as before 
    biosequence.features.each do |feature|
      # next iteration if the name of the feature is not "myrepeat"
      next unless feature.feature == "myrepeat"
      source = "BioRuby"
      qualifiers = feature.assoc
      type = qualifiers.key('CTTCTT')
      score = '.'
      strand = qualifiers['strand']
      phase = '.'
      attributes = '.'
      locations = feature.locations
      location = locations.locations
      location.each do |loc|
        # start position of CTTCTT inside the gene 
        start_in_gene = loc.from.to_i
        # start position of CTTCTT inside the chromosome
        start_in_chr = start_gene.to_i + start_in_gene
        # end position of CTTCTT inside the gene
        finish_in_gene = loc.to.to_i
        # start position of CTTCTT inside the chromosome
        finish_in_chr = start_gene.to_i + finish_in_gene
        line.puts "Chr#{chr}\t#{source}\t#{type}\t#{start_in_chr}\t#{finish_in_chr}\t#{score}\t#{strand}\t#{phase}\t#{attributes}\n"
      end 
    end
  end
end
end

# 6. writing a report file with genes of the ArabidopsisSubNetwork_GeneList.txt list that do NOT have exons with the CTTCTT repeat
def report(file)
  # open a file and write it line by line
  File.open(file, 'w') do |line|
  line.puts "Genes (represented by AGI locus code) within ArabidopsisSubNetwork_GeneList.txt that do NOT have exons with the 'CTTCTT' repeat:"
  # iterating over hash
  @all_biosequence.each do |gene, biosequence|
    # creating a new array to store features 
    features = Array.new
    # iterating over the features 
    biosequence.features.each do |feature|
      # saving feature in array
      features << feature.feature 
    end
    # next iteration if the array contains "myrepeat" feature 
    next if features.include?("myrepeat")
    # writing the AGI locus code of those genes that don't have "myrepeat" as feature
    line.puts "#{gene}"
  end
  end

end

##### RUNNING THE PROGRAM #####
##### ruby assignment3.rb ArabidopsisSubNetwork_GeneList.txt results_gene.gff3 results_chr.gff3 report.txt #####
add_repeat_features(ARGV[0])
gff3_file_genes(ARGV[1])
gff3_file_chr(ARGV[2])
report(ARGV[3])
puts "'CTTCTT' repeat has been searched over each exon of #{ARGV[0]} list genes"
puts "Found repeats have been added to Bio::Sequence as Bio::Feature object"
puts "GFF3 (#{ARGV[1]}) file with 'CTTCTT' features has been generated --> with AGI locus code as seqname and positions within the gene"
puts "GFF3 (#{ARGV[2]}) file with 'CTTCTT' features has been generated --> with Chr as seqname and positions within entire chromosome"
puts "Plain text (#{ARGV[3]}) file has been created with those genes that do not have CTTCTT repeats in their exons"
