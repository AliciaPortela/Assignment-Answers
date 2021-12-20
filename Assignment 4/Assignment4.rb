##### Assignment 4 ##### 
##### This script runs a Best Reciprocal Hits (BRH) analysis

##### BLAST 1:
##### 1) guessing BLAST program --> tblastn
##### 2) each protein in S.pombe is used as a query against all Arabidopsis genome
##### 3) S.pombe is a FlatFile object, Arabidopsis as database
##### 4) only retrieve best hit
##### BLAST 2:
##### 1) guessing BLAST program --> blastx
##### 2) each sequence in Arabidopsis (that was the best hit) is used as a query against all S.pombe proteome
##### 3) S.pombe is a database, Arabidopsis as a FlatFile
##### 4) only retrieve best hit
##### 5) if the best hit is equal to the protein of S.pombe --> orthologues candidates

##### PARAMETERS
##### For choosing parameters, the majority of papers or websites I've read, it is said that
##### E-value has not a defined cutoff. It is calculated by the formula E = m*n*2^-S, being 'm' the effective
##### length of the query and 'n' the effective length of the database
##### Source: https://calculator.academy/e-value-calculator/; https://moviecultists.com/how-to-calculate-e-value; https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
##### But for most of the websites and papers I have read, it is said that E-value > 1e-03 is sufficient for searching orthologues
##### Finally, I chose these parameters cutoff --> E-value = 1e-06 and coverage = 0.5 based on this paper:
##### Source: Gabriel Moreno-Hagelsieb, Kristen Latimer, Choosing BLAST options for better detection of orthologs as reciprocal best hits, Bioinformatics, Volume 24, Issue 3, 1 February 2008, Pages 319–324, https://doi.org/10.1093/bioinformatics/btm585



# to use BioRuby classes 
require 'bio'

######## METHODS ########
# 1. Guess BLAST program based on the type of sequencies in query and database: nucleotide/protein
# BLAST programs needed for this assignment
### blastx --> query:translatednucleotide (Arabidopsis); db:protein (S.pombe)
### tblastn --> query:protein (S.pombe); db:translatednucleotide (Arabidopsis)
# requires Bio::FlatFile objects as arguments
def blast_program(flat_fasta_query,flat_fasta_db)
  # passing the first sequences of query and db to Bio::Sequence objects
  # next_entry method of Bio::FlatFile takes the next entry of the file
  # to_biosequence to convert it into Bio::Sequence object --> apply guess to know what class belongs to
  query_first_seq = flat_fasta_query.next_entry.to_biosequence
  db_first_seq = flat_fasta_db.next_entry.to_biosequence
  # what BLAST program to use in each case?
  if query_first_seq.guess.equal? Bio::Sequence::NA and db_first_seq.guess.equal? Bio::Sequence::AA
    blast_program = 'blastx'
  elsif query_first_seq.guess.equal? Bio::Sequence::AA and db_first_seq.guess.equal? Bio::Sequence::NA
    blast_program = 'tblastn'
  else
    puts "Do not know which program to use"
  end 
  return blast_program
end

### BLAST requires:
# A sequence to be used as the query
# A BLAST-formatted database

# 2. Create a Bio::FlatFile object to: 1) guess blast program; 2) iterate over each query
def create_flat_file(fasta_file)
  flat_file = Bio::FlatFile.auto(fasta_file)
  return flat_file
end

# 3. Create a database 
def create_database(program_blast,fasta_file)
  # removing the '.fa' extension
  name_db = fasta_file.split('.')[0]
  if program_blast == 'blastx'
    `makeblastdb -in #{fasta_file} -dbtype 'prot' -out #{name_db}`
  elsif program_blast == 'tblastn'
    `makeblastdb -in #{fasta_file} -dbtype 'nucl' -out #{name_db}`
  end
end

# 4. Doing a BLAST and retrieving the best hit with parameters 
def blast_best_hit(prog,db,query)
  # defining parameters: evalue, coverage
  evalue_threshold = 1e-06 
  coverage_threshold = 0.5
  # 1. create a BLAST factory
  # local_blast_factory = Bio::Blast.local('blastn','/path/to/db')
  factory = Bio::Blast.local(prog,"#{File.dirname(db)}/#{File.basename(db,".fa")}") 
  # 2. run actual BLAST by querying the factory
  report = factory.query(query)
  if report.hits[0] != nil 
    if prog == "tblastn"
      coverage = report.hits[0].overlap.to_f/query.length.to_f
    elsif prog == "blastx"
      entry_length = query.length/3
      coverage = report.hits[0].overlap.to_f/entry_length.to_f
    end
    
    evalue = report.hits[0].evalue.to_f
    
    if coverage >= coverage_threshold && evalue <= evalue_threshold
      return report.hits[0].definition.split("|")[0].strip
    end
  end
end

# 5.Find reciprocal best hits 
def reciprocal_best_hits(file_1,file_2)
  # SETTING FILES TO GET QUERIES 
  # convert input files to Bio::FlatFile objects
  # Bio::FlatFile is a helper and wrapper class to read a biological data file
  # file 1
  flat_file_1 = create_flat_file(file_1)
  # file 2
  flat_file_2 = create_flat_file(file_2)
  
  # SETTING TYPE OF BLAST with 'blast_program' function
  # blast 1: flat_file_1 (S.pombe) as query and flat_file_2 (Arabidopsis) as db 
  blast_1 = blast_program(flat_file_1,flat_file_2)
  # blast 2: flat_file_2 (Arabidopsis) as query and flat_file_2 (S.pombe) as db 
  blast_2 = blast_program(flat_file_2,flat_file_1)
   
  # BLAST 1: query--> S.pombe; database-->Arabidopsis
  # new empty hash to store best hits in BLAST 1
  best_hits = Hash.new
  # new empty Hash to store orthologs candidates 
  orthologues_candidates = Hash.new
  # create Arabidopsis database with 'create_database' function
  create_database(blast_1,file_2)
  # rewind the file to set it in the beginning
  flat_file_1.rewind 
  # iterating over each entry in file_1 --> entry = query 
  flat_file_1.each do |query|
    # print it to the screen in order to see what's happening while running (query name)
    puts "BLAST 1 (S.pombe queries against Arabidopsis): query #{query.entry_id}..."
    # retrieving the best hit using 'blast_best_hit' function
    best_hit = blast_best_hit(blast_1,file_2,query)
    if best_hit != nil
      puts "The best hit is #{best_hit}"
      best_hits[query.entry_id] = best_hit
    elsif best_hit == nil
      puts "There are not best hits for this query"
    end
  end 
  
  # BLAST 2: query--> Arabidopsis; database-->S.pombe
  # create S.pombe database with 'create_database' function
  create_database(blast_2,file_1)
  flat_file_2.rewind 
  flat_file_2.each do |query|
    # doing BLAST only if Arabidopsis query is in the list of best hits (Hash values)
    # 'entry_id' now is the same as hits[0].definition of before 
    if best_hits.values.include? query.entry_id 
      puts  "BLAST 2 (Arabidopsis (best hits) queries against S.pombe): query #{query.entry_id}..."
      best_hit = blast_best_hit(blast_2,file_1,query)
      if best_hit != nil && best_hits[best_hit] == query.entry_id
        puts "The best hit is #{best_hit} --> #{best_hit} and #{query.entry_id} are orthologues candidates"
        # store in the dictionary pairs of reciprocal hits
        orthologues_candidates[best_hit] = query.entry_id
      end
    end
  end
  return orthologues_candidates 
end

##### RUNNING THE PROGRAM #####
##### ruby Assignment4.rb Spombe.fa Arabidopsis.fa #####
rbh = reciprocal_best_hits(ARGV[0],ARGV[1])
puts "RBH search has been completed"
# writting a report with the pairs of orthologues candidates 
File.open('report.txt', 'w+') do |line|
   line.puts "Pairs of orthologues candidates:"
   rbh.each do |key,value|
      line.puts "\t #{key} and #{value}"
   end
end
puts "Report with pairs of orthologues candidates has been generated"

