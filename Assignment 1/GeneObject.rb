###### DEFINING GENE CLASS ######
###### defining attributes and methods ######


###### ATTRIBUTES ######
class Gene
  # creating a variable class with @@ shared by all Gene instances 
  # this variable is an empty array to store all instances when created
  @@genes = Array.new
  # creating attribute accesors 
  attr_accessor :geneid
  attr_accessor :genename 
  attr_accessor :phenotype
  # creating a new attribute for later 
  attr_accessor :linkedto

  
  # initializing variables 
  def initialize (params = {})
    # adding each object to the shared hash class variable
    @@genes << self
    @geneid = params.fetch(:geneid, "000")
    @genename = params.fetch(:genename, 'name')
    @phenotype = params.fetch(:phenotype, 'Some Phenotype')
    # linked will be nil because by now is an empty value
    @linked = nil
  
  end

  
###### CLASS METHODS ######
###### 1. creates new instances in Gene Class from gene_information.tsv ######
  def Gene.new_genes(file_name) # this function takes one argument, the file name 
    first_line = true
    # creating an empty hash to hold Gene instances 
    genes = Hash.new
    # checking valid path with file? method of File object
    unless File.file?(file_name)
      return false
    end
    # opening the file and reading it line by line
    genefile = File.open(file_name, "r")
    genefile.each_line do |line|
      if (first_line)
        # if matches the header line
        if (line =~ /Gene_ID\tGene_name\tmutant_phenotype\n/) 
          # skipping the header and going to the next line
          first_line = false
          next
        else
          return false
        end
      else
        # getting 5 arrays, 1 per record, with 3 elements each one, 1 per column
        splitted_gene_table = line.split("\t")
        # creating a new instance for the Gene class
        instance_gene = Gene.new(
            :geneid => splitted_gene_table[0],
            :genename => splitted_gene_table[1],
            :phenotype => splitted_gene_table[2]
            )
        # storing new instance into genes hash with the gene ID being the key 
        genes[splitted_gene_table[0]] = instance_gene
      end
    end
    return genes
  end
  
###### 2. returns each object stored in the genes array class variable ######
# necessary to input each gene object into properties of stock instances 
  def Gene.return_object(geneid)
    @@genes.each do |gene|
    return gene if gene.geneid == geneid
  end
  end
  
  
end