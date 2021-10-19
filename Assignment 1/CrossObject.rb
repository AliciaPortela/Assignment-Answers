###### DEFINING CROSS CLASS ######
###### with attributes and methods ######
# for this purpose Gene and Stock classes are required
require './GeneObject.rb'
require './StockObject.rb'


###### ATTRIBUTES ######
class Cross
  # creating a variable class with @@ shared by all Gene instances 
  # this variable is an empty array to store all instances when created
  @@crosses = Array.new
  # creating attribute accesors
  # parent 1 and parent 2 cross attributes are stock objects 
  attr_accessor :parent1 
  attr_accessor :parent2
  attr_accessor :F2wild 
  attr_accessor :F2P1
  attr_accessor :F2P2 
  attr_accessor :F2P1P2
  
  # initializing variables
  def initialize (params = {})
    # adding each object to the shared hash class variable
    @@crosses << self
    @parent1 = params.fetch(:parent1, "0000")
    @parent2 = params.fetch(:parent2, "0000")
    @F2wild = params.fetch(:F2wild, "000")
    @F2P1 = params.fetch(:F2P1, "000")
    @F2P2 = params.fetch(:F2P2, "000")
    @F2P1P2 = params.fetch(:F2P1P2, "000")
  end


###### METHODS ######
###### 1. INSTANCE METHOD: for checking if there is linckage between genes ######
  def checking_linckage()
    # must apply a chi-square test to check the linckage between genes
    # calculating the total real number of descendents 
    total_sum = @F2wild.to_f + @F2P1.to_f + @F2P2.to_f + @F2P1P2.to_f
    
    # calculating expected observations
    expected_F2wild = 9.0/16.0*total_sum.to_f
    expected_F2P1 = 3.0/16.0*total_sum.to_f
    expected_F2P2 = 3.0/16.0*total_sum.to_f
    expectedF2P1P2 = 1.0/16.0*total_sum.to_f
    
    # chi-square value (following chi-square formula)
    chi_square = ((@F2wild.to_f - expected_F2wild)**2/expected_F2wild)+((@F2P1.to_f - expected_F2P1)**2/expected_F2P1)+((@F2P2.to_f - expected_F2P2)**2/expected_F2P2)+((@F2P1P2.to_f - expectedF2P1P2)**2/expectedF2P1P2)
    # chi-square value por p < 0.05
    if chi_square > 7.815
        puts "Recording: #{@parent1.geneid.genename} is genetically linked to #{@parent2.geneid.genename} with chisquare score #{chi_square}"
    end
  end
    
###### 2. CLASS METHOD: creating new instances in Cross Class from cross_data.tsv and from Stock class ######
  # this function takes 2 arguments, the list of stock objects and the file name
  def Cross.new_crosses(file_name)
    # opening file with open method
    # skipping header line
    first_line = true
    # creating an empty array in order to store instances
    crosses = Hash.new 
    crossfile = File.open(file_name, "r")
    crossfile.each_line do |line|
      if (first_line)
        # if matches the header line
        if (line =~ /Parent1\tParent2\tF2_Wild\tF2_P1\tF2_P2\tF2_P1P2\n/) 
            first_line = false
            # goes to the next line
            next
        else
          return false
        end
      else
        # getting 5 arrays, 1 per record, with 6 elements each one, 1 per column
        splitted_cross_table = line.split("\t")
        stock1 = Stock.return_object(splitted_cross_table[0])
        stock2 = Stock.return_object(splitted_cross_table[1])
        # creating new instances for the Cross class
        instance_cross = Cross.new(
                # assigning the object from stocks list as a parent 1 and parent 2 properties in cross objects
                # with the value of parent 1 and parent 2 in cross table
                :parent1 => stock1,
                :parent2 => stock2,
                :F2wild => splitted_cross_table[2], 
                :F2P1 => splitted_cross_table[3],
                :F2P2 => splitted_cross_table[4],
                :F2P1P2 => splitted_cross_table[5]
                )
        # storing new instance into cross array with the parent1 value as name 
        crosses[splitted_cross_table[0]] = instance_cross
      end
    end
    return crosses
  end
  
end