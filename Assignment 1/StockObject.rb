###### DEFINING STOCK CLASS ######
###### with attributes and methods ######
# for this purpose Gene class is required
require './GeneObject.rb'


###### ATTRIBUTES ######
class Stock
  # creating a variable class with @@ shared by all Gene instances 
  # this variable is an empty array to store all instances when created
  @@stocks = Array.new
  # creating attribute accesors 
  attr_accessor :seedstock
  # geneid stock attribute will be a gene object
  attr_accessor :geneid 
  attr_accessor :lastplanted 
  attr_accessor :storage
  attr_accessor :gramsremaining
 
  # initializing variables
  def initialize (params = {})
    @seedstock = params.fetch(:seedstock, "A000")
    @geneid = params.fetch(:geneid, "AT0000000")
    @lastplanted = params.fetch(:lastplanted, "00/00/0000")
    @storage = params.fetch(:storage, "cama00")
    @gramsremaining = params.fetch(:gramsremaining, "00")
    # converting gramsremaining into float in order to operate with this values
    @gramsremaining = @gramsremaining.to_f 
    # adding each object to the shared hash class variable
    @@stocks << self
  end


###### METHODS ######
###### 1. INSTANCE METHOD: planting 7 gr. of seeds for each stock record ######
# this function takes 1 argument: the number of grams of seeds wanted to plant
  def planting_seeds(grams_to_plant)
    # if having more grams remaining than wanted to plant, then plant
    if @gramsremaining > grams_to_plant
      @gramsremaining = @gramsremaining - grams_to_plant
    # if having less or equal grams remaining than wanted to plant, cannot plant
    else
      @gramsremaining = 0
      # warning message with interpolated attribute
      puts "WARNING: we have run out of Seed Stock #{@seedstock}"#code
    end
    return @gramsremaining
    
  end
###### 2. CLASS METHOD: creating new instances in Stock Class from seed_stock_data.tsv and from Gene class ######
  # this function takes 2 arguments, the list of gene objects and the file name
  def Stock.new_stocks(file_name) 
    # opening file with open method
    # skipping header line
    first_line = true
    # creating an empty hash in order to store instances
    stocks = Hash.new
    stockfile = File.open(file_name, "r")
    stockfile.each_line do |line|
      if (first_line)
        # if matches the header line
        if (line =~ /Seed_Stock\tMutant_Gene_ID\tLast_Planted\tStorage\tGrams_Remaining\n/) 
            first_line = false
            # goes to the next line
            next
        else
          return false
        end
      else
        # getting 5 arrays, 1 per record, with 5 elements each one, 1 per column
        splitted_stock_table = line.split("\t")
        # returning gene object
        gene = Gene.return_object(splitted_stock_table[1])
        # creating stock instances 
        instance_stock = Stock.new(
                :seedstock => splitted_stock_table[0],
                # assigning gene object to a stock attribute
                :geneid => gene, 
                :lastplanted => splitted_stock_table[2],
                :storage => splitted_stock_table[3],
                :gramsremaining => splitted_stock_table[4]
                )
        # storing new instance into stock hash with the stock ID as key 
        stocks[splitted_stock_table[0]] = instance_stock
      end
    end
    return stocks
  end
###### 3. CLASS METHOD: returns each object stored in the stocks array class variable ######
# necessary to input each stock object into properties of cross instances
  def Stock.return_object(stockid)
    @@stocks.each do |stock|
    return stock if stock.seedstock == stockid
  end
  end
    
end