###### CREATING OBJECTS FROM FILES (gene_information.tsv, seed_stock_data.tsv, cross_data.tsv) ######
###### PLANTING 7 GRAMS OF SEEDS IN EACH OBJECT OF STOCK CLASS AND WRITING RESULT IN NEW .tsv FILE ######
###### CHECKING LINKED GENES IN CROSS CLASS ######
# for this purpose Gene, Stock and Cross classes are required
require './GeneObject.rb'
require './StockObject.rb'
require './CrossObject.rb'

###### CREATING OBJECTS FROM FILES (gene_information.tsv, seed_stock_data.tsv, cross_data.tsv) ######
###### 1. CREATING Gene OBJECTS ######
genes = Gene.new_genes("./gene_information.tsv")
###### 2. CREATING Stock OBJECTS ######
stocks = Stock.new_stocks("./seed_stock_data.tsv")
###### 3. CREATING Cross OBJECTS ######
crosses = Cross.new_crosses("./cross_data.tsv")
puts "This is Gene Class\n"
puts "#{genes}\n"
puts "This is Stock Class, properties are Gene objects\n"
puts "#{stocks}\n"
puts "This is Cross Class, properties are Stock objects\n"
puts "#{crosses}\n"


###### PLANTING 7 GRAMS OF SEEDS IN EACH OBJECT OF STOCK CLASS AND WRITING RESULT IN NEW .tsv FILE ######
# iterating over each value of stocks hash and applying the planting.seed method to each one
File.open("new_file.tsv", "w") do |file|
  # writing the header of the new_file
  file.write("Seed_Stock\tMutant_Gene_ID\tLast_Planted\tStorage\tGrams_Remaining\n")
  stocks.values.each do |seed|
    # writing de rest of the content of the file 
    file.write (seed.seedstock + "\t")
    file.write (seed.geneid.geneid + "\t")
    file.write (Time.now.strftime("%d/%m/%Y" + "\t"))
    file.write (seed.storage + "\t")
    file.write (seed.planting_seeds(7).to_s + "\n")
  end
end




