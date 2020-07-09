import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

outputpath = r'C:\Users\jdeed\OneDrive - University of Exeter\Peroxisome-Targeted BPs\Bacterial sequence files'
outputname = '2020-06-05 Bacterial_proteins_positive_SKL (Processed)'

df = pd.read_csv(r'C:/Users/jdeed/OneDrive - University of Exeter/Peroxisome-Targeted BPs/Bacterial sequence files/2020-06-05 Bacterial_proteins_positive_SKL (Original).csv')
description = df["Description"]
# To remove the gene name in the description column I used i.split with the
# maxsplit parameter set to '1', and removed the first character
# as designated by the square brackets, 
# then iterated through each row, assigning them to a new list then
# replacing the original Description column with the 'newdescription' list
newdescription = []
for i in description:
    i = i.split(' ',1)[1]
    newdescription.append(i)
df['Description'] = newdescription
# extract bacterial names into new column - I used this solution https://stackoverflow.com/questions/4894069/regular-expression-to-return-text-between-parenthesis/4894156#4894156
# to solve this - essentially finding whatever lies between the square brackets
# assigning that to a new series and adding that new series (descriptionseries)
# to the end of the df
descriptionseries = df['Description']
descriptionseries = descriptionseries.apply(lambda st: st[st.find('[')+1:st.find("]")])
df['Species'] = descriptionseries

# remove bacterial names from description to make for better reading
# https://stackoverflow.com/questions/20894525/how-to-remove-parentheses-and-all-data-within-using-pandas-python
descriptionseries1 = df['Description']
descriptionseries1 = descriptionseries1.str.replace(r"\[.*\]","")
df['Description'] = descriptionseries1
#Remove the "MULTISPECIES: " mark
descriptionseries2 = df['Description']
descriptionseries2 = descriptionseries2.str.replace("MULTISPECIES: ","")
# Remove all rows with duplicate values - I want to do this based on 
# the gene nae, nothing else
df.sort_values("Identifyer",inplace = True)
df.drop_duplicates(subset = "Identifyer", keep = False, inplace = True)
print(df)

df.to_csv('{}/{}.csv'.format(outputpath,outputname), index = None, header = True,)