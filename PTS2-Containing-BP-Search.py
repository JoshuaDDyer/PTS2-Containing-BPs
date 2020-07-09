import gzip
from Bio import SeqIO
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import glob
import re

#input folder
inputfolder = r'E:\Bacterial Protein Sequences\ncbi-genomes-2020-06-05'
# output folder 
outputfolder = r'E:\Bacterial Protein Sequences\ncbi-genomes-2020-06-05'
#output name
outputname = '2020-07-07 PTS2-containing BPs'
# motif
motif = r"[RK][LVIQ]..[LVIHQ][LSGAK].[HQ][LAF]"
# make an empty column to append PTS2 positive results to
collateddfcolumns = ["Identifyer",
                     "Description",
                     "Sequence",
                     "2-45 N-Terminal Seq",
                     "Identified PTS2 Motif"]
collated_df = pd.DataFrame(columns = collateddfcolumns)
#retrieve list
filelist = glob.glob(inputfolder + "/*protein.faa.gz")
for i,file in enumerate(filelist):
    print ('Processing file {} of {}'.format(i+1, len(filelist)))
    # open a list for the identifyer and for the sequence
    identifyer = []
    sequence = []
    description = []
    # open faa.gz file
    with gzip.open(file, "rt") as handle:
        # now append the sequence and identifyer to the open lists
        for record in SeqIO.parse(handle, "fasta"):
            identifyer.append(record.id)
            # note that I have converted the fasta to a string here (else it
            # adds a comma between each amino acid)
            sequence.append(str(record.seq))
            description.append(record.description)
        df = pd.DataFrame(data = [identifyer,description,sequence]).T
        df.columns = ['Identifyer','Description','Sequence']
        # Now I create a new column which contains the sequence between
        # aa positions 2 - 45
        df['2-45 N-Terminal Seq'] = df['Sequence'].str[2:45]
        # now identify any sequences that contain the PTS2motif within the 
        # '2-45 N-terminal Seq' column
        PTS2mask = df['2-45 N-Terminal Seq'].str.contains(motif)
        PTS2df = df[PTS2mask]
        # for confirmations I will now loop through this df and find the 
        # motif
        search = []
        for values in PTS2df['2-45 N-Terminal Seq']:
            search.append(re.search(motif,values).group())
        PTS2df['Identified PTS2 Motif'] = search
        # now append to a collated df that will collect all PTS2 positive 
        # sequences
        collated_df = collated_df.append(PTS2df) 
print(collated_df)

collated_df.to_csv('{}/{}.csv'.format(outputfolder, outputname), index = None, header = True)











"""
# Test for slicing data ; this shows that I can slice some sequence between 
# my areas of interest, and if the sequence is shorter than the cut off point,
# then .str will simply slice as far as possible. Perfect!
testdf = pd.DataFrame(["ACTGTSDKDJFLSDACTGTSDKDJFLSDACTGTSDKDJFLSDACTGTSDKDJFLSD",
                        "ACTGTSDKDJFLSD"])
print(testdf)
testdf["new_sample"] = testdf[0].str[2:45]
print(testdf)

"""
"""
# a rough model for how to search for motifs using regex: anything
# inside square brackets is allowed for that position, '.' means that 
# ANY character - number or special inc - can be used. 
positivetestsequence = 'TTTTTTTTRVTTVLTQFTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
negativetestsequence = 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
# example of a peroxisome targeted protein peroxisomal 3-keto-acyl-CoA thiolase
# [Hevea brasiliensis] GenBank: AFJ74324.1
realtestsequence = 'MEKAINRQRVLLDHLRPSSSSSHNYESSLSASACLAGDSAAYHRTSVYGDDVVIVAAYRTPLCKSKRGGFKDTHADDLLAPVLKAVIEKTNLNPSEVGDIVVGTVLAPGSQRASECRMAAFYAGFPETVPIRTVNRQCSSGLQAVADVAAAIKAGFYDIGIGAGLESMTSNPMAWDGDVNPKVKAFEQAQNCLLPMGVTSENVAHRFGVTRQEQDQAAVESHRKAAAATASGKFKNEIIPVATKIVDPKTGHEKPVTISVDDGIRPNTSLSELGKLKPVFKKDGTTTAGNSSQVTDGAGAVLLMKRSVAMRKGLPILGVFRTFAAVGVDPAIMGIGPAVAIPAAVKAAGLELADIDLFEINEAFASQFVYCRKKLELDPEKINVNGGAMAIGHPLGATGARCLATLLHEMKRRGRDCRFGVVSMCIGTGMGAAAVFERGDAADDLCNARRVETNDLLSKDAM'
consensus = r"[RK][LVIQ]..[LVIHQ][LSGAK].[HQ][LAF]"
sequenceslist = [positivetestsequence,negativetestsequence, realtestsequence]
for seq in sequenceslist:
    if re.search(consensus, seq):
        print(seq, "perox site found")
    else:
        print (seq, "perox site not found")
"""


