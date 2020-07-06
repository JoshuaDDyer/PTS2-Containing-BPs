import re

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