# a script to generate a barcode whitelist for Parse Bio
## this may take a couple minutes ince there are a lot of barcode combinations to make

from sys import argv
import pandas as pd

infolder = argv[1]
protocol = argv[2].lower() # mini, normal, or mega
outsuf = argv[3]

# get the right bc1 table
match protocol:
  case "mini":
    bc1tab = pd.read_csv(infolder+"/bc_data_n24_v4.csv", header=0)
    
  case "normal":
    bc1tab = pd.read_csv(infolder+"/bc_data_n96_v4.csv", header=0)
    
  case "mega":
    bc1tab = pd.read_csv(infolder+"/bc_data_n192_v4.csv", header=0)

# get bc2 and bc3 tables
## yes they are the same table
bc2tab = pd.read_csv(infolder+"/bc_data_v1.csv", header=0)
bc3tab = pd.read_csv(infolder+"/bc_data_v1.csv", header=0)

# bc3-bc2-bc1 (not how the variables are named!)
bc1bc2bc3 = []
w1w2w3 = []
bci1 = []
bci2 = []
bci3 = []
for bc1i in range(bc1tab.shape[0]):
  for bc2i in range(bc2tab.shape[0]):
    for bc3i in range(bc3tab.shape[0]):
      bc1bc2bc3.append(bc3tab.sequence[bc3i]+bc2tab.sequence[bc2i]+bc1tab.sequence[bc1i])
      w1w2w3.append(bc3tab.well[bc3i]+";"+bc2tab.well[bc2i]+";"+bc1tab.well[bc1i])
      bci1.append(bc1tab.bci[bc1i])
      bci2.append(bc2tab.bci[bc2i])
      bci3.append(bc3tab.bci[bc3i])
      
whitelist = pd.DataFrame(bc1bc2bc3)
wellref = pd.DataFrame({"barcode": bc1bc2bc3, "well_combo": w1w2w3, 
                        "bci1": bci1, "bci2": bci2, "bci3": bci3})

whitelist.to_csv("ParseBio_whitelist_"+protocol+"_"+outsuf+".csv", index=False, header=False)
wellref.to_csv("ParseBio_wells_"+protocol+"_"+outsuf+".csv", index=False)
