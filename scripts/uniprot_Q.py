import os
import requests

from ribctl import ASSETS, ASSETS_PATH

# Define the UniProt entry name

# File not found: /home/rtviii/dev/riboxyz/ribctl/assets/fasta_proteins_mitochondrial/bS1m.fasta
# File not found: /home/rtviii/dev/riboxyz/ribctl/assets/fasta_proteins_mitochondrial/uS2m.fasta
# File not found: /home/rtviii/dev/riboxyz/ribctl/assets/fasta_proteins_mitochondrial/mS26.fasta


uniref_clusters ={
    "bS1m" : "UniRef90_Q9Y2Q9",
    "uS2m" : "UniRef90_Q9Y399",
    "uS3m" : "UniRef90_P02381",
    "uS4m" : "UniRef90_P27929",
    "uS5m" : "UniRef90_P82675",
    "bS6m" : "UniRef90_P82932",
    "uS7m" : "UniRef90_Q9Y2R9",
    "uS8m" : "UniRef90_Q03799",
    "uS9m" : "UniRef90_P82933",
    "uS10m": "UniRef90_P82670",
    "uS11m": "UniRef90_P82912",
    "uS12m": "UniRef90_P53732",
    "uS13m": "UniRef90_P08977",
    "uS14m": "UniRef90_O60783",
    "uS15m": "UniRef90_P82914",
    "bS16m": "UniRef90_Q9Y3D3",
    "uS17m": "UniRef90_Q03246",
    "bS18m": "UniRef90_B3LS63",
    "uS19m": "UniRef90_P39697",
    "bS21m": "UniRef90_P82921",
    "mS22" : "UniRef90_P82650",
    "mS23" : "UniRef90_Q9Y3D9",
    "mS25" : "UniRef90_P82663",
    "mS26" : "UniRef90_Q9BYN8",
    "mS27" : "UniRef90_P12686",
    "mS29" : "UniRef90_P51398",
    "mS31" : "UniRef90_P82925",
    "mS33" : "UniRef90_P53305",
    "mS34" : "UniRef90_P82930",
    "mS35" : "UniRef90_P82673",
    "mS37" : "UniRef90_Q96BP2",
    "mS38" : "UniRef90_P32344",
    "mS39" : "UniRef90_Q96EY7",
    "mS40" : "UniRef90_Q9Y676",
    "mS41" : "UniRef90_P38783",
    "mS42" : "UniRef90_P47141",
    "mS43" : "UniRef90_P10662",
    "mS44" : "UniRef90_P12686",
    "mS45" : "UniRef90_P53292",
    "mS46" : "UniRef90_Q03430",
    "mS47" : "UniRef90_P28817",

    #mLSU
    "uL1m"  :"UniRef90_Q9BYD6",
    "uL2m"  :"UniRef90_P32611",
    "uL3m"  :"UniRef90_P09001",
    "uL4m"  :"UniRef90_P51998",
    "uL5m"  :"UniRef90_P36519",
    "uL6m"  :"UniRef90_P32904",
    "bL9m"  :"UniRef90_Q9BYD2",
    "uL10m" :"UniRef90_Q7Z7H8",
    "uL11m" :"UniRef90_Q9Y3B7",
    "bL12m" :"UniRef90_P52815",
    "uL13m" :"UniRef90_Q9BYD1",
    "uL14m" :"UniRef90_Q6P1L8",
    "uL15m" :"UniRef90_Q9P015",
    "uL16m" :"UniRef90_Q9NX20",
    "bL17m" :"UniRef90_P22353",
    "uL18m" :"UniRef90_Q9H0U6",
    "bL19m" :"UniRef90_P25626",
    "bL20m" :"UniRef90_Q9BYC9",
    "bL21m" :"UniRef90_Q8L9A0",
    "uL22m" :"UniRef90_Q9NWU5",
    "uL23m" :"UniRef90_Q16540",
    "uL24m" :"UniRef90_Q96A35",
    "bL27m" :"UniRef90_Q9P0M9",
    "bL28m" :"UniRef90_Q13084",
    "uL29m" :"UniRef90_Q9HD33",
    "uL30m" :"UniRef90_P20084",
    "bL31m" :"UniRef90_Q7Z7F7",
    "bL32m" :"UniRef90_Q9BYC8",
    "bL33m" :"UniRef90_O75394",
    "bL34m" :"UniRef90_Q9BQ48",
    "bL35m" :"UniRef90_Q9NZE8",
    "bL36m" :"UniRef90_Q9P0J6",
    "mL37"  :"UniRef90_Q9BZE1",
    "mL38"  :"UniRef90_Q96DV4",
    "mL39"  :"UniRef90_Q9NYK5",
    "mL40"  :"UniRef90_Q9NQ50",
    "mL41"  :"UniRef90_Q8IXM3",
    "mL42"  :"UniRef90_Q9Y6G3",
    "mL43"  :"UniRef90_Q8N983",
    "mL44"  :"UniRef90_Q9H9J2",
    "mL45"  :"UniRef90_Q9BRJ2",
    "mL46"  :"UniRef90_Q9H2W6",
    "mL48"  :"UniRef90_Q96GC5",
    "mL49"  :"UniRef90_Q13405",
    "mL50"  :"UniRef90_Q8N5N7",
    "mL51"  :"UniRef90_Q4U2R6",
    "mL52"  :"UniRef90_Q86TS9",
    "mL53"  :"UniRef90_Q96EL3",
    "mL54"  :"UniRef90_Q6P161",
    "mL57"  :"UniRef90_P36523",
    "mL58"  :"UniRef90_P22354",
    "mL59"  :"UniRef90_P23369",
    "mL60"  :"UniRef90_P14063",
    "mL61"  :"UniRef90_P32388",
    "mL62"  :"UniRef90_Q14197",
    "mL63"  :"UniRef90_Q9BQC6",
    "mL64"  :"UniRef90_Q8TAE8",
    "mL65"  :"UniRef90_Q9NP92",
    "mL66"  :"UniRef90_Q9NVS2",
    "mL67"  :"UniRef90_Q06630",
}

os.path.join(ASSETS_PATH,'')
for (mrp, uniref_cluster_id) in uniref_clusters.items():
    if mrp in [ "bS1m", "uS2m", "mS26"]:
        url               = "https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=accession%2Corganism_id%2Csequence&format=tsv&query=%28%28uniref_cluster_90%3A{}%29+NOT+%28accession%3AQ9Y2Q9%29%29&size=500".format(uniref_cluster_id)
        response          = requests.get(url)
        savepath = os.path.join(ASSETS['fasta_proteins_mitochondrial'], '{}.fasta'.format(mrp))
        if response.status_code == 200:
            with open(savepath, "w") as file:
                file.write(response.text)
        else:
            print(f"Error: Unable to retrieve data. Status code: {response.status_code}")
#+ -----------------------------------------------------------------------------------------------------------------------+
