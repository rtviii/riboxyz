# After https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4358319/
from typing import Literal

list_PolymericFactorClass = [
    "Initiation Factor",
    "Translation Factor",
    "Biogenesis Factor",
    "Release Factor",
    "Elongation Factor",
    "Termination Factor",
    "Nascent Chain",
    "Ribozyme",
    "Riboswitch",
    "RNP","mRNA","tRNA",
    "snoRNA", "lncRNA", "Other"
]

#TODO: This is a goddamn mess. tRNAs should be partioned into their own categories, likewise with factors.
    # Ribosome Assembly Factors
    #     RbfA
    #     RimM
    #     Era
    #     GTPases: EngA, ObgE

    # Initiation Factors:
    #     IF1 (Initiation Factor 1)
    #     IF2 (Initiation Factor 2)
    #     IF3 (Initiation Factor 3)

    # Elongation Factors:
    #     EF-Tu (Elongation Factor Tu)
    #     EF-G (Elongation Factor G)
    #     EF-P (Elongation Factor P)
    #     EF-Ts (Elongation Factor Ts)
    #     EF-1α (Elongation Factor 1α)
    #     EF-1β (Elongation Factor 1β)

    # Release Factors:
    #     RF1 (Release Factor 1)
    #     RF2 (Release Factor 2)
    #     RF3 (Release Factor 3)

    # Ribosome Recycling Factors:
    #     RRF (Ribosome Recycling Factor)
    #     RliB (Ribosome Large Subunit Interface B)
    #     IF6 (Initiation Factor 6)

    # RNA Helicases:
    #     RhlB (RNA Helicase B)
    #     DbpA (ATP-dependent RNA Helicase DbpA)
    #     DbpB (ATP-dependent RNA Helicase DbpB)
    #     Prp43 (Pre-mRNA Processing Factor 43)

    # GTPases:
    #     Era
    #     Obg (YjeQ)
    #     LepA

    # Non-Ribosomal Factors:
    #     HflX
    #     YjeQ
    #     YihA

    # Quality Control Factors:
    #     tmRNA (Transfer-messenger RNA)
    #     ArfA (Alternative Ribosome-Rescue Factor A)
    #     ArfB (Alternative Ribosome-Rescue Factor B)

    # Maturation Factors:
    #     RbfA (Ribosome Binding Factor A)
    #     RimM (Ribosome Maturation Factor M)

    # Trans-Translation Factors:
    #     SmpB (Small Protein B)
    #     tmRNA (Transfer-messenger RNA)

PolymericFactorClass = Literal[
    "Elongation Factor",
    "Initiation Factor",
    "Translation Factor",
    "Biogenesis Factor",
    "Release Factor",
    "Termination Factor",
    "Nascent Chain",
    "Ribozyme",
    "Riboswitch",
    "RNP",
    "mRNA",
    "tRNA",
    "snoRNA",
    "lncRNA",
    "Other"
]


NonpolymericLigandClass = Literal[
    "Antibiotic",
    "Ion",
    "Cofactor",
    "Analog",
    "Other"
]
list_NonpolymericLigandClass = [
    "Antibiotic", "Ion", "Cofactor", "Analog", "Other"]


RNAClass = Literal[
    "5SrRNA",
    "5.8SrRNA",
    "12SrRNA",
    "16SrRNA",
    "21SrRNA",
    "23SrRNA",
    "25SrRNA",
    "28SrRNA",
    "35SrRNA",
]

SSUProteinClass = Literal["bS1", "eS1", "uS2", "uS3", "uS4", "eS4", "uS5", "bS6", "eS6", "uS7", "eS7",
                          "uS8", "eS8", "uS9", "uS10", "eS10", "uS11", "uS12", "eS12", "uS13",
                          "uS14", "uS15", "bS16", "uS17", "eS17", "bS18", "uS19", "eS19", "bS20",
                          "bS21", "bTHX", "eS21", "eS24", "eS25", "eS26", "eS27", "eS28", "eS30",
                          "eS31", "RACK1"]

LSUProteinClass = Literal["uL1", "uL2", "uL3", "uL4", "uL5", "uL6", "eL6", "eL8", "bL9", "uL10",
                          "uL11", "bL12", "uL13", "eL13", "uL14", "eL14", "uL15", "eL15", "uL16",
                          "bL17", "uL18", "eL18", "bL19", "eL19", "bL20", "eL20", "bL21", "eL21",
                          "uL22", "eL22", "uL23", "uL24", "eL24", "bL25", "bL27", "eL27", "bL28",
                          "eL28", "uL29", "eL29", "uL30", "eL30", "bL31", "eL31", "bL32", "eL32",
                          "bL33", "eL33", "bL34", "eL34", "bL35", "bL36", "eL36", "eL37", "eL38",
                          "eL39", "eL40", "eL41", "eL42", "eL43", "P1P2"]

# ----------------------------------------------------------------------------------------------
# These list_ constructs are just a hack for purposes of iteration (given that typing.Literal behaves poorly in some useful contexts)
list_LSUProteinClass = ["uL1", "uL2", "uL3", "uL4", "uL5", "uL6", "eL6", "eL8", "bL9", "uL10", "uL11", "bL12", "uL13", "eL13", "uL14", "eL14", "uL15", "eL15", "uL16", "bL17", "uL18", "eL18", "bL19", "eL19", "bL20", "eL20", "bL21", "eL21", "uL22", "eL22",
                        "uL23", "uL24", "eL24", "bL25", "bL27", "eL27", "bL28", "eL28", "uL29", "eL29", "uL30", "eL30", "bL31", "eL31", "bL32", "eL32", "bL33", "eL33", "bL34", "eL34", "bL35", "bL36", "eL36", "eL37", "eL38", "eL39", "eL40", "eL41", "eL42", "eL43", "P1P2"]
list_SSUProteinClass = ["bS1", "eS1", "uS2", "uS3", "uS4", "eS4", "uS5", "bS6", "eS6", "uS7", "eS7", "uS8", "eS8", "uS9", "uS10", "eS10", "uS11", "uS12", "eS12", "uS13",
                        "uS14", "uS15", "bS16", "uS17", "eS17", "bS18", "uS19", "eS19", "bS20", "bS21", "bTHX", "eS21", "eS24", "eS25", "eS26", "eS27", "eS28", "eS30", "eS31", "RACK1"]
list_ProteinClass = [*list_LSUProteinClass, *list_SSUProteinClass]
list_RNAClass = ["5SrRNA", "5.8SrRNA", "12SrRNA", "16SrRNA", "21SrRNA",
                 "23SrRNA", "25SrRNA", "28SrRNA", "35SrRNA", "mRNA", "tRNA"]
