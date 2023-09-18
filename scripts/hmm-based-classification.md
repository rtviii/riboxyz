## Classification of polypeptides using HMMs

### Gist

For each established polypeptide nomenclature class (i.e. uL2/23SrRNA/eF-Tu/eIF-5/etc) an HMM exists. When a new sequence is received (at ETL time), it is compared against every HMM. The HMM with the highest score (lowest e-value) is selected as the corresponding class.


### Details

Tools involved:
- Biopython SeqIO
- HMMER3
- pyhmmer

So far we had hmm's built from Ribovision's per-class MSAs. Ideally, the HMMs are constructed on the fly:

- for each candidate polypeptide class:
    - from a repo of sequences of the given class we pick **N**==10 that are phylogenetically closest to the organism being classified. 
    - we construct an HMM based on these N sequences 
    


------------------

## Prelim. progress:

classifying a structure via hmms:
- given the organism in Q, grab the taxid
- for 
