from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


a = SeqRecord(Seq("AAAACGT"), id="Alpha")
b = SeqRecord(Seq("AAA-CGT"), id="Beta")
c = SeqRecord(Seq("AAAAGGT"), id="Gamma")
align = MultipleSeqAlignment([a, b, c],
                             annotations={"tool": "demo"},
                             column_annotations={"stats": "CCCXCCC"})

# read in the fast file:
from Bio import AlignIO
aln:MultipleSeqAlignment = AlignIO.read("./ribovision-3j7z.fas", 'fasta')

# for i in aln:
#     print(i.seq)
  
print(aln.__len__())
print(aln[67])


