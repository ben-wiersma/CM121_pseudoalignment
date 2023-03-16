from Bio import SeqIO
for seq_transcriptome in SeqIO.parse("chr11_transcriptome.fasta", "fasta"):
    print(seq_transcriptome.id)
    print(repr(seq_transcriptome.seq))
    print(len(seq_transcriptome))



def pseudoalignment(transcriptome, k):
    print("hello world")

