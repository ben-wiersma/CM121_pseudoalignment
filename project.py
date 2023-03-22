from Bio import SeqIO




'''
for seq_transcriptome in SeqIO.parse("chr11_transcriptome.fasta", "fasta"):
    #print(seq_transcriptome.id)
    print(segment(seq_transcriptome.seq, 10))
    #print(len(seq_transcriptome))
'''



    


print("running")


def pseudoalignment(transcriptome, k):

    transcriptome_dict = dict()

    for seq_transcriptome in SeqIO.parse("chr11_transcriptome.fasta", "fasta"):
        for i in range(0, len(seq_transcriptome) - (k + 1)):

            seq_str = str(seq_transcriptome.seq)[i:(k + i)]

            if seq_str not in transcriptome_dict:
                transcriptome_dict[seq_str] = {seq_transcriptome.id}
            else:
                transcriptome_dict[seq_str].add(seq_transcriptome.id)


    print(list(transcriptome_dict.items())[:10])
    
pseudoalignment("bruh", 30)

'''

output = list()
string = "hellothere"

for i in range(0, len(string) - 2):
    output.append(string[i:(i+3)])

print(output)

'''