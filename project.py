from Bio import SeqIO
import pandas as pd

print("running")


def pseudoalignment(k):

    #define global variables

    transcriptome_dict = dict()

    equivalence_dict = dict()

    data = pd.DataFrame(columns=["counts","number of items", "isoforms"])


    #create transcriptome hash table

    for seq_transcriptome in SeqIO.parse("chr11_transcriptome.fasta", "fasta"):

        for i in range(0, len(seq_transcriptome) - k + 1):

            seq_str = str(seq_transcriptome.seq)[i:(k + i)] # create k-mer

            seq_comp = str(seq_transcriptome.seq.reverse_complement())[i:(k + i)] # create reverse complement k-mer

            if seq_str not in transcriptome_dict:
                transcriptome_dict[seq_str] = {seq_transcriptome.id}
            else:
                transcriptome_dict[seq_str].add(seq_transcriptome.id)

            if seq_comp not in transcriptome_dict:
                transcriptome_dict[seq_comp] = {seq_transcriptome.id}
            else:
                transcriptome_dict[seq_comp].add(seq_transcriptome.id)

    
    # map equivalence classes
    
    for seq_read in SeqIO.parse("reads.fasta", "fasta"):

        equivalence = set()

        for i in range(0, len(seq_read) - k + 1):

            seq_str = str(seq_read.seq)[i:(k + i)]

            if seq_str in transcriptome_dict:
                
                if(len(equivalence) == 0 ):
                    equivalence = transcriptome_dict[seq_str]
                else:
                    equivalence = equivalence.intersection(transcriptome_dict[seq_str])
                
        equivalence_frozen = frozenset(equivalence)
        
        if equivalence_frozen not in equivalence_dict:
            equivalence_dict[equivalence_frozen] = 1
        else:
            equivalence_dict[equivalence_frozen] += 1


    # add data to dataframe
    
    for isoform, count in equivalence_dict.items():
        data.loc[len(data)] = [count, len(isoform), ','.join(isoform)]

    # save data to tsv

    data.sort_values(by=['number of items']).to_csv('data.tsv', sep='\t', header=False, index=False)





    
pseudoalignment(30)
