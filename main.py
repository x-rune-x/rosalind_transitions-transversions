class FastaObj:
    def __init__(self, fasta_id, sequence):
        self.fasta_id = fasta_id
        self.sequence = sequence

    def get_id(self):
        return self.fasta_id

    def get_seq(self):
        return self.sequence

    def get_length(self):
        return len(self.sequence)


def create_fasta_list(file_name):
    file = open(file_name, "r")

    current_id = ""
    current_seq = ""
    fasta_list = []

    for line in file:
        if line[0] == ">" and current_id == "":
            current_id = line[1:].rstrip()
        elif line[0] == ">" and current_id != "":
            fasta_list.append(FastaObj(current_id, current_seq))
            current_id = line[1:].rstrip()
            current_seq = ""
        else:
            current_seq += line.rstrip()
    else:
        fasta_list.append(FastaObj(current_id, current_seq))

    file.close()
    return fasta_list


def transition_transversion(seq1, seq2):
    transition_num = 0
    transversion_num = 0
    base_types = {
        'A': 'purine',
        'G': 'purine',
        'C': 'pyrimidine',
        'T': 'pyrimidine'
    }
    for position in range(len(seq1)):
        if seq1[position] != seq2[position]:
            if base_types[seq1[position]] != base_types[seq2[position]]:
                transversion_num += 1
            else:
                transition_num += 1

    return transition_num/transversion_num


fastas = create_fasta_list('rosalind_tran.txt')

s1 = fastas[0].get_seq()
s2 = fastas[1].get_seq()

print(transition_transversion(s1, s2))
