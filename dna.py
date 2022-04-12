'''
Header:
    bits:       bool
    rna:        bool

    desc size:  uint
Content:
    desc:       str
    data:       bytes
'''

from Bio import SeqIO
from Bio.Seq import Seq

key = {
    '?': 0b0000, 
    'A': 0b0001, 
    'T': 0b0010, 
    'G': 0b0011, 
    'C': 0b0100, 
    '[': 0b0101, 
    '/': 0b0110, 
    ']': 0b0111,
    0:  0b1000,
    1:  0b1001,
    2:  0b1010,
    3:  0b1011,
    4:  0b1100,
    5:  0b1101,
    6:  0b1110,
    7:  0b1111,
}

def seq_to_dna(sequence: Seq, filename: str):
    '''
    Bio.Seq -> *.dna
    '''

    from collections import Counter

    pairs = Counter([sequence[i:i+2] for i in range(len(sequence) + 1)])
    common = [pair for pair, _ in pairs.most_common(8) if len(pair) > 1]

    # count positions of pairs, then collapse them later?

    data = []
    # carry = XXXX????
    i = 0; carry = None
    while i < len(sequence) - 1:
        pair = sequence[i:i+2]

        try:
            value = key[common.index(pair)]
            i += 1
        except(ValueError):
            value = key[sequence[i]]
        
        if carry is not None:
            data.append((carry << 4) | value)
            carry = None
        else:
            carry = value

        i += 1

    # handle last base, use ] as end sentinel
    if i == len(sequence) - 1:
        if carry is not None:
            data.append((carry << 4) | key[sequence[i]])
            data.append((key[']'] << 4))
        else:
            data.append((key[sequence[i]] << 4) | key[']'])

    else:
        if carry is not None:
            data.append((carry << 4) | key[']'] << 4)


    with open(filename, 'wb') as f:
        # write mapping
        for i in range(8):
            if i < len(common):
                f.write(bytes(common[i]))
            else:
                f.write(bytes('**', 'utf-8'))

        f.write(bytes(data))

        print(f'dna bytes:      {f.tell()}')


def dna_to_seq(filename: str):
    '''
    *.dna -> Bio.Seq
    '''

    data = []
    with open(filename, 'rb') as f:
        # read mapping
        ikey = {v: k for k, v in key.items()}
        ikey.update({0b1000 + i: f.read(2).decode('utf-8') for i in range(8)})

        byte = f.read(1)
        while byte:
            byte = int.from_bytes(byte, 'big')
            data.append(ikey[(byte >> 4) & 0b1111])
            data.append(ikey[byte & 0b1111])
            byte = f.read(1)

        if data[-1] != ']':
            data.pop()
        data.pop()

    return Seq(''.join(data))



for seq_record in SeqIO.parse("covid.fasta", "fasta"):
    sequence = seq_record.seq
    # sequence = Seq('ATTAAAAAA')
    print(f'fasta bytes:    {len(sequence)}')
    seq_to_dna(sequence, 'covid.dna')
    decoded = dna_to_seq('covid.dna')
    print(f'lossless:       {sequence == decoded}')

    print(repr(sequence))
    print(repr(decoded))