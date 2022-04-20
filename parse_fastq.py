'''
0101:
0110:
0111:
1000:
1001:
1010:
1011:
1100:
1101:
1110:
1111:
description \n
varint read length
[
    XXXX[SSSSSSSS]?
    ...
]
description \n
varint read length
[
    ...
]
'''

import sys

lookup = {
    '?':    0b0000, 
    'A':    0b0001, 
    'T':    0b0010, 
    'G':    0b0011, 
    'C':    0b0100, 
    0b0000: '?', 
    0b0001: 'A',
    0b0010: 'T',
    0b0011: 'G',
    0b0100: 'C'
}

def fastq_to_dna(infile: str, outfile: str):
    from collections import Counter
    from itertools import zip_longest

    seqids = []
    reads = []

    pairs = Counter()

    # via https://stackoverflow.com/a/1657385
    with open(infile) as f:
        for seqid, bases, _, scores in zip_longest(*[f]*4):
            seqids.append(seqid)
            seq = list(zip(bases[:-1], scores[:-1]))
            reads.append(seq)
            pairs.update(seq)

        print(f'uncompressed filesize:  {f.tell()} bytes')

    table = {pair[0]: 0b0101 + i for i, pair in enumerate(pairs.most_common(10)) if len(pair[0]) > 1}
    # print(pairs.most_common(10))

    # parse reads
    data = []
    for read in reads:
        read_data = []
        read_data.append(len(read))
        # TODO Change thesee to read data

        carry = None
        for base, score in read:
            if (base, score) in table:
                value = table[(base,score)]
                if carry is None:
                    carry = value
                else:
                    read_data.append((carry << 4) | value)
                    carry = None
            else:
                value = lookup[base]
                score = ord(score[0])
                if carry is None:
                    read_data.append((value << 4) | (score >> 4))
                    carry = score & 0b1111
                else:
                    read_data.append((carry << 4) | value)
                    read_data.append(score)
                    carry = None

        if carry is not None:
            read_data.append(carry << 4)

        data.append(read_data)

    with open(outfile, 'wb') as f:
        # write mapping
        for base, score in table:
            f.write(bytes(base + score, 'utf-8'))

        # fill empty mapping spaces (if necessary)
        for i in range(len(table) - 10):
            f.write(bytes('**', 'utf-8'))

        # write reads
        for seqid, read in zip(seqids, data):
            f.write(bytes(seqid, 'utf-8'))
            f.write(bytes(read))

        print(f'compressed filesize:    {f.tell()} bytes')


def dna_to_fastq(infile: str, outfile: str):
    with open(infile, 'rb') as f:
        # read mapping
        table = {0b0101 + i: f.read(2).decode('utf-8') for i in range(10)}
        
        # decompress
        seqids = []
        reads = []

        seqid = f.readline()
        # FIXME: we do not have a set size we are writing for this read length
        length = int.from_bytes(f.read(1), 'big')
        read = []
        carry = None
        for _ in range(length):
            if carry is None:
                byte = int.from_bytes(f.read(1), 'big')
                
                if (byte >> 4) in table:
                    read.append(tuple(table[byte >> 4]))
                    carry = byte & 0b1111
                else:
                    base = lookup[(byte >> 4)]
                    nb = int.from_bytes(f.read(1), 'big')

                    value = chr((byte << 4) | (nb >> 4))
                    read.append((base, value))
                    carry = nb & 0b1111
            else:
                if carry in table:
                    read.append(tuple(table[carry]))
                else:
                    base = lookup[carry]
                    nb = int.from_bytes(f.read(1), 'big')

                    read.append((base, chr(nb)))
                
                carry = None

        for base, score in read:
            # print(base, score)
            pass

fastq_to_dna(sys.argv[1], sys.argv[2])
dna_to_fastq(sys.argv[2], 'test.fastq')