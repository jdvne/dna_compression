## Header:
'''
?   -> 0000 
A   -> 0001 
T/U -> 0010 
G   -> 0011 
C   -> 0100 
[   -> 0101 
/   -> 0110 
]   -> 0111
'''

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

from Bio import SeqIO
from collections import Counter

class dna():
    bits: bool
    rna: bool
    pairs: Counter
    data: bytes
    description_size: int
    description: str

for seq_record in SeqIO.parse("covid.fasta", "fasta"):
    sequence = seq_record.seq
    print(repr(sequence))
    print(len(sequence))

pairs = Counter()

last_base = ''
for base in sequence:
    if(last_base != ''):
        pair = ''.join([last_base, base])
        # print(pair)
        pairs.update([last_base + base])
    last_base = base

mc = [t[0] for t in pairs.most_common(8)]
print(mc)

data = bytearray()

# carry = XXXX????
carry = None
i = 0
while i < len(sequence) - 1:
    pair = ''.join(sequence[i:i+2])
    value = 0b0

    if pair in mc:
        value = key[mc.index(pair)]
        i += 1
    else:
        value = key[sequence[i]] 

    if carry is not None:
        data.append((carry << 4) | value)
        carry = None
    else:
        carry = value
    
    i += 1

# TODO remember to handle possible floating carry bit
# always end the file on an extra bracket ']'
# unmatched close-bracket = EOF
print(len(data)) # 69% increase in space efficiency... nice

with open('covid.dna', 'wb') as f:
    f.write(bytes(data))

