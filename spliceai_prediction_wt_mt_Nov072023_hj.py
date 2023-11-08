from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode
import numpy as np
import pandas as pd
import sys

### usage ###
#python spliceai_prediction_wt_mt_Nov072023_hj.py [chrom] [gene] [position] [ref] [alt] [txstart from gencode] [txend from gencode]

# ---------------------------------------------------------
# argument.
chrom = str(sys.argv[1])
gene = str(sys.argv[2])
mut_pos = int(sys.argv[3])
ref = str(sys.argv[4])
alt = str(sys.argv[5])
txstart = int(sys.argv[6])
txend = int(sys.argv[7])

# ----------------------------------------------------------
# file name.
fasta_file = "./" + gene + "/" + "human_g1k_v37_chr" + chrom + "_rmnewline.fasta"
wt_file = "./" + gene + "/" + "spliceai_wt_" + gene + "_hj.txt"
mt_file = "./" + gene + "/" + "spliceai_mt_" + gene + "_hj.txt"

# ----------------------------------------------------------
# WT sequence
# ----------------------------------------------------------
# open file.
with open(fasta_file) as f:
    chrom = f.read()

chrom = list(chrom) # make list.

txstart_extend = txstart - 5000
txend_extend = txend + 5000
input_sequence = chrom[txstart_extend - 1:txend_extend]
input_sequence = ''.join(input_sequence)

context = 10000
paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(resource_filename('spliceai', x)) for x in paths]
x = one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :]
y = np.mean([models[m].predict(x) for m in range(5)], axis=0)

acceptor_prob = y[0, :, 1]
donor_prob = y[0, :, 2]
acceptor_prob = acceptor_prob.tolist()
donor_prob = donor_prob.tolist()

sequence = list(input_sequence)
position = list(range(txstart_extend, txend_extend+1))

tb = pd.DataFrame({'position' : position, 'sequence' : sequence, 'acceptor_prob' : acceptor_prob, 'donor_prob' : donor_prob})
tb.to_csv(wt_file, index = False, mode = 'w', header = True, sep = "\t")

# ----------------------------------------------------------
# MT sequence 
# ----------------------------------------------------------
# open file.
with open(fasta_file) as f:
    chrom = f.read()

chrom = list(chrom) # make list.

# sanity check using ref. 
if chrom[mut_pos-1] != ref :
    print("reference nucleotides from argument and fasta file are not consistent.")
    exit()

chrom[mut_pos-1] = alt # substitution.

txstart_extend = txstart - 5000
txend_extend = txend + 5000
input_sequence = chrom[txstart_extend - 1:txend_extend]
input_sequence = ''.join(input_sequence)

context = 10000
paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(resource_filename('spliceai', x)) for x in paths]
x = one_hot_encode('N'*(context//2) + input_sequence + 'N'*(context//2))[None, :]
y = np.mean([models[m].predict(x) for m in range(5)], axis=0)

acceptor_prob = y[0, :, 1]
donor_prob = y[0, :, 2]
acceptor_prob = acceptor_prob.tolist()
donor_prob = donor_prob.tolist()

sequence = list(input_sequence)
position = list(range(txstart_extend, txend_extend+1))

tb = pd.DataFrame({'position' : position, 'sequence' : sequence, 'acceptor_prob' : acceptor_prob, 'donor_prob' : donor_prob})
tb.to_csv(mt_file, index = False, mode = 'w', header = True, sep = "\t")

