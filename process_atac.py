import pyBigWig
import numpy as np
import pickle,os
from scipy.sparse import csr_matrix,save_npz
from optparse import OptionParser

def atac_processing():
    usage = 'usage: %prog [options] <bigWig_file> <prefix>'
    parser = OptionParser(usage)
    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.error('Please justify the DNase-seq bigwig file path')
    else:
        dnase_file = args[0]
    if len(args)==2:
        prefix=args[1]
    else:
        prefix=''
    os.mkdir('atacseq')
    bw = pyBigWig.open(dnase_file)
    for chrom, length in bw.chroms().items():
        try:
            if chrom == 'chrX':
                chr = 'X'
            else:
                chr = int(chrom[3:])
        except Exception:
            continue
        temp = np.zeros(length)
        intervals = bw.intervals(chrom)
        for interval in intervals:
            temp[interval[0]:interval[1]] = interval[2]
        seq_length = length // 1000 * 1000
        save_npz('atacseq/%s_atac_chr%s.npz'%(prefix,chr),csr_matrix(temp.astype('float32')[:seq_length]))
if __name__ == '__main__':
    atac_processing()