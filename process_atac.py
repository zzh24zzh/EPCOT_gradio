import pyBigWig
import numpy as np
import pickle,os,argparse
from scipy.sparse import csr_matrix,save_npz
from optparse import OptionParser

def parser_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputFileName', '-i',
                          help='Input prediction file name',
                          required=True)
    parser.add_argument('--prefix', '-p',default='',help='prefix')
    parser.add_argument('--out_dir', '-o',
                       help='Output directory',
                       required=True)
    args = parser.parse_args()
    return args


def atac_processing():
    args=parser_args()
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)

    bw = pyBigWig.open(args.inputFileName)
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
        save_npz(os.path.join(args.out_dir,'%s_atac_chr%s.npz'%(args.prefix,chr)),csr_matrix(temp.astype('float32')[:seq_length]))
if __name__ == '__main__':
    atac_processing()