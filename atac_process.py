import argparse,os
import uuid
import pyBigWig
import numpy as np
import pickle
from scipy.sparse import csr_matrix
from optparse import OptionParser

def parser_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam',required = True)
    parser.add_argument('-p','--num_processers',default=6)
    args = parser.parse_args()
    return args

def atac_bwtonpz(atac_file):
    bw = pyBigWig.open(atac_file)
    signals = {}
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
        signals[chr] = csr_matrix(temp.astype('float32')[:seq_length])
    with open('atac_'+atac_file.replace('bigWig','pickle'),'wb') as f:
        pickle.dump(signals,f)

    print('The processed ATAC-seq file has been saved to ' + os.path.abspath(
        'atac_' + atac_file.replace('bigWig', 'pickle')))




def main():
    args=parser_args()

    fid=str(uuid.uuid4())
    os.system('samtools index %s'%args.bam)
    os.system('bamCoverage --bam %s -o %s --outFileFormat bigwig --normalizeUsing RPGC '
              '--effectiveGenomeSize 2913022398 --Offset 1 --binSize 1 --numberOfProcessors %s '
              '--blackListFileName data/black_list.bed'%(args.bam,fid+'.bigWig',args.num_processers))
    print('---------\n')

    atac_bwtonpz(fid+'.bigWig')
if __name__=="__main__":
    main()

