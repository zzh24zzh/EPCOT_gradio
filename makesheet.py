import argparse,os
from zipfile import ZipFile
import xlsxwriter
import shutil
def parser_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputFileName', '-i',
                          help='Input prediction file name',
                          required=True)
    parser.add_argument('--outFileName', '-o',
                       help='Output file name',
                       required=True)
    args = parser.parse_args()
    return args

def bedpe2hic(filename,out_dir):
    res=1000 if filename.startswith('microc') else 5000
    tmp_file=os.path.join(out_dir,filename)
    os.system('java -jar data/juicer_tools.jar pre -n -r '+str(res)+
              ' -d %s %s hg38'%(tmp_file,tmp_file.replace('.bedpe','.hic')))




def main():
    args = parser_args()
    out_dir=args.inputFileName.replace('.zip','')
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
        print('directory has been overwrited')
    os.mkdir(out_dir)
    with ZipFile(args.inputFileName, 'r') as zObject:
        zObject.extractall(path=out_dir)


    workbook = xlsxwriter.Workbook(args.outFileName)
    sheet1 = workbook.add_worksheet('Config')

    sheet1.write(0,0,'Var')
    sheet1.write(0, 1, 'Value')
    sheet1.write(1,0,'root')
    sheet1.write(1, 1, os.path.abspath(out_dir))
    sheet2=workbook.add_worksheet('Index')
    lines = [['Genome', 'Title', 'Type', 'Name', 'Values'],
             ['hg38','Epigenomic feature','track','A','B,C,D'],
             ['hg38','CAGE-seq','track','A','B,C,D'],
             ['hg38','Chromatin organization','track','A','B,C,D']]
    for i in range(len(lines)):
        for j in range(len(lines[0])):
            sheet2.write(i,j,lines[i][j])
    sheet3 = workbook.add_worksheet('Epigenomic feature')
    sheet4 = workbook.add_worksheet('CAGE-seq')
    sheet5 = workbook.add_worksheet('Chromatin organization')
    header=['shortLabel','uri','metaLink','longLabel']
    for i in range(len(header)):
        sheet3.write(0,i,header[i])
        sheet4.write(0, i, header[i])
        sheet5.write(0, i, header[i])
    epi_id,co_id=1,1
    for f in os.listdir(out_dir):

        if f.startswith('cage'):
            tmpline=['CAGE-seq',os.path.abspath(os.path.join(out_dir,f)),'None','EPCOT predicted CAGE-seq track']
            for i in range(len(header)):
                sheet4.write(1,i,tmpline[i])
        elif f.endswith('bedpe'):
            bedpe2hic(f, out_dir)
            tmpline = [f.replace('.bedpe',''), os.path.abspath(os.path.join(out_dir, f.replace('.bedpe','.hic'))), 'None', 'EPCOT predicted '+f.replace('.bedpe','')+' contact maps']
            for i in range(len(header)):
                sheet5.write(co_id, i, tmpline[i])
            co_id+=1
        else:
            tmpline = [f.replace('.bigWig', ''), os.path.abspath(os.path.join(out_dir, f)), 'None',
                       'EPCOT predicted ' + f.replace('.bigWig', '') + ' track']
            for i in range(len(header)):
                sheet3.write(epi_id, i, tmpline[i])
            epi_id += 1

    workbook.close()

main()