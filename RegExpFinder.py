import fastaparser
import re
from os import listdir
import subprocess


def search():
    pattern = input('Input the regular expression\n')
    pattern = re.compile(pattern)
    output = {}
    file_list = listdir()
    file_list = [f for f in file_list if '.fasta' in f]
    for fasta_file in file_list:
        count = 0
        with open(fasta_file) as f:
            genome = fastaparser.Reader(f, parse_method='quick')
            for seq in genome:
                if 'chromosome' in seq.header:
                    chromosome =  seq.header.split('chromosome ')[1].split(',')[0]
                    point = 0
                elif 'mitochondrion' in seq.header:
                    chromosome = 'M'
                    point = 0
                elif '(.)' in seq.header:
                    chromosome = seq.header.split('chr')[1][0]
                    point = int(seq.header.split(':')[1].split('-')[0])
                for t in re.finditer(pattern, seq.sequence):
                    count = count + 1
                    if chromosome in output:
                        if t.group() in output[chromosome]:
                            output[chromosome][t.group()].append(t.start()+point)
                        else:
                            output[chromosome][t.group()] = [t.start()+point]
                    else:
                        output[chromosome] = {t.group():[t.start()+point]}
        print(str(count)+' matches in '+fasta_file)  
    f = open('position.bed', 'w')
    t = open('matches.txt', 'w')
    for chromosome in output:
        for seq in output[chromosome]:
            for p in output[chromosome][seq]:
                f.write('chr'+chromosome+'\t'+str(p)+'\t'+str(p+len(seq))+'\n')
                t.write(seq+'\n')
    f.close()
    t.close()
    return output 


out = search()
subprocess.run(['Rscript', 'Annotation-Go.R'])
