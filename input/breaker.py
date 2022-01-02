#python3 breaker.py lhg22L20MC5x.fa 2
import sys
inputFastaFile = sys.argv[1]
n = int(sys.argv[2])
fhand = open(inputFastaFile,'r')
seqs = fhand.read().split('>')
seq_len = int(len(seqs)/ n) + 1 

for i in range(n):
    fhand2 = open(str(inputFastaFile+'_'+str(i)+'.fa'),'w')
    for j in range(seq_len*i,seq_len*(i+1),1):
        if j == len(seqs):
            break
        if seqs[j].strip() == '':
            continue
        fhand2.write('>'+seqs[j])
    fhand2.close()
fhand.close()