from Bio import SeqIO
import subprocess, os

def get_read(log_file):
    with open(log_file, 'rb') as f:
        seqs = f.readlines()[-3:-1]

    seq_file = log_file.replace('.log', '.fa')
    with open(seq_file, 'w+') as g:
        g.write('>Simulated_' + log_file.replace('./logs/log_', '').replace('.log', '') + '\n')
        g.write(seqs[0][:-1] + seqs[1][:-1])

    ART_cmd = ['/Users/Xiang/GitFolders/CSC530_FinalProject/ART/art_illumina', '-sam', '-i', seq_file, '-p', '-l', '150' ,'-ss', 'HS25',
               '-f', '29', '-m', '200', '-s', '10', '-o', seq_file.replace('.fa', '_')]
    subprocess.check_output(ART_cmd)

def reassemble(read_files):
    PEAR_cmd = ['/usr/local/bin/pear', '-f', read_files[0], '-r', read_files[1], '-o', read_files[0].replace('1.fq', 'pear')]
    try:
        subprocess.check_output(PEAR_cmd)
        SeqIO.convert(read_files[0].replace('1.fq', 'pear.assembled.fastq'), 'fastq',
                  read_files[0].replace('1.fq', 'pear.assembled.fasta'), 'fasta')
    except:
        print "PEAR failed"

if __name__ == '__main__':
    IGC_g = 0.07
    IGC_q = 0.025
    IGC_threshold = 0.8

    for replicate in range(1, 51):
        log_file = '/Users/Xiang/GitFolders/CSC530_FinalProject/logs/log_g_' + str(IGC_g) + '_q_' + str(IGC_q) + '_threshold_' + str(IGC_threshold) + '_rep_' + str(replicate) + '.log'
        get_read(log_file)
        read_files = [log_file.replace('.log', '_1.fq'), log_file.replace('.log', '_1.fq')]
        reassemble(read_files)

