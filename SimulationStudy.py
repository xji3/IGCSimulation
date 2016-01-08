from IGCsimulator import *
import argparse

def main(args):
    IGC_threshold = args.IGC_threshold
    IGC_g_list = [0.01, 0.05, 0.1, 0.5]
    IGC_q_list = [0.01, 0.02, 0.05, 0.1]
    IGC_threshold_list = [i * 0.1 for i in range(10)]

    blen = 50.0

    num_exon = 100
    num_intron = 300
    # x_exon and x_intron from MG94_YBR191W_YPL079W_nonclock_save.txt
    x_exon = [0.49249355302375575, 0.60985035555249456, 0.42155795722934408, 8.1662933909645563, 0.092804167727196338]
    x_IGC = [0.05, 0.02, 0.8]  # These values vary for the simulation study
    x_intron = [0.49249355302375575, 0.60985035555249456, 0.42155795722934408, 8.1662933909645563]


    pi = [0.3003473684188488, 0.21394338106802796, 0.19214618460490698, 0.2935630659082163]
    div_limit = 1 - sum([i **2 for i in pi])
    for IGC_g in IGC_g_list:
        for IGC_q in IGC_q_list:
            for replicate in range(1, 6):
                print 'IGC_g = ', IGC_g, 'IGC_q = ', IGC_q, 'IGC_threshold = ', IGC_threshold, ' replicate = ', replicate
                x_IGC = [IGC_g, IGC_q, IGC_threshold]
                log_file = './logs/log_g_' + str(IGC_g) + '_q_' + str(IGC_q) + '_threshold_' + str(IGC_threshold) + '_rep_' + str(replicate) + '.log'
                div_file = './logs/div_g_' + str(IGC_g) + '_q_' + str(IGC_q) + '_threshold_' + str(IGC_threshold) + '_rep_' + str(replicate) + '.log'
                test = OneBranchIGCSimulator(blen, num_exon, num_intron, x_exon, x_IGC, x_intron, log_file, div_file)
                try:
                    test.sim_one_branch(test.initial_seq, blen)
                except:
                    print "Failed"
                    test.add_final_seq()
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--IGCthreshold', dest = 'IGC_threshold', required = True, help = 'IGC threshold')
    main(parser.parse_args())

##    IGC_threshold_list = [i * 0.1 for i in range(10)]
##    sh_line = 'sbatch -o CSC530-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/'
##
##    with open('./Run_simulation.sh', 'w+') as f:
##        f.write('#!/bin/bash' + '\n')
##        for IGC_threshold in IGC_threshold_list:
##            f.write(sh_line + 'IGC_threshold_' + str(IGC_threshold) + '.sh\n' )
##            with open('./ShFiles/' + 'IGC_threshold_' + str(IGC_threshold) + '.sh', 'w+') as g:
##                g.write('#!/bin/bash' + '\n')
##                g.write('python SimulationStudy.py --IGCthreshold ' + str(IGC_threshold) + '\n')
                

