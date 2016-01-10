# IGC Simulation on a Tree
# Imports simulation on a branch from my CSC530 project with modification
# Xiang Ji
# xji3@ncsu.edu


# Unlike the one branch simulator, this one is quite simple
# consider only one exon with codon model

from IGCsimulator import OneBranchIGCSimulator, draw_from_distribution
from CodonGeneconFunc import *

class TreeIGCSimulator:
    def __init__(self, num_exon, newick_tree,
                 x_exon, log_folder, div_folder):
        self.newicktree   = newick_tree
        self.num_exon     = num_exon
        self.edge_to_blen = None
        self.node_to_num  = None
        self.num_to_node  = None
        self.edge_list    = None

        # Node sequence
        self.node_to_sequence = None
        self.node_to_sim      = {}

        # OneBranchIGCSimulator Related Parameters
        self.x_exon       = x_exon        # parameter values for exon model
        self.Model   = 'MG94'
        self.num_paralog  = 2

        self.distn         = None
        self.mut_Q         = None

        self.OneBranchSimulator = None

        # Constants for Sequence operations
        bases = 'tcag'.upper()
        codons = [a+b+c for a in bases for b in bases for c in bases]
        amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        
        self.nt_to_state    = {a:i for (i, a) in enumerate('ACGT')}
        self.state_to_nt    = {i:a for (i, a) in enumerate('ACGT')}
        self.codon_table    = dict(zip(codons, amino_acids))
        self.codon_nonstop  = [a for a in self.codon_table.keys() if not self.codon_table[a]=='*']
        self.codon_to_state = {a.upper() : i for (i, a) in enumerate(self.codon_nonstop)}
        self.state_to_codon = {i:a.upper() for (i, a) in enumerate(self.codon_nonstop)}
        if self.Model == 'MG94':
            self.pair_to_state  = {pair:i for i, pair in enumerate(product(self.codon_nonstop, repeat = 2))}
        self.state_to_pair  = {self.pair_to_state[pair]:pair for pair in self.pair_to_state}


        self.initiate()
        


    def initiate(self):
        self.get_tree()
        self.unpack_x()

    def unpack_x(self):
        if self.Model == 'MG94':
            # get stationary nucleotide distribution of codon model of exonic region
            pi_codon = self.unpack_x_exon()[0]
            distn_codon = [ reduce(mul, [pi_codon['ACGT'.index(b)]  for b in codon], 1) for codon in self.codon_nonstop ]
            distn_codon = np.array(distn_codon) / sum(distn_codon)
            self.distn = distn_codon
            self.mut_Q = self.get_MG94()
            
        self.node_to_sequence = {node:[] for node in self.node_to_num.keys()}

    def unpack_x_exon(self, log = False):
        if log:
            x_exon = np.exp(self.x_exon)
        else:
            x_exon = self.x_exon
            
        if self.Model == 'MG94':
            assert(len(self.x_exon) == 5)
            # %AG, %A, %C, kappa, omega
            pi_a = x_exon[0] * x_exon[1]
            pi_c = (1 - x_exon[0]) * x_exon[2]
            pi_g = x_exon[0] * (1 - x_exon[1])
            pi_t = (1 - x_exon[0]) * (1 - x_exon[2])
            pi = [pi_a, pi_c, pi_g, pi_t]

            kappa = x_exon[3]
            omega = x_exon[4]
            return [pi, kappa, omega]

    def get_MG94(self):
        Qbasic = np.zeros((61, 61), dtype = float)
        pi, kappa, omega = self.unpack_x_exon()
        for ca in self.codon_nonstop:
            for cb in self.codon_nonstop:
                if ca == cb:
                    continue
                Qbasic[self.codon_to_state[ca], self.codon_to_state[cb]] = get_MG94BasicRate(ca, cb, pi, kappa, omega, self.codon_table)
        expected_rate = np.dot(self.distn, Qbasic.sum(axis = 1))
        Qbasic = Qbasic / expected_rate
        return Qbasic

    def get_tree(self):
        tree = Phylo.read( self.newicktree, "newick")
        #set node number for nonterminal nodes and specify root node
        numInternalNode = 0
        for clade in tree.get_nonterminals():
            clade.name = 'N' + str(numInternalNode)
            numInternalNode += 1
        tree_phy = tree.as_phyloxml(rooted = 'True')
        tree_nx = Phylo.to_networkx(tree_phy)

        triples = ((u.name, v.name, d['weight']) for (u, v, d) in tree_nx.edges(data = True)) # data = True to have the blen as 'weight'
        T = nx.DiGraph()
        edge_to_blen = {}
        for va, vb, blen in triples:
            edge = (va, vb)
            T.add_edge(*edge)
            edge_to_blen[edge] = blen

        self.edge_to_blen = edge_to_blen

        # Now assign node_to_num
        leaves = set(v for v, degree in T.degree().items() if degree == 1)
        self.leaves = list(leaves)
        internal_nodes = set(list(T)).difference(leaves)
        node_names = list(internal_nodes) + list(leaves)
        self.node_to_num = {n:i for i, n in enumerate(node_names)}
        self.num_to_node = {self.node_to_num[i]:i for i in self.node_to_num}

        # Prepare for generating self.tree so that it has same order as the self.x_process
        nEdge = len(self.edge_to_blen)  # number of edges
        l = nEdge / 2 + 1               # number of leaves
        k = l - 1   # number of internal nodes. The notation here is inconsistent with Alex's for trying to match my notes.

        leaf_branch = [edge for edge in self.edge_to_blen.keys() if edge[0][0] == 'N' and str.isdigit(edge[0][1:]) and not str.isdigit(edge[1][1:])]
        out_group_branch = [edge for edge in leaf_branch if edge[0] == 'N0' and not str.isdigit(edge[1][1:])] [0]
        internal_branch = [x for x in self.edge_to_blen.keys() if not x in leaf_branch]
        assert(len(internal_branch) == k-1)  # check if number of internal branch is one less than number of internal nodes

        leaf_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        internal_branch.sort(key = lambda node: int(node[0][1:]))  # sort the list by the first node number in increasing order
        edge_list = []
        for i in range(len(internal_branch)):
            edge_list.append(internal_branch[i])
            edge_list.append(leaf_branch[i])
        for j in range(len(leaf_branch[i + 1:])):
            edge_list.append(leaf_branch[i + 1 + j])

        self.edge_list = edge_list


    def sim(self):
        self.sim_root()
        for edge in self.edge_list:
            print edge
            if edge in self.outgroup:
                rate_mat = self.Basic_mat
            else:
                rate_mat = self.Geneconv_mat

            IGC_mat = self.IGC_mat

            end_seq = []
            num_IGC = 0
            num_All = 0
            for site in self.node_to_sequence[edge[0]]:
                #print site
                site_seq, add_num_IGC, add_num_All = self.sim_one_branch(site, rate_mat, IGC_mat, self.edge_to_blen[edge])
                num_IGC += add_num_IGC
                num_All += add_num_All
                end_seq.append(site_seq)
            self.node_to_sim[edge[1]] = [end_seq, num_IGC, num_All]
            self.node_to_sequence[edge[1]] = end_seq

        for node in self.node_to_sim.keys():
            seq1 = ''.join([self.state_to_pair[i][0] for i in self.node_to_sim[node][0]])
            seq2 = ''.join([self.state_to_pair[i][1] for i in self.node_to_sim[node][0]])
            self.node_to_sequence[node] = (seq1, seq2)


    def sim_root(self):
        if self.Model == 'MG94':
            seq = draw_from_distribution(self.distn, self.num_exon, self.codon_nonstop)

        #self.node_to_sequence['N0'] = np.array([self.pair_to_state[(i, i)] for i in seq])
        self.node_to_sequence['N0'] = [''.join(seq), ''.join(seq)]
        self.node_to_sim['N0'] = [self.node_to_sequence['N0'], 0, 0]



if __name__ == '__main__':
    paralog1 = 'YDR418W'
    paralog2 = 'YEL054C'

    paralog = [paralog1, paralog2]
    newicktree = './YeastTree.newick'
    num_exon = 163
    log_folder = './log/'
    div_folder = './div/'

    x_exon = [0.49249355302375575, 0.60985035555249456, 0.42155795722934408, 8.1662933909645563, 0.092804167727196338]
    
    test = TreeIGCSimulator(num_exon, newicktree, x_exon, log_folder, div_folder)
    self = test

    test.sim_root()

    
