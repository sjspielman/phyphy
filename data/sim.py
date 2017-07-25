from pyvolve import *


t = read_tree(file = "test.tre")
m = Model("MG94", {"omega":  np.arange(0.1, 1.4, 0.2)})
p = Partition(size = 10, models = m)
e = Evolver(tree = t, partitions=p)
e(seqfile = "seqs.fasta", infofile=None, ratefile=None)

m2 = Model("WAG")
p2 = Partition(size = 10, models = m2)
e2 = Evolver(tree = t, partitions=p2)
e2(seqfile = "aa.fasta", infofile=None, ratefile=None)