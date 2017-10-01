import sys
sys.path.append("src/")
from phyphy_parser import *


p = HyPhyParser("json/BUSTED_1partition_testbg.json")
print( p.extract_model_tree("Nucleotide GTR", partition=0))

# print()
# 
# p = HyPhyParser("json/RELAX_testref_all.json", "RELAX")
# print( p.extract_model_rate_distributions("MG94xREV with separate rates for branch sets"))
# 
