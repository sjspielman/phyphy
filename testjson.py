import sys
sys.path.append("src/")
from phyphy_parser import *
from phyphy_runner import *



h = FEL(data = "json/lysin.nex", genetic_code = "Invertebrate mtDNA")
h.run_analysis()
p = HyPhyParser("json/FEL_1partition_testbg.json")
print( p.reveal_fields())

# print()
# 
# p = HyPhyParser("json/RELAX_testref_all.json", "RELAX")
# print( p.extract_model_rate_distributions("MG94xREV with separate rates for branch sets"))
# 
