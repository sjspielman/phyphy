import sys
sys.path.append("src/")
from phyphy_parser import *


p = HyPhyParser("json/ABSREL_alltest.json")
print( p.extract_model_tree("Rate Distributions", partition=0))

# print()
# 
# p = HyPhyParser("json/RELAX_testref_all.json", "RELAX")
# print( p.extract_model_rate_distributions("MG94xREV with separate rates for branch sets"))
# 
