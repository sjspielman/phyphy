import sys
sys.path.append("src/")
from phyphy_parser import *


p = HyPhyParser("json/slac.json")
p.extract_csv("slac.csv")