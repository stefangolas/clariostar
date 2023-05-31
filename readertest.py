# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 12:46:25 2023

@author: stefa
"""

from clariostar import ClarioStar, PlateData

def sim_test():
    reader = ClarioStar(simulating = True)
    output = reader.run_protocols(['a','b'])
    print(output[0].path)