# -*- coding: utf-8 -*-
"""
Created on Wed May  3 10:11:59 2023

@author: stefa
"""

from clariostar import ClarioStar
from clariostar import PlateData

reader_int = ClarioStar(simulating = True)
reader_int.plate_out(block=True)

plate_datas = reader_int.run_protocols(['luminescence'], plate_id_1='plate_1')

plate_data = plate_datas[0]

print(filename = plate_data.path)
print(plate_id = plate_data.header.plate_ids[0])
print(plate_data.value_at(1,1))