#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 11:04:39 2019

@author: rosalyn
"""

import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt

fp1_table = pd.read_table("./fp1/output/tables/experiments_results_table.csv",
                          dtype ={"user_id": float})
fp2_table = pd.read_table("./fp2/output/tables/experiments_results_table.csv",
                          dtype ={"user_id": float})


