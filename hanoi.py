# -*- coding: utf-8 -*-
"""
Created on Mon May 20 19:24:11 2024

@author: joanv
"""

import numpy as np 


def hanoi(n, r1 = 'A',r2 ='B', r3='C'):
    rods = []
    if n==1:
        step = [[n,r3]]
        return step
    else:
        rods = {r1,r2,r3}
        before = hanoi(n-1,r1,r2,r3)
        last_rod = before[-1][1]
        diff = rods - {r1,last_rod}
        mid_rod = diff.pop()
        midstep = [[n,mid_rod]]
        if n%2==0:
            after = hanoi(n-1,last_rod,r1,mid_rod) 
        else:
            after = hanoi(n-1,last_rod,mid_rod,r1)
            
        return before + midstep + after
              