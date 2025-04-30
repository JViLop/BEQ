# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 10:34:39 2024

@author: joanv
"""

import multiprocessing as mp
import numpy as np
def square(x, q):
    q.put(x * x)


if __name__ == "__main__":
    
    qout = mp.Queue()
    processes = [mp.Process(target=square, args=(i, qout))
                 for i in [1,3,4,5,10,22,7]]
    
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    
    result = [qout.get() for p in processes]
    print(result)