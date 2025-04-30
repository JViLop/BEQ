# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 19:22:18 2024

@author: joanv
"""

import multiprocessing as mp

def square(x, q):
    q.put(x * x)

qout = mp.Queue()
processes = [mp.Process(target=square, args=(i, qout))
             for i in range(2, 10)]

for p in processes:
    p.start()

for p in processes:
    p.join()

result = [qout.get() for p in processes]
print(result)
# if __name__ == "__main__":
#     qout = mp.Queue()
#     processes = [mp.Process(target=square, args=(i, qout))
#                  for i in range(2, 10)]

#     for p in processes:
#         p.start()

#     for p in processes:
#         p.join()

#     result = [qout.get() for p in processes]
#     print(result)

        
    