# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 14:04:46 2024

@author: joanv
"""

import matplotlib.pyplot as plt
import numpy as np
# Create the x and y values for the step function
patchsize = 0.5
x = np.array([0, 1, 2, 3, 4, 5])
y = np.array([0, 1, 2, 3, 2, 1])
dy = np.array([0.5,0.2,0.4,0.25,0.8,0.1])

dxleft = np.abs(x[0] - x[1])/2
dxright = np.abs(x[-1] - x[-2])/2 
xleft = np.array([x[0]-dxleft,x[0]])
xright = np.array([x[-1],x[-1] + dxright])
yleft = np.array([y[0],y[0]])
yright = np.array([y[-1],y[-1]])


# Create the step function
plt.step(x, y, where='mid',color='black')
plt.plot(xleft, yleft)
plt.plot(xright, yright)
plt.fill_between(x, y-dy,y+dy,step='mid' )
plt.fill_between(xleft, yleft-dy[0],yleft+dy[0],step='pre' )
plt.fill_between(xright, yright-dy[-1],yright+dy[-1],step='post' )
plt.scatter(x,y,color='black')
# Set the title and labels for the plot
plt.title('Step Function with Different Widths')
plt.xlabel('X')
plt.ylabel('Y')

# Show the plot
plt.show()