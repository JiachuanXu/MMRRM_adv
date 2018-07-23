#!/usr/bin/python3.5

import numpy as np
import matplotlib.pyplot as plt
 
xData = np.arange(0, 10, 1)
yData1 = xData.__pow__(2.0)
yData2 = np.arange(15, 61, 5)
er = []
er.append(3)
er.append(3)
er.append(3)
er.append(3)
er.append(3)
er.append(3)
er.append(3)
er.append(3)
er.append(3)
er.append(3)
plt.figure(num=1, figsize=(8, 6))
plt.title('Plot 1', size=14)
plt.xlabel('x-axis', size=14)
plt.ylabel('y-axis', size=14)
print(xData)
print(er)
plt.errorbar(xData, yData1, yData1, color='b', linestyle='--', marker='o', label='y1 data')
plt.errorbar(xData, yData2, er, color='r', linestyle='-', label='y2 data')
plt.legend(loc='upper left')
plt.savefig('images/plot1.png', format='png')
