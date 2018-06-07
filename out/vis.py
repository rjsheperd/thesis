from matplotlib.pylab import *
import matplotlib.animation as animation
from numpy import outer
import csv
import numpy as np
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

A = []
with open(dir_path+'/final_tests.csv', 'rb') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',')
	for row in spamreader:
		T = []
		for item in row:
			if item != '':
				item = int(item, 10)
				T.append(float(item))
		A.append(T)

for row in range(len(A)):
	#print row
	for item in range(len(A[0])): 
		#print item, type(item)
		if A[row][item] == 32767:
			A[row][item] = 0	
		#print item, type(item)
Max = 0.0
Min = 1000
for row in A:
	for item in row:
		if item > Max and item != 32767:
			Max = item	
		if item != 0 and item < Min:
			Min = item

print ("min: " + str(Min)) 

A2 = np.zeros(np.array(A).shape)#np.random.rand(np.array(A).shape[0],np.array(A).shape[1])
fig = figure()

#im = imshow(A2, interpolation='none')
im = imshow(A2, interpolation='none',vmin=0,vmax=Max)
#pcolor(vmin=0, vmax=max)
def animate(i):
	#print i
	for row in range(len(A)):
		for col in range(len(A[0])):
			if A[row][col] == i and i != 0:
				A2[row][col] = A[row][col] * 2.0
				#print "A2: " + str(A2[row][col])
	
	im.set_array(A2) 

	return im,
#print A
print Max

#imshow(A2, interpolation='none')

ani = animation.FuncAnimation(fig, animate,interval=1,frames=600)
#imshow(A)
#grid(False)
show()
#print A
