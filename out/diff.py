from matplotlib.pylab import *
from numpy import outer
import itertools
import csv


GPU = []
with open('BD_test.csv', 'rb') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=' ')
	for row in spamreader:
		T = []
		for item in row:
			if item != '':
				item = int(item, 10)
				T.append(float(item))
				# T.append(item)
		GPU.append(T)


SEQ = []
with open('BD_test_Base.csv', 'rb') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=' ')
	for row in spamreader:
		T = []
		for item in row:
			if item != '':
				item = int(item, 10)
				# T.append(float(item))
				T.append(item)
		SEQ.append(T)
counter = 0.0
distErr = 0.0
for index in range(len(SEQ)):
	for colIndex in range(len(SEQ[0])):
		# row[item] = abs(row[item] - rowG[item])
		SEQ[index][colIndex] = abs(SEQ[index][colIndex] - GPU[index][colIndex])
		if(SEQ[index][colIndex] != 0):
			distErr += abs(SEQ[index][colIndex] - GPU[index][colIndex])
			counter += 1

# print SEQ
if counter != 0:
	print "Number of incorrect Cells:"
	print counter
	print "Error: "
	print counter / (len(SEQ)*len(SEQ))
	print "Dist Err:"
	print distErr / counter
	figure(3)
	imshow(SEQ, interpolation='none')
	imshow(SEQ)
	grid(True)
	show()
else:
	print "There is no error in this simulation"
