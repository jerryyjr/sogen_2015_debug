"""
Creates snapshots of cell tissue colored based on concentration levels
Copyright (C) 2012 Ahmet Ay, Jack Holland, Adriana Sperlea

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import Image, ImageDraw
import struct
import sys
import shared

def main():
	print 'Reading command-line arguments...'
	args = sys.argv[1:]
	if len(args) == 3:
		cons_fname1 = args[0]
		cons_fname2 = args[1]
		directory = args[2]
	else:
		usage()
	
	print 'Reading concentrations file 1...'
	min_con1 = float('inf')
	max_con1 = 0
	cons_data1 = []
	if cons_fname1.endswith('.cons'): # Read ASCII file
		cons_file1 = shared.openFile(cons_fname1, 'r')
		width, height = map(lambda num: shared.toInt(num), cons_file1.readline().split(' ')) # The first line contains the width and height
		checkSize(width, height)
		for line in cons_file1:
			cons = map(lambda num: shared.toFlo(num), line.split(' ')[1:-1]) # Remove the time step column and newline when taking the concentrations
			for con in cons:
				min_con1 = min(min_con1, con)
				max_con1 = max(max_con1, con)
			
			cons_data1.append(cons)
	elif cons_fname1.endswith('.bcons'): # Read binary file
		cons_file1 = shared.openFile(cons_fname1, 'rb') # Read the file as a binary
		# The first two ints are the width and height
		width, = struct.unpack('i', cons_file1.read(4))
		height, = struct.unpack('i', cons_file1.read(4))
		checkSize(width, height)
		size = width * height
		cons1 = []
		cons_length1 = 0
		while True:
			con_str1 = cons_file1.read(4)
			if con_str1 == '': # While not EOF
				break;
			else:
				# There are width * height concentration floats per time step
				con, = struct.unpack('f', con_str1)
				min_con1 = min(min_con1, con)
				max_con1 = max(max_con1, con)
				cons1.append(con)
				cons_length1 += 1
				if cons_length1 == height:
					cons_data1.append(cons)
					cons1 = []
	else:
		usage()
		
	print 'Reading concentrations file 2...'
	min_con2 = float('inf')
	max_con2 = 0
	cons_data2 = []
	if cons_fname2.endswith('.cons'): # Read ASCII file
		cons_file2 = shared.openFile(cons_fname2, 'r')
		width, height = map(lambda num: shared.toInt(num), cons_file2.readline().split(' ')) # The first line contains the width and height
		checkSize(width, height)
		for line in cons_file2:
			cons = map(lambda num: shared.toFlo(num), line.split(' ')[1:-1]) # Remove the time step column and newline when taking the concentrations
			for con in cons:
				min_con2 = min(min_con2, con)
				max_con2 = max(max_con2, con)
			cons_data2.append(cons)
	elif cons_fname2.endswith('.bcons'): # Read binary file
		cons_file2 = shared.openFile(cons_fname2, 'rb') # Read the file as a binary
		# The first two ints are the width and height
		width, = struct.unpack('i', cons_file2.read(4))
		height, = struct.unpack('i', cons_file2.read(4))
		checkSize(width, height)
		size = width * height
		cons2 = []
		cons_length2 = 0
		while True:
			con_str2 = cons_file2.read(4)
			if con_str2 == '': # While not EOF
				break;
			else:
				# There are width * height concentration floats per time step
				con, = struct.unpack('f', con_str2)
				min_con2 = min(min_con2, con)
				max_con2 = max(max_con2, con)
				cons2.append(con)
				cons_length2 += 1
				if cons_length2 == height:
					cons_data2.append(cons)
					cons2 = []
	else:
		usage()
		
	
	print 'Creating the directory if necessary...'
	directory = shared.ensureDir(directory)
	if (directory[-1] != '/'):
		directory = directory + '/'
	
	cons_data = combine_cons(cons_data1,cons_data2, max_con1, min_con1, max_con2, min_con2)
	
	print 'Creating snapshots...'
	edge, size = findSizes(width, height) # Configure the hexagon edge and window size based on the grid size
	index = 0
	for line in cons_data:
		if (index % 10 == 0 and index >= 21000):
			
			plotHexagons(directory, size, index, line, edge, width, height)
		index += 1
	
	print 'Done. Your snapshots are stored in ' + directory

def checkSize(width, height):
	if width < 4 or height < 4 or (width % 2 != 0) or (height % 2 != 0):
		print 'The size of the tissue must be at least 4x4 and its width and height must be even numbers.'
		exit(2)
		
'''
# Finds the minimum and maximum concentrations given a properly formatted file
def findMinMax(filename):
	f = shared.openFile(filename, 'r')
	minCon = float('inf')
	maxCon = 0
	first = True
	for line in f:
		if not first:
			lineNums = line.split()
			for i in range(len(lineNums)):
				lineNums[i] = shared.toFlo(lineNums[i])
			num = min(lineNums[1:])
			if (num < minCon):
				minCon = num
			num = max(lineNums[1:])
			if (num > maxCon):
				maxCon = num
		first = False
	f.close()
	return minCon, maxCon

# Extracts the width and height from the concentrations file and parses each line into an array
def readData(filename):
	f = shared.openFile(filename, 'r')
	parsedFile = []
	index = 0
	width = height = 0
	for line in f:
		line = line.strip()
		if (index == 0):
			width, height = shared.widthAndHeight(line.split(' '), filename)
			if width < 4 or height < 4 or (width % 2 != 0) or (height % 2 != 0):
				print 'The size of the tissue must be at least 4x4 and its width and height must be even numbers.'
				exit(2)
		else:
			aux = line.split(' ')
			parsedFile.append(aux)
		index += 1
	f.close()
	return width, height, parsedFile
'''

# Finds the edge length and size of a hexagon, properly formatted for 300 pixels per inch (journal specifications)
def findSizes(width, height):
	edge = (float) (2 * 960) / (3 * width + 1)
	size = (960, int(edge * 1.73 * height))
	return edge, size
	
def combine_cons(cons_data1, cons_data2, max_con1, min_con1, max_con2, min_con2):
	cons_data=[]
	if (len(cons_data1) != len(cons_data2)):
		print "error\n"
		exit(2)
	else:
		for index in range(len(cons_data1)):
			scores=[]
			for r in range(len(cons_data1[index])):
				
				score1=int(((cons_data1[index][r] - min_con1) / (max_con1 - min_con1))*255)
				score2=int(((cons_data2[index][r] - min_con2) / (max_con2 - min_con2))*255)
				scores+=[format(score1*65536 + score2*256, '06X')]
			
			cons_data.append(scores)
		return cons_data
			
	
	


# Plot the hexagons for the given line of data
def plotHexagons(output, size, index, protLevels, edge, width, height): #bug
	center = (edge, float(edge * 1.73 / 2.0)) # place the center of the first hexagon
	im = Image.new('RGB', size, (255, 255, 255))
	draw = ImageDraw.Draw(im)
	drawGrid(draw, center, edge, width, height, protLevels)
	fix = fixDigits((index - 21000) / 10, 4)
	im.save(output + fix + '.png', 'PNG')
	del im

# Draws a grid of hexagons with the first hexagon's center at 'begin'
def drawGrid(draw, begin, edge, width, height, protLevels):
	center = begin	
	
	rad = float(edge * 1.73 / 2.0)
	for i in range(height):
		for j in range(width):
			# draw a half of the last cell
			if (i == 0 and j % 2 == 1):
				drawCell(draw, protLevels, (height - 1) * width + j, edge, (center[0], 0))
		
			cell = i * width + j
			drawCell(draw, protLevels, cell, edge, center)			

			if (j % 2 == 0):
				center = (center[0] + edge + edge / 2.0, center[1] + edge * 1.73 / 2.0) # move center to the right by 1.5 * edge and down by edge * sqrt(3) / 2
			else:
				center = (center[0] + edge + edge / 2.0, center[1] - edge * 1.73 / 2.0) # move center to the right by the same amount, but up instead of down
				
		center = (begin[0], begin[1] + edge * 1.73 * (i + 1))

# Converts a triplet of decimal integers to a triplet of hexadecimal integers
hexDigits = '0123456789abcdef'
def rgb(triplet):
	
    triplet = triplet.lower()
	
    return (hexDigits.index(triplet[0]) * 16 + hexDigits.index(triplet[1]),
            hexDigits.index(triplet[2]) * 16 + hexDigits.index(triplet[3]),
            hexDigits.index(triplet[4]) * 16 + hexDigits.index(triplet[5]))

# The gradients of red to use
redShades = [rgb('FEB4EF'), rgb('FEB4EF'), rgb('FE5A77'), rgb('FE2D3B'), rgb('FF0000'), rgb('BF0000'), rgb('7F0000'), rgb('3F0000'), rgb('000000'), rgb('FFFFFF')]
numShades = len(redShades)


# Draws an individual hexagon
def drawCell(draw, protLevels, cell, edge, center):
	 

	#if protLevel == 0:
	#	colorindex = numShades - 1
	#else:
	
	colorindex = protLevels[cell]
	#print colorindex + " color"
		#if colorindex < 0 or colorindex == 9:
			#print colorindex, protLevel, minCon, maxCon
	
	hex = hexagon(center, edge, rgb('fcc75e')) # border
	draw.polygon(hex.draw(), fill = hex.color)
	del hex
	
	hex = hexagon(center, edge - 2, rgb(colorindex)) # actual hexagon
	draw.polygon(hex.draw(), fill = hex.color)
	del hex

# Fix the digits of a number so that they are padded with 0's up to n
def fixDigits(index, n):
	zero = ''
	fix = str(index)
	while (len(fix) + len(zero) < n):
		zero += '0'
	fix = zero + fix
	return fix

# Contains information required to draw a hexagon
class hexagon():
	# Creates a hexagon with the center at self.center and the edge of size self.edge
	def __init__(self, center, edge, color): 
		self.center = center
		self.edge = edge
		self.color = color
	
	# Draws the hexagon centered at self.center
	def draw(self):
		half = float(self.edge / 2.0)
		rad = float(self.edge * 1.73 / 2.0)
		i = self.center[0]
		j = self.center[1]
		lines = [(i - half, j - rad), (i + half, j - rad), (i + self.edge, j), (i + half, j + rad), (i - half, j + rad), (i - self.edge, j)]
		return lines

def usage():
	print 'Usage: python plot-tissue-snapshots.py <file with concentration levels> <directory to store images>'
	exit(0)

main()

