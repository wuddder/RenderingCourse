import sys
from PIL import Image


####################
# class definition #
####################
class Point:
	def __init__(self, theX, theY):
		self.x = theX
		self.y = theY
	def getX(self):
		return self.x
	def getY(self):
		return self.y

#####################
# utility functions #
#####################

def print_usage():
	print >> sys.stderr, "Usage:\n$python " + __file__ + " <txt file path> [background color]"
	print >> sys.stderr, "(valid background color: BLACK or WHITE)"

######################
# main program entry #
#######################
if __name__ == "__main__":
	
	# check the arguments
	if len(sys.argv) not in [2, 3]:
		print_usage()
		sys.exit(1)

	# get background color if given
	FOREGROUND, BACKGROUND = None, None
	if len(sys.argv) == 3:
		if sys.argv[2] == "WHITE":
			FOREGROUND = 0
			BACKGROUND = 255
	if (FOREGROUND is None) or (BACKGROUND is None):
		FOREGROUND = 255
		BACKGROUND = 0
	
	# open the file
	try:
		f = open(sys.argv[1])
	except:
		print >> sys.stderr, "fatal: cannot open file at: " + sys.argv[1]
		sys.exit(1)
	
	# parse the file: 5 steps.
	# 1) setup variables
	point_count = -1
	points = []
	# 2) read the first line, which contains the number of points
	try:
		lineBuf = f.readline()
		point_count = int(lineBuf)
	except:
		print >> sys.stderr, "fatal: cannot get number of points at first line"
		sys.exit(1)
	# 3) read all the points
	lineBuf = f.readline()
	while lineBuf != "":
		parsed_list = lineBuf.lstrip("\t\r\n(").rstrip("\t\r\n)").split(", ")
		if len(parsed_list) == 2:
			try:
				theX = float(parsed_list[0])
				theY = float(parsed_list[1])
			except:
				print >> sys.stderr, "cannot convert to floating point number, this line: " + lineBuf.rstrip("\n")
			new_point = Point(theX, theY)
			points.append(new_point)
		else:
			print >> sys.stderr, "cannot recognize this line: " + lineBuf.rstrip("\n")
		lineBuf = f.readline()
	# 4) checking
	if len(points) != point_count:
		print >> sys.stderr, "number of points incorrect: " + len(points) + " / " + point_count
	# 5) close file
	f.close()
	
	# drawing
	im = Image.new("L", (1024, 1024), BACKGROUND) # Image.new(mode, size[, color])
	imPixels = im.load()
	for p in points:
		# draw the point
		drawX = int(p.getX() * 1024)
		drawY = int(p.getY() * 1024)
		imPixels[drawX, drawY] = FOREGROUND # draw as white (since default is black)
	im.save("Poisson.bmp")
	
	

















