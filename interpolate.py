import sys
import math
import Image

def main():
	file1 = raw_input("First file: ")
	file2 = raw_input("Second file: ")
	interpolationValue = input("Where in between? (between 0 and 1) ")

	im1 = Image.open(file1)
	im2 = Image.open(file2)

	check(im1, im2)

	pyramid = [(im1, im2)] # each level will have 2 images
	for i in range(1, math.log(im1.size[0])+2):
		print "Building layer", i
		currSize = im1.size[0]/(2**i)
		t1 = Image.new("L", (currSize, currSize))
		t2 = Image.new("L", (currSize, currSize))
		layer = (t1, t2)
		prevLayer = pyramid[-1]
		p1, p2 = prevLayer
		
		newpix1 = t1.load()
		newpix2 = t2.load()
		
		pix1 = p1.load()
		pix2 = p2.load()
		
		for x in range(currSize):
			for y in range(currSize):
				newpix1[x,y] = (pix1[x*2, y*2] + pix1[x*2+1, y*2] + pix1[x*2, y*2+1] + pix1[x*2+1, y*2+1])/4
				newpix2[x,y] = (pix2[x*2, y*2] + pix2[x*2+1, y*2] + pix2[x*2, y*2+1] + pix2[x*2+1, y*2+1])/4
		pyramid += [layer]
		t1.save("a"+str(currSize)+".bmp")
		t2.save("b"+str(currSize)+".bmp")
	
	# Pyramid is built, now to do the hard work
		
		
def check(im1, im2):
	# make sure dimensions match
	if im1.size != im2.size:
		print "images don't have same size"
		sys.exit()

	s = im1.size

	# checking for square image
	if (s[0] != s[1]):
		print "please use square images"
		sys.exit()

	# checking for power of 2
	if (s[0] & -s[0]) != s[0]:
		print "please use images with power of 2 width and height"
		sys.exit()

if __name__ == "__main__":
    main()