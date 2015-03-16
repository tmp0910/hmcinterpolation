About the path implementation described in the paper.

The method is based on there being a path at every pixel--and describes the transition from a pixel p in an image A to the same pixel p in a second image B.  Then we can move the pixels along a specified path and copy the gradient values at the pixel position in each image.  What we want to aim to do is to define in where one pixel is in two images so that we may be able to define a path between the two positions.  The path framework described in the paper is described as "an inverse optical flow," which allows for the avoidance of holes, and instead goes through each intermediate pixel, ensuring that all possible holes are avoided.  It also assumes that the movement of the point takes place is always a straight line--the limitation, thus, is that curved motion cannot be accounted for; however, this shouldn't be too much of a problem as long as the images we deal with are not overly complex.