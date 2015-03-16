# Introduction #
We begin by first constructing the Gaussian pyramid (in which our original images are at the bottom of the pyramid and increasingly smaller and scaled versions are layered on top) for each of the input images. Starting with the topmost (smallest) layer, we begin our path construction for each pixel in our desired interpolated image.

For every pixel p in the interpolated image, we consider the path w denoted by the general motion vectors m\_A and m\_B (from input image A and B respectively) in which m\_A/(|m\_A|) = m\_B/(|m\_B|). (Because we begin with the topmost layer, the number of possible vectors are limited.) The transition points p\_A and p\_B are defined as p\_A = p + m\_A and p\_B = p - m\_B.

The cost C of a path w at pixel p is defined to be how well transition points p\_A and p\_B match, taking into account both intensity and gradients.  The equation for the matching formula is described by equations 4,5, and 6 in the main paper (equation 7 gives an alternative formula).  Because the cost C is determined on a per-pixel basis, we need a coherency cost V that takes into account the similarity of neighboring pixels (defined to be the pixels to the left, right, top, and bottom), given by equation 8.  We define the global energy function E to be the sum of the paths cost C for all pixels in the image and the pairwise coherency cost V for every pair of neighboring pixels (equation 9).

We minimize the energy function E with a fast iterative expansion algorithm based on graph cuts (see Boykov, Veksler, and Zabih for pseudocode) In actuality, there is no graph representation required. Each pixel has some number of potential paths; we begin by arbitrarily assigning one of these paths to each pixel. Then, in one iteration of the minimization algorithm, we go through each pixel p and select the "best" path w so that E(p, w) is minimized. (how is the neighboring pixel path p\_N decided for the correspondence function?) After going through all pixels, we find the total energy cost of the image and compare it to our original cost (before running the algorithm). If the new cost is less, then we run another iteration of the algorithm; otherwise we stop.

Once we choose the correct paths, we repeat the process for the next layer of the pyramid, doubling the length of each path and allowing the shift of one pixel (?) to account for the expanded size of the new layer. (See Section 4.1)

The forward and backward flows of a path are defined to be v\_A = p\_B - p\_A and v\_B = p\_A - p\_B respectively.  Using these flows, it is possible to detect inconsistent flow (occluded) regions in our images (if the forward and backward flows do not match).  Treating the occluded region as one, we see which flow should be assigned the entire region (based upon the flows of the included pixels in the region) and compute corrected paths.

After computing the correct path, we solve the 3D Poisson equation (equation 10) to compute the pixel values for the interpolated image.

Note: this way is slower, and also cannot handle more than 2 layers of occlusion.


# Details #

Add your content here.  Format your content with:
  * Text in **bold** or _italic_
  * Headings, paragraphs, and lists
  * Automatic links to other wiki pages