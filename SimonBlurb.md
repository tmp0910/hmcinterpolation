# Blurb #

We first create a gaussian pyramid for each of the input images. This allows us to find paths with arbitrary transition points without having an extremely large set of paths to consider. After generating the gaussian pyramids, we would start finding paths starting from the smallest layer of the pyramid. We pick paths for pixels based on a few metrics. A path defines two transition points between the two images for a pixel location. We constrain the two vectors that go towards the transition points to have the parallel, but opposite, directions. This is to ensure that movement of the transition points is uniform throughout the interpolation range. We want to minimize both correspondence (difference in the transition points in terms of gradient and color) and coherency (the difference in motion between the neighboring pixel paths). The paper proposes a graph cut method of minimizing the energy. Once we pick the correct path for a layer, we can move on to the next layer.

After completing in finding the path, we can then easily interpolate the new image by following the paths to find where the original pixel value should be. The paper proposes to use gradients instead of intensity for better results, especially in preserving edges.


# Details #

Add your content here.  Format your content with:
  * Text in **bold** or _italic_
  * Headings, paragraphs, and lists
  * Automatic links to other wiki pages