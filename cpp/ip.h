#ifndef IP_H
#define IP_H


#include "common.h"
#include "image.h"


/*
 * IMPORTANT - DO NOT CHANGE THE INTERFACES DEFINED HERE - IMPORTANT
 *
 * It's very important that the interface to your program does not
 * change from what is provided, so that automated testing scripts
 * will work.  If you add additional requests for input, these scripts
 * will no longer work and it will negatively affect your grade.  Each
 * method has sufficient information to carry out the required task.
 *
 * The only method you are may get more user input in is ip_warp().
 * This option will not be tested automatically, since everyone's
 * will be different, so you may ask for whatever input is necessary.
 * To support multiple warps, have ip_warp() print a menu of options.
 *
 * Of course, you may add whatever other functions you may need.
 */


/*
 * definitions for the required sampling techniques
 */
#define I_NEAREST	0
#define I_BILINEAR	1
#define I_GAUSSIAN	2
#define KERNEL(W,H) (k[(W)+ksize*(H)])

Image* ip_interpolate (const char* imageName1, const char* imageName2, double inter);
double gradientx(Image* img, int x, int y);
double gradienty(Image* img, int x, int y);
double correspondence(Image* A, Image* B, int p1_x, int p1_y, int p2_x, int p2_y);
double stdDev(Image* image, int x, int y);
Image*  interpolate (Image* i1, Image* i2, double alpha);
Image*  convolute (Image* src, double* k, int ksize);
double  resample (Image* src, double w, double h, int chan, int mode, int size, double sigma);
Image*	ip_blur_box (Image* src, int size);
Image*	ip_blur_gaussian (Image* src, int size, double sigma);
Image*	ip_blur_triangle (Image* src, int size);

Image*	ip_brighten (Image* src, double alpha);
Image*	ip_contrast (Image* src, double alpha);
Image*	ip_composite (Image* src, const char* imageName, 
					  const char* maskName);
Image*	ip_crop (Image* src, int x0, int y0, int x1, int y1);
Image*	ip_edge_detect (Image* src);
Image*	ip_extract (Image* src, int channel);
Image*	ip_grey (Image* src);
Image*	ip_grey_chan (Image* src);
Image*	ip_invert (Image* src);
Image*	ip_noisify (Image* src, double alpha);

Image*	ip_quantize_simple (Image* src, int bitsPerChannel);
Image*	ip_quantize_ordered (Image* src, int bitsPerChannel);
Image*	ip_quantize_fs (Image* src, int bitsPerChannel);

Image*	ip_rotate (Image* src, double theta, int x, int y, int mode, 
		int size, double sigma);
Image*	ip_saturate (Image* src, double alpha);
Image*	ip_scale (Image* src, double x, double y, int mode, int size, 
		double sigma);
Image*	ip_threshold (Image* src, double cutoff);
Image*	ip_warp (Image* src, int warpF);


#endif // IP_H
