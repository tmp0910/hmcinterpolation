#include "ip.h"
#include "main.h"
#include <algorithm>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <queue>

Image* ip_interpolate (const char* imageName1, const char* imageName2, double inter)
{
	cout << imageName1 << endl;
	cout << imageName2 << endl;
	cout << inter << endl;
	
	Image* i1 = new Image(imageName1);
	Image* i2 = new Image(imageName2);
	
	// Create gaussian pyramid
	vector<Image*> p1;
	vector<Image*> p2;
	p1.push_back(i1);
	p2.push_back(i2);
	for (int i = 1; i < (log(i1->getHeight())/log(2)); i++) {
		cout << i << endl;
		int currSize = i1->getHeight()/pow(2, i);
		Image* t1 = new Image(currSize, currSize, 1);
		Image* t2 = new Image(currSize, currSize, 1);
		
		Image* prev1 = p1.back();
		Image* prev2 = p2.back();
		
		for (int x = 0; x < currSize; x++) {
			for (int y = 0; y < currSize; y++) {
				t1->setPixel(x, y, 0, (prev1->getPixel(x*2, y*2, 0)+prev1->getPixel(x*2+1, y*2, 0)+prev1->getPixel(x*2, y*2+1, 0)+prev1->getPixel(x*2+1, y*2+1, 0))/4);
				t2->setPixel(x, y, 0, (prev2->getPixel(x*2, y*2, 0)+prev2->getPixel(x*2+1, y*2, 0)+prev2->getPixel(x*2, y*2+1, 0)+prev2->getPixel(x*2+1, y*2+1, 0))/4);
			}
		}
		
		p1.push_back(t1);
		p2.push_back(t2);
		
		char buffer1[30];
		sprintf(buffer1,"t1_%d.bmp",currSize);
		t1->writeBMP(buffer1);
		
		char buffer2[30];
		sprintf(buffer2,"t2_%d.bmp",currSize);
		t2->writeBMP(buffer2);
	}
	
	cout << "should be perfect!\n result is " << correspondence(i1, i2, 3, 3, 29, 29) << endl;
	cout << "no match!\n result is " << correspondence(i1, i1, 3, 3, 0, 0) << endl;
	
	Image* result = new Image(100,100);
	return result;
}

#define ij(X,Y) ((X)+SIZE*(Y))

// findPath takes a source and destination image, the paths of the layer one smaller, and a vector of the same size
void findPath(Image* src, Image* dst, vector<Path*>* smallerMap, vector<Path*>* newMap)
{
	int SIZE = src->getHeight();
	// Assume newMap is initialized for now to random points
	
	// Add all pixels to calculation queue
	queue<XY*> calcQueue;
	XY* tmpXY;
	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			tmpXY = new XY();
			tmpXY->x = i;
			tmpXY->y = j;
			calcQueue.push(tmpXY);
		}
	}
	
	// Create calculation storage data structure
	//   Each pixel has 4 possible end locations, so we need to store that
	vector<double> calcStore(SIZE*SIZE*4);
	
	while (true) {
		// Empty queue to calculate and find the best move
		XY bestSrc;
		XY bestDst;
		while (!calcQueue.empty())
		{
			XY* pt = calcQueue.front();
			calcQueue.pop();
			
			int x = pt->x;
			int y = pt->y;
			
			// First find the baseline for this pixel
			double baseline = energy(src, dst, x, y, (*newMap)[ij(x,y)]);
			
			// Next calculate what kind of changes we can achieve with each of the 4 choices
		}
		
		// If there is a good move, make the change using that move
		
		// Add affected pixels to the calculation queue
	}
}

double energy(Image* src, Image* dst, int x, int y, Path* p)
{
	return 0;
}

double correspondence(Image* A, Image* B, int p1_x, int p1_y, int p2_x, int p2_y)
{
	int pA_x = p1_x;
	int pA_y = p1_y;
	int pB_x = p2_x;
	int pB_y = p2_y;
	
	// function
	double gradientCost = pow(sqrt(pow( gradientx(A, pA_x, pA_y) - gradientx(B, pB_x, pB_y),2) + pow( gradienty(A, pA_x, pA_y) - gradienty(B, pB_x, pB_y),2)),2);
	double intensityCost = 0.5*pow(sqrt( pow(A->getPixel(pA_x, pA_y, 0) - B->getPixel(pB_x, pB_y, 0),2 )),2);
	double totalCost = sqrt( (gradientCost + intensityCost) / ( stdDev(A, pA_x, pA_y) * stdDev(B, pB_x, pB_y) ) );
	cout << gradientCost << endl;
	cout << intensityCost << endl;
	cout << totalCost << endl;
	return exp(totalCost)-1;
}

double gradientx(Image* img, int x, int y)
{
	return img->getPixel_(x-1, y, 0) * (-0.5) + img->getPixel_(x+1, y, 0) * 0.5;
}


double stdDev(Image* image, int x, int y)
{
	double stdDev = 0;
	double mean = 0;
	double numNeighbors = 0;
	double coord[] = {x, y-1, x-1, y, x+1, y, x, y+1};
	
	for (int i =0; i<8; i+=2) {
		if (coord[i] >= 0 and coord[i] < image->getWidth() and coord[i+1] >= 0 and coord[i+1] < image->getHeight()) {
			mean += image->getPixel(coord[i], coord[i+1], 0);
			numNeighbors += 1;
		}
		else {
			coord[i] = -1;
			coord[i+1] = -1;
		}
		
	}
	mean /= numNeighbors;
	for (int i=0; i<8; i+=2) {
		if (coord[i]>=0) {
			stdDev += pow( image->getPixel(coord[i], coord[i+1], 0)-mean,2);
		}
	}
	stdDev /= numNeighbors;
	stdDev = sqrt(stdDev); // got the STDDEV of neighbors
	return (image->getPixel(x, y, 0)-mean)/stdDev;
}	

double gradienty(Image* img, int x, int y)
{	
	return img->getPixel_(x, y-1, 0) * (-0.5) + img->getPixel_(x, y+1, 0) * 0.5;
}

/*
 * Inter/extrapolate function
 */
Image* interpolate (Image* i1, Image* i2, double alpha)
{
	Image* result = new Image(i1->getWidth(), i1->getHeight(), i1->getChannels(), i1->getBits());
	int chan, w, h;

	for (chan=0; chan < result->getChannels(); chan++)
	{
		for (w=0; w < result->getWidth(); w++)
		{
			for (h=0; h < result->getHeight(); h++)
			{
				result->setPixel_(w, h, chan, alpha*i1->getPixel(w, h, chan)+(1-alpha)*i2->getPixel(w, h, chan));
			}
		}
	}
	return result;
}

/*
 * Convolution function
 */
Image* convolute (Image* src, double* k, int ksize){
	Image* result = new Image(src->getWidth(), src->getHeight(), src->getChannels(), src->getBits());
	int chan, w, h, kw, kh;

	for (chan=0; chan < result->getChannels(); chan++)
	{
		for (w=0; w < result->getWidth(); w++)
		{
			for (h=0; h < result->getHeight(); h++)
			{
				double r = 0.0;
				for (kw=-(ksize/2); kw <= ksize/2; kw++)
				{
					for (kh=-(ksize/2); kh <= ksize/2; kh++)
					{
						r += src->getPixel_(w+kw, h+kh, chan) * KERNEL(kw+ksize/2, kh+ksize/2);
					}
				}
				result->setPixel_(w, h, chan, r);
			}
		}
	}
	return result;
}

/*
 * Resample function
 */
double resample (Image* src, double w, double h, int chan, int mode, int size, double sigma)
{
	//if (w < 0 || h < 0 || w>=src->getWidth() || h>=src->getHeight())
	//	return 0.0;
	if (mode==I_NEAREST) {
		int x = (int) (w+0.5);
		int y = (int) (h+0.5);
		return src->getPixel_(x, y, chan);
	} else if (mode==I_BILINEAR) {
		int w1 = (int) w;
		int w2 = (int) (w+1);
		int h1 = (int) h;
		int h2 = (int) (h+1);
		double result1 = (w-w1)*src->getPixel_(w1, h1, chan)+(w2-w)*src->getPixel_(w2, h1, chan);
		double result2 = (w-w1)*src->getPixel_(w1, h2, chan)+(w2-w)*src->getPixel_(w2, h2, chan);
		return (h-h1)*result1 + (h2-h)*result2;
	} else if (mode==I_GAUSSIAN) {
		double result=0.0;
		int x=((int)w)-size/2+1, y=((int)h)-size/2+1;
		double* k = new double[size*size];
		int i, j;
		double sum = 0.0;
		for (i = 0; i < size; i++)
		{
			for (j = 0; j < size; j++)
			{
				k[i+j*size] = pow(M_E, -( (sqr(w-x+i)+sqr(h-y+j)) / (2*sqr(sigma)) ));
				sum += k[i+j*size];
			}
		}
		for (i=0; i<size*size; i++)
		{
			//printf("%f\n", k[i]);
			k[i] /= sum;
		}
		int xx, yy;
		for (xx=0; xx<size; xx++)
		{
			for (yy=0; yy<size; yy++)
			{
				result += k[xx+yy*size]*src->getPixel_(xx+x, yy+y, chan);
			}
		}
		return result;
	}
	return 0.0;
}

/*
 * convolve with a box filter
 */
Image* ip_blur_box (Image* src, int size)
{
	double* k = new double[size*size];
	int i;
	for (i = 0; i < size*size; i++)
		k[i] = 1.0/(size*size);
	Image* result = convolute(src, k, size);
	delete[] k;
	return result;
}


/*
 * convolve with a gaussian filter
 */
Image* ip_blur_gaussian (Image* src, int size, double sigma)
{
	double* k = new double[size*size];
	int i, j;
	double sum = 0.0;
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			k[i+j*size] = pow(M_E, -( (sqr(i-size/2)+sqr(j-size/2)) / (2*sqr(sigma)) ));
			sum += k[i+j*size];
		}
	}
	for (i=0; i<size*size; i++)
		k[i] /= sum;
	Image* result = convolute(src, k, size);
	delete[] k;
	return result;
}


/*
 * convolve with a triangle filter
 */
Image* ip_blur_triangle (Image* src, int size)
{
  return NULL;
}


/*
 * interpolate with a black image
 */
Image* ip_brighten (Image* src, double alpha)
{
	Image* black = new Image(src->getWidth(), src->getHeight(), src->getChannels());
	int chan, w, h;

	// Make black
	for (chan=0; chan < src->getChannels(); chan++)
	{
		for (w=0; w < src->getWidth(); w++)
		{
			for (h=0; h < src->getHeight(); h++)
			{
				black->setPixel(w, h, chan, 0.0);
			}
		}
	}
	return interpolate(src, black, alpha);
}


/*
 * interpolate with the average intensity of the src image
 */
Image* ip_contrast (Image* src, double alpha)
{
	Image* grey = new Image(src->getWidth(), src->getHeight(), src->getChannels());
	int chan, w, h;
	double avg = 0.0, final[3];
	for (chan=0; chan < src->getChannels(); chan++)
	{
		avg = 0.0;
		for (w=0; w < src->getWidth(); w++)
		{
			for (h=0; h < src->getHeight(); h++)
			{
				avg += src->getPixel(w, h, chan);
			}
		}
		final[chan] = avg / (src->getWidth()*src->getHeight());
	}
	double cr = .2126,  cg = .7152,  cb = .0722;
	avg = final[RED]*cr+final[GREEN]*cg+final[BLUE]*cb;

	for (chan=0; chan < src->getChannels(); chan++)
	{
		for (w=0; w < src->getWidth(); w++)
		{
			for (h=0; h < src->getHeight(); h++)
			{
				grey->setPixel(w, h, chan, avg);
			}
		}
	}

	return interpolate(src, grey, alpha);
}


/*
 * use a mask image for a per-pixel alpha value to perform
 * interpolation with a second image
 */
Image* ip_composite (Image* src, const char* imageName, 
                     const char* maskName)
{
	Image* result = new Image(src->getWidth(), src->getHeight(), src->getChannels());
	Image* src2 = new Image(imageName);
	Image* mask = new Image(maskName);
	int chan, w, h;

	for (chan=0; chan < src->getChannels(); chan++)
	{
		for (w=0; w < src->getWidth(); w++)
		{
			for (h=0; h < src->getHeight(); h++)
			{
				result->setPixel(w, h, chan, src->getPixel(w, h, chan)*mask->getPixel(w, h, chan)+src2->getPixel(w, h, chan)*(1-mask->getPixel(w, h, chan)));
			}
		}
	}
	return result;
}


/*
 * cut away all but a subset of the image
 */
Image* ip_crop (Image* src, int x0, int y0, int x1, int y1)
{
	int x = (x0-x1>0) ? x0-x1 : x1-x0;
	int y = (y0-y1>0) ? y0-y1 : y1-y0;
	int xx = min(x0, x1);
	int yy = min(y0, y1);
	Image* result = new Image(x, y, src->getChannels());
	int w, h, chan;
	for (chan=0; chan < src->getChannels(); chan++)
	{
		for (w=0; w < x; w++)
		{
			for (h=0; h < y; h++)
			{
				result->setPixel(w, h, chan, src->getPixel(xx+w, yy+h, chan));
			}
		}
	}
	return result;
}


/*
 * convolve with an edge detection kernel
 */
Image* ip_edge_detect (Image* src)
{
	double* k = new double[9];
	k[0] = -1;
	k[1] = -1;
	k[2] = -1;
	k[3] = -1;
	k[4] = 8;
	k[5] = -1;
	k[6] = -1;
	k[7] = -1;
	k[8] = -1;

	Image* result = convolute(src, k, 3);
	delete[] k;
	return result;
}


/*
 * create a new image with color values from one channel of src
 */
Image* ip_extract (Image* src, int channel)
{
	Image* result = new Image(src->getWidth(), src->getHeight(), src->getChannels());

	int w, h;
	for (w=0; w < src->getWidth(); w++)
	{
		for (h=0; h < src->getHeight(); h++)
		{
			result->setPixel(w, h, channel, src->getPixel(w, h, channel));
		}
	}
	return result;
}


/*
 * create a new image with values equal to the psychosomatic intensities
 * of the source image
 */
Image* ip_grey (Image* src)
{
	Image* result = new Image(src->getWidth(), src->getHeight(), 1);
	
	double cr = .2126,  cg = .7152,  cb = .0722;

	int w, h;
	for (w=0; w < src->getWidth(); w++)
	{
		for (h=0; h < src->getHeight(); h++)
		{
			result->setPixel(w, h, 0, cr*src->getPixel(w,h,RED)+cg*src->getPixel(w,h,GREEN)+cb*src->getPixel(w,h,BLUE));
		}
	}
	return result;
}

Image* ip_grey_chan (Image* src)
{
	Image* result = new Image(src->getWidth(), src->getHeight(), src->getChannels());
	
	double cr = .2126,  cg = .7152,  cb = .0722;

	int w, h, chan;
	for (chan = 0; chan < src->getChannels(); chan++){
	for (w=0; w < src->getWidth(); w++)
	{
		for (h=0; h < src->getHeight(); h++)
		{
			result->setPixel(w, h, chan, cr*src->getPixel(w,h,RED)+cg*src->getPixel(w,h,GREEN)+cb*src->getPixel(w,h,BLUE));
		}
	}
	}
	return result;
}


/*
 * subtract the image from a white image
 */
Image* ip_invert (Image* src)
{
	Image* grey = new Image(src->getWidth(), src->getHeight(), src->getChannels());
	int chan, w, h;

	for (chan=0; chan < src->getChannels(); chan++)
	{
		for (w=0; w < src->getWidth(); w++)
		{
			for (h=0; h < src->getHeight(); h++)
			{
				grey->setPixel(w, h, chan, 0.5);
			}
		}
	}
	return interpolate(src, grey, -1);
}


/*
 * interpolate with random noise
 */
Image* ip_noisify (Image* src, double alpha)
{
	Image* noise = new Image(src->getWidth(), src->getHeight(), src->getChannels());
	int chan, w, h;

	// Make black
	for (chan=0; chan < src->getChannels(); chan++)
	{
		for (w=0; w < src->getWidth(); w++)
		{
			for (h=0; h < src->getHeight(); h++)
			{
				noise->setPixel(w, h, chan, ((float)rand())/RAND_MAX);
			}
		}
	}
	return interpolate(src, noise, alpha);
}


/*
 * round each pixel to the nearest value in the new number of bits
 */
Image* ip_quantize_simple (Image* src, int bitsPerChannel)
{
	Image* result = new Image(src->getWidth(), src->getHeight(), src->getChannels(), bitsPerChannel);
	int chan, w, h;

	for (chan=0; chan < src->getChannels(); chan++)
	{
		for (w=0; w < src->getWidth(); w++)
		{
			for (h=0; h < src->getHeight(); h++)
			{
				result->setPixel_(w, h, chan, src->getPixel(w, h, chan));
			}
		}
	}
	return result;
}


/*
 * dither each pixel to the nearest value in the new number of bits
 * using a static 4x4 matrix
 */
Image* ip_quantize_ordered (Image* src, int bitsPerChannel)
{
	if (bitsPerChannel == 8)
		return NULL;
	int matrix[4][4] = {{15,7,13,5},
						{3,11,1,9},
						{12,4,14,6},
						{0,8,2,10}};

	Image* result = new Image(src->getWidth(), src->getHeight(), src->getChannels(), src->getBits());
	int chan, w, h;
	for (chan=0; chan < src->getChannels(); chan++)
	{
		for (w=0; w < src->getWidth(); w++)
		{
			for (h=0; h < src->getHeight(); h++)
			{
				result->setPixel(w, h, chan, clamp( src->getPixel(w, h, chan)+(0.5-(2*matrix[w%4][h%4]+1)/(32.0)), 0,1 ));
			}
		}
	}
	return ip_quantize_simple(result, bitsPerChannel);
}


/*
 * dither each pixel to the nearest value in the new number of bits
 * using error diffusion
 */
Image* ip_quantize_fs (Image* src, int bitsPerChannel)
{
	Image* result = new Image(src->getWidth(), src->getHeight(), src->getChannels(), bitsPerChannel);
	int chan, w, h;
	double a=7.0/16, b=3.0/16, c=5.0/16, d=1.0/16, val, diff;

	for (chan=0; chan < src->getChannels(); chan++)
	{
		double* vals = new double[src->getWidth()*src->getHeight()];

		for (w = 0; w < src->getWidth()*src->getHeight(); w++)
			vals[w] = src->getPixel(w%src->getWidth(), w/src->getWidth(), chan);

		for (h=0; h < src->getHeight(); h++)
		{
			for (w=0; w < src->getWidth(); w++)
			{
				result->setPixel_(w, h, chan, vals[w+h*src->getWidth()]);
				val = result->getPixel(w, h, chan);
				diff = vals[w+h*src->getWidth()]-val;
				if (w+1 + h*src->getWidth() < src->getWidth()*src->getHeight()) vals[w+1 + h*src->getWidth()] += (diff*a);
				if (h+1 < src->getHeight()){
					if (w-1 >= 0) vals[w-1 + (h+1)*src->getWidth()] += (diff*b);
					vals[w + (h+1)*src->getWidth()] += (diff*c);
					if (w+1 < src->getWidth()) vals[w+1 + (h+1)*src->getWidth()] += (diff*d);
				}
			}
		}
	}
	return result;
}

bool close (double d, int i)
{
	return ((int)floor(d-0.2) == i && (int)floor(d+0.8) == i);
}

/*
 * rotate image using one of three sampling techniques
 */
Image* ip_rotate (Image* src, double theta, int x, int y, int mode, 
                  int size, double sigma)
{
	Image* result = new Image(src->getWidth(), src->getHeight(), src->getChannels());
	int w, h, chan;
	for (w=0; w < src->getWidth(); w++)
	{
		for (h=0; h < src->getHeight(); h++)
		{
			// Shift coord
			int dx = w-x, dy=h-y;
			// Find new location
			double thetar = deg2rad(-theta);
			double dxp = dx*cos(thetar)+dy*sin(thetar);
			double dyp = dy*cos(thetar)-dx*sin(thetar);
			// Restore coord
			double newx = dxp+x, newy = dyp+y;
			for (chan = 0; chan < src->getChannels(); chan++)
			{
				if (close(newx, w) && close(newy, h))
					result->setPixel_(w, h, chan, src->getPixel(w, h, chan));
				else
					result->setPixel_(w, h, chan, resample(src, newx, newy, chan, mode, size, sigma));
			}
		}
	}
	return result;
}


/*
 * interpolate with the greyscale version of the image
 */
Image* ip_saturate (Image* src, double alpha)
{
	Image* grey = ip_grey_chan(src);
	return interpolate(src, grey, alpha);
}


/*
 * scale image using one of three sampling techniques
 */
Image* ip_scale (Image* src, double xFac, double yFac, int mode, 
                 int size, double sigma)
{
	Image* result = new Image(src->getWidth(), src->getHeight(), src->getChannels());
	int w, h, chan;
	for (chan = 0; chan < src->getChannels(); chan++)
	{
		for (w=0; w < src->getWidth(); w++)
		{
			for (h=0; h < src->getHeight(); h++)
			{
				result->setPixel_(w, h, chan, resample(src, w/xFac, h/yFac, chan, mode, size, sigma));
			}
		}
	}
	return result;
}


/*
 * create a new one bit/channel image with the intensity 0 or 1
 * depending on whether the input value is above or below the 
 * threshold
 */
Image* ip_threshold (Image* src, double cutoff)
{
	Image* result = new Image(src->getWidth(), src->getHeight(), src->getChannels(), 1);
	int chan, w, h;
	for (chan=0; chan<src->getChannels(); chan++)
	{
		for (w=0; w < src->getWidth(); w++)
		{
			for (h=0; h < src->getHeight(); h++)
			{
				result->setPixel(w, h, chan, int(src->getPixel(w, h, chan) > cutoff));
			}
		}
	}
	return result;
}

/*
 * Spherical
 */
Image* ip_warp (Image* src, int warpFactor)
{
	Image* result = new Image(src->getWidth(), src->getHeight(), src->getChannels());
	int w, h, chan, x=src->getWidth()/2, y=src->getHeight()/2;
	for (w=0; w < src->getWidth(); w++)
	{
		for (h=0; h < src->getHeight(); h++)
		{
			int dx = w-x, dy=h-y;
			double theta = atan2((double)dy, dx);
			double r = sqrt((double)dx*dx+dy*dy);
			double rr = sqrt(x*cos(theta)*x*cos(theta)+y*sin(theta)*y*sin(theta));
			double newx = asin(r/rr)*2/M_PI*r*cos(theta)+x;
			double newy = asin(r/rr)*2/M_PI*r*sin(theta)+y;
			for (chan = 0; chan < src->getChannels(); chan++)
			{
				result->setPixel(w, h, chan, resample(src, newx, newy, chan, I_NEAREST, 4, 2));
			}
		}
	}
	return result;
}
