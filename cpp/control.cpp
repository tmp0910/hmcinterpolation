#include "control.h"
#include "ip.h"
#include "main.h"
#include <stdlib.h>


/*
 * IMPORTANT - DO NOT CHANGE THIS FILE - IMPORTANT
 */


enum {
	M_QUIT = 0,
	M_HELP = 1,
	
	M_FILE_OPEN =2,
	M_FILE_SAVE =3,
	M_FILE_INFO =4,
	M_FILE_REVERT =5,
	
	M_PROCESS_BLUR_BOX =6,
	M_PROCESS_BLUR_GAUSSIAN =7,
	M_PROCESS_BLUR_TRIANGLE=8,
	
	M_PROCESS_BRIGHTEN=9,
	M_PROCESS_CONTRAST=10,
	M_PROCESS_COMPOSITE=11,
	M_PROCESS_CROP=12,
	M_PROCESS_EDGE_DETECT=13,
	M_PROCESS_EXTRACT=14,
	M_PROCESS_GREY=15,
	M_PROCESS_INVERT=16,
	M_PROCESS_NOISIFY=17,
	
	M_PROCESS_QUANTIZE_SIMPLE=18,
	M_PROCESS_QUANTIZE_ORDERED=19,
	M_PROCESS_QUANTIZE_FLOYD_STEINBERG=20,
	
	M_PROCESS_ROTATE=21,
	M_PROCESS_SATURATE=22,
	M_PROCESS_SCALE=23,
	M_PROCESS_THRESHOLD=24,
	M_PROCESS_WARP=25,
	M_PROCESS_INTERPOLATE=26,
	
	M_LAST_ENUM
} MENU_ITEMS;


int make_menu ()
{
	int file = glutCreateMenu(menu_func);
	glutAddMenuEntry( "Open...",		M_FILE_OPEN);
	glutAddMenuEntry( "Save...",		M_FILE_SAVE);
	glutAddMenuEntry( "Get Image Info",		M_FILE_INFO);
	glutAddMenuEntry( "Revert",		M_FILE_REVERT);
	
//	int blur = glutCreateMenu(menu_func);
//	glutAddMenuEntry( "Box...",		M_PROCESS_BLUR_BOX);
//	glutAddMenuEntry( "Gaussian...",	M_PROCESS_BLUR_GAUSSIAN);
//	glutAddMenuEntry( "Triangle...",	M_PROCESS_BLUR_TRIANGLE);
//	
//	int quantize = glutCreateMenu(menu_func);
//	glutAddMenuEntry( "Simple...",	M_PROCESS_QUANTIZE_SIMPLE);
//	glutAddMenuEntry( "Ordered...",	M_PROCESS_QUANTIZE_ORDERED);
//	glutAddMenuEntry( "Floyd-Steinberg...", M_PROCESS_QUANTIZE_FLOYD_STEINBERG);
	
	int process = glutCreateMenu(menu_func);
	glutAddMenuEntry( "INTERPOLATE...",	M_PROCESS_INTERPOLATE);
//	glutAddSubMenu(   "Blur",		blur);
//	glutAddMenuEntry( "Brighten...",	M_PROCESS_BRIGHTEN);
//	glutAddMenuEntry( "Contrast...",	M_PROCESS_CONTRAST);
//	glutAddMenuEntry( "Composite...",	M_PROCESS_COMPOSITE);
//	glutAddMenuEntry( "Crop...",		M_PROCESS_CROP);
//	glutAddMenuEntry( "Edge Detect",	M_PROCESS_EDGE_DETECT);
//	glutAddMenuEntry( "Extract...",	M_PROCESS_EXTRACT);
//	glutAddMenuEntry( "Grey",		M_PROCESS_GREY);
//	glutAddMenuEntry( "Invert",		M_PROCESS_INVERT);
//	glutAddMenuEntry( "Noisify...",	M_PROCESS_NOISIFY);
//	glutAddSubMenu(   "Quantize",		quantize);
//	glutAddMenuEntry( "Rotate...",	M_PROCESS_ROTATE);
//	glutAddMenuEntry( "Saturate...",	M_PROCESS_SATURATE);
//	glutAddMenuEntry( "Scale...",		M_PROCESS_SCALE);
//	glutAddMenuEntry( "Threshold...",	M_PROCESS_THRESHOLD);
//	glutAddMenuEntry( "Warp",		M_PROCESS_WARP);
	
	int main = glutCreateMenu(menu_func);
	glutAddSubMenu(   "File",		file);
	glutAddSubMenu(   "Process",		process);
	glutAddMenuEntry( "Help",		M_HELP);
	glutAddMenuEntry( "Quit",		M_QUIT);
	
	glutAttachMenu(GLUT_RIGHT_BUTTON);
	
	return main;
}


static inline void checkStream (const istream& in)
{
	if (in.fail())
	{
		cerr << "Fatal error: stream failed!" << endl;
		exit(-1);
	}
}


void menu_func (int value)
{
	// variables used in the switch statement
	char filename1[MAX_LINE], filename2[MAX_LINE];
	double alpha, sigma, theta, xFactor, yFactor;
	alpha=theta=xFactor=yFactor=0;
	sigma=1;
	int filtersize, x0, x1, y0, y1, channel, bitsPerChannel, mode;
	filtersize=x0=x1=y0=y1=channel=bitsPerChannel=mode=0;
	Image* resultImage = NULL;
	
	switch (value)
	{
		case M_QUIT:
			exit(0);
			break;
			;;
			
			
		case M_HELP:
			menu_help();
			break;
			;;
			
			
		case M_FILE_OPEN:   // enum #2
			if (!quietMode)
				cerr << "Open file (string - no spaces try /Users/lokisnake/Downloads/AWESOME.bmp ) : ";
			cin  >> filename1;
			checkStream(cin);
			image_load(filename1);
			break;
			
			
		case M_FILE_SAVE:   // enum #3
			if (!quietMode)
				cerr << "Save as (string - no spaces) : ";
			cin  >> filename1;
			checkStream(cin);
			image_save(filename1);
			break;
			
			
		case M_FILE_INFO:
			image_print_info();
			break;
			
			
		case M_FILE_REVERT:  // enum #5
			image_revert();
			break;
			
			

		#pragma mark INTERPOLATE
		case M_PROCESS_INTERPOLATE:
//			if (!currentImage)
//			{
//				cerr << "No Image!" << endl;
//				return;
//			}
			
			if (!quietMode)
				cerr << "First image file (string - no spaces) : ";
			cin  >> filename1;
			if (!quietMode)
				cerr << "Composite with image file (string - no spaces) : ";
			cin  >> filename2;
			if (!quietMode)
				cerr << "Interpolation : ";
			double inter;
			cin  >> inter;
			
			checkStream(cin);
			
			resultImage = ip_interpolate(filename1, filename2, inter);
			break;
			
			
		default:
			break;
			;;
	}
	
	if (resultImage != NULL)
	{
		delete currentImage;
		currentImage = resultImage;
		
		if (currentImage->getWidth()  != window_width    ||
			currentImage->getHeight() != window_height)
			reshape(currentImage->getWidth(), currentImage->getHeight());
		
		if (!quietMode)
			cerr << "done!" << endl;
		
		if (!textMode)
			glutPostRedisplay();
	}
}


void keyboard_func (unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'H':
		case 'h':
			menu_help();
			break;
			;;
			
		case 'Q':
		case 'q':
			exit(0);
			break;
			;;
	}
}


void menu_help ()
{
	cerr << endl
	<< "hmc cs155 image processor" << endl
	<< "please see the ip manual for usage and algorithm information" << endl
	<< "http://www.cs.hmc.edu/courses/2002/fall/cs155/proj1/doc/ip_manual.html"
	<< endl << endl;
}


#define MENUOP(num, tag)	cerr << " " << num << ") " << tag << endl;

void printMenu ()
{
    if (quietMode)
		return;
	
    MENUOP(M_QUIT, "quit")
    MENUOP(M_HELP, "help")
    cerr << endl
	<< "file" << endl;
    MENUOP(M_FILE_OPEN,			"open...")
    MENUOP(M_FILE_SAVE,			"save...")
    MENUOP(M_FILE_INFO,			"get info")
    MENUOP(M_FILE_REVERT,		"revert")
    cerr << endl
	<< "process" << endl;
    MENUOP(M_PROCESS_INTERPOLATE,		"INTERPOLATE...")
    MENUOP(M_PROCESS_BLUR_BOX,		"box blur...")
    MENUOP(M_PROCESS_BLUR_GAUSSIAN,	"gaussian blur...")
    MENUOP(M_PROCESS_BLUR_TRIANGLE,	"triangle blur...")
    MENUOP(M_PROCESS_BRIGHTEN,		"brighten...")
    MENUOP(M_PROCESS_CONTRAST,		"contrast...")
    MENUOP(M_PROCESS_COMPOSITE,		"composite...")
    MENUOP(M_PROCESS_CROP,		"crop...")
    MENUOP(M_PROCESS_EDGE_DETECT,	"edge detect")
    MENUOP(M_PROCESS_EXTRACT,		"extract...")
    MENUOP(M_PROCESS_GREY,		"grey")
    MENUOP(M_PROCESS_INVERT,		"invert")
    MENUOP(M_PROCESS_NOISIFY,		"noisify...")
    MENUOP(M_PROCESS_QUANTIZE_SIMPLE,	"simple quantize...")
    MENUOP(M_PROCESS_QUANTIZE_ORDERED,	"ordered quantize...")
    MENUOP(M_PROCESS_QUANTIZE_FLOYD_STEINBERG, "floyd-steinberg quantize...")
    MENUOP(M_PROCESS_ROTATE,		"rotate...")
    MENUOP(M_PROCESS_SATURATE,		"saturate...")
    MENUOP(M_PROCESS_SCALE,		"scale...")
    MENUOP(M_PROCESS_THRESHOLD,		"threshold...")
    MENUOP(M_PROCESS_WARP,		"warp")
}


void textMenuLoop ()
{
	char command[MAX_LINE];
	
	printMenu();
	
	while (true)
	{
		if (!quietMode)
			cerr << endl
			<< "selection > " << flush;
		cin  >> command;
		
		switch (command[0])
		{
			case '\n':
			case '\0':
				printMenu();
				break;
				
			case 'Q':
			case 'q':
				menu_func(M_QUIT);
				break;
				
			case 'H':
			case 'h':
				menu_func(M_HELP);
				break;
				
			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9':
				menu_func(atoi(command));
				break;
				
			default:
				printMenu();
				break;
		}
	}
}


void image_load (const char* filename)
{
	if (currentImage)
		delete currentImage;
	if (originalImage)
		delete originalImage;
	currentImage  = NULL;
	originalImage = NULL;
    
	originalImage = new Image();
	originalImage->read(filename);
	
	if (originalImage->good())
	{  
		currentImage = new Image(*originalImage);
		reshape(currentImage->getWidth(), currentImage->getHeight());
	}
	else
	{
		delete originalImage;  
		originalImage = NULL;
		cerr << "Couldn't load image " << filename << "!" << endl;
		return;
	}
	
	if (!textMode)
		glutPostRedisplay();
	
	if (!quietMode)
		cerr << "done!" << endl;
}  


void image_save (const char* filename)
{
	if (currentImage)
	{
		if (currentImage->write(filename) == 0)
		{
			//delete originalImage;
			//originalImage = new Image(*currentImage);
		}
	}  
	else if (originalImage)
	{
		originalImage->write(filename);
	}
	else
	{
		cerr << "No image!" << endl;
		return;
	}
	
	if (!quietMode)
		cerr << "done!" << endl;
}


void image_print_info ()
{  
	if (currentImage) {
		cerr << "width:    " << currentImage->getWidth() << endl
		<< "height:   " << currentImage->getHeight() << endl
		<< "channels: " << currentImage->getChannels() << endl
		<< "bits:     " << currentImage->getBits() << endl;
	}
	cerr << "done!" << endl;
}


void image_revert ()
{
	if (currentImage)
		delete currentImage;
	
	if (originalImage)
	{
		currentImage = new Image(*originalImage);
		
		if (window_width  != currentImage->getWidth() ||
			window_height != currentImage->getHeight())
			reshape(currentImage->getWidth(), currentImage->getHeight());
	}
	else
	{
		cerr << "No image!" << endl;
		return;
	}
	
	if (!textMode)
		glutPostRedisplay();
	
	if (!quietMode)
		cerr << "done!" << endl;
}  

