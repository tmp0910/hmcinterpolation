#ifndef MAIN_H
#define MAIN_H


#include "common.h"
#include "image.h"


/*
 * IMPORTANT - DO NOT CHANGE THIS FILE - IMPORTANT
 */


extern int window_width;
extern int window_height;

extern Image* currentImage;
extern Image* originalImage;

extern bool quietMode;
extern bool textMode;


int  main (int argc, char** argv);
char* init (int argc, char** argv);
void usage ();
void display ();
void unreshape (int width, int height);
void reshape (int width, int height);



#endif // MAIN_H
