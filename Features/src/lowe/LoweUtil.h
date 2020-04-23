#ifndef _LOWE_UTIL_H_
#define _LOWE_UTIL_H_

/************************************************************************
  Copyright (c) 2003. David G. Lowe, University of British Columbia.
  This software is being made freely available for research purposes
  only (see file LICENSE.txt for conditions).  This notice must be
  retained on all copies.
*************************************************************************/

/* From the standard C libaray: */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*------------------------- Macros and constants  -------------------------*/

/* Following defines TRUE and FALSE if not previously defined. */
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/* Value of PI, rounded up, so orientations are always in range [0,PI]. */
#ifndef PI
#define PI 3.1415927
#endif

#define ABS(x)    (((x) > 0) ? (x) : (-(x)))
#define MAX(x,y)  (((x) > (y)) ? (x) : (y))
#define MIN(x,y)  (((x) < (y)) ? (x) : (y))

/* Given the name of a structure, NEW allocates space for it in the
   given pool (see util.c) and returns a pointer to the structure.
*/
#define NEW(s,pool) ((struct s *) MallocPool(sizeof(struct s),pool))
#define NEWOBJ(s,pool) ((s *) MallocPool(sizeof(s),pool))

/* Assign a unique number to each pool of storage needed for this application. 
*/
#define PERM_POOL  0     /* Permanent storage that is never released. */
#define IMAGE_POOL 1     /* Data used only for the current image. */
#define KEY_POOL   2     /* Data for set of keypoints. */

/* These constants specify the size of the index vectors that provide
   a descriptor for each keypoint.  The region around the keypoint is
   sampled at OriSize orientations with IndexSize by IndexSize bins.
   VecLength is the length of the resulting index vector.
*/

#define OriSize 8
#define IndexSize 4
#define VecLength (IndexSize * IndexSize * OriSize)

//extern "C" {
/*---------------------------- Structures -------------------------------*/

/* Data structure for a float image.
*/
typedef struct LoweImageSt 
{
	int rows, cols;          /* Dimensions of image. */
	float **pixels;          /* 2D array of image pixels. */
	struct LoweImageSt *next;    /* Pointer to next image in sequence. */
} *LoweImage;

/* Data structure for a keypoint.  Lists of keypoints are linked
 *    by the "next" field.
 *    */
typedef struct LoweKeypointSt 
{
	float row, col;             /* Subpixel location of keypoint. */
	float scale, ori;           /* Scale and orientation (range [-PI,PI]) */
	unsigned char *ivec;        /* Vector of descriptor values */
	struct LoweKeypointSt *next;    /* Pointer to next keypoint in list. */
} *LoweKeypoint;

/*------------------------------- Externals -----------------------------*/

//const int MagFactor = 3;
//const float GaussTruncate = 4.0f;
extern const float GaussTruncate;
extern const int MagFactor;

/*-------------------------- Function prototypes -------------------------*/
/* The are prototypes for the external functions that are shared
   between files.
*/

/* Only interface needed to key.c. */
LoweKeypoint GetLoweKeypoints(LoweImage image,
                              int firstOctave,int numOctaves,int numScales,
                              float edgeThresh, float peakThresh);

/* Following are from util.c */
void *MallocPool(int size, int pool);
void FreeStoragePool(int pool);
float **AllocMatrix(int rows, int cols, int pool);
LoweImage CreateLoweImage(int rows, int cols, int pool);
LoweImage CopyLoweImage(LoweImage image, int pool);
LoweImage DoubleSize(LoweImage image);
LoweImage HalfLoweImageSize(LoweImage image);
void SubtractLoweImage(LoweImage im1, LoweImage im2);
void GradOriLoweImages(LoweImage im, LoweImage grad, LoweImage ori);
void GaussianBlur(LoweImage image, float sigma);
void SolveLeastSquares(float *solution, int rows, int cols, float **jacobian,
		       float *errvec, float **sqarray);
void SolveLinearSystem(float *solution, float **sq, int size);
float DotProd(float *v1, float *v2, int len);

void FatalError(char *fmt, ...);
LoweImage ReadPGMFile(char *filename);
LoweImage ReadPGM(FILE *fp);
void WritePGM(FILE *fp, LoweImage image);
void DrawLine(LoweImage image, int r1, int c1, int r2, int c2);
LoweKeypoint ReadKeyFile(char *filename);
LoweKeypoint ReadKeys(FILE *fp);

void SkipComments(FILE *fp);


#endif
