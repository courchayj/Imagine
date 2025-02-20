/************************************************************************
  Copyright (c) 2003. David G. Lowe, University of British Columbia.
  This software is being made available for research purposes only
  (see file LICENSE for conditions).  This notice must be retained on
  all copies or modifications of this software.
*************************************************************************/

/* util.c
   This file contains a range of utility routines to support keypoint
   detection.
*/

#include "LoweUtil.h"
#include <stdarg.h>


//extern "C" {

/*---------------------- Local function prototypes -----------------------*/

void ConvHorizontal(LoweImage image, float *kernel, int ksize);
void ConvVertical(LoweImage image, float *kernel, int ksize);
void ConvBuffer(float *buffer, float *kernel, int rsize, int ksize);
void ConvBufferFast(float *buffer, float *kernel, int rsize, int ksize);


/*-------------------- Pooled storage allocator ---------------------------*/

/* The following routines allow for the efficient allocation of
     storage in small chunks from a named pool.  Rather than requiring
     each structure to be freed individually, an entire pool of
     storage is freed at once.
   This method has two advantages over just using malloc() and free().
     First, it is far more efficient for allocating small objects, as
     there is no overhead for remembering all the information needed
     to free each object.  Second, the decision about how long to keep
     an object is made at the time of allocation by assigning it to a
     pool, and there is no need to track down all the objects to free
     them.  In practice, this leads to code with little chance of
     memory leaks.

   Example of how to use the pooled storage allocator:
     Each pool is given a name that is a small integer (in header file):
       #define IMAGE_POOL 2
     Following allocates memory of "size" bytes from pool, IMAGE_POOL:
       mem = MallocPool(size, IMAGE_POOL);
     Following releases all memory in pool IMAGE_POOL for reuse:
       FreeStoragePool(IMAGE_POOL);
*/

/* We maintain memory alignment to word boundaries by requiring that
   all allocations be in multiples of the machine wordsize.  WORDSIZE
   is the maximum size of the machine word in bytes (must be power of 2).
   BLOCKSIZE is the minimum number of bytes requested at a time from
   the system (should be multiple of WORDSIZE).
*/
#define WORDSIZE 8  
#define BLOCKSIZE 2048

/* Following defines the maximum number of different storage pools. */
#define POOLNUM 100

/* Pointers to base of current block for each storage pool (C automatically
   initializes static memory to NULL).
*/
static char *PoolBase[POOLNUM];

/* Pointers to base of current block "to delete" for each storage pool (C automatically
   initializes static memory to NULL).
*/
static char *PoolBase_del[POOLNUM];

/* Number of bytes left in current block for each storage pool (initialized
   to 0). */
static int PoolRemain[POOLNUM];


/* Returns a pointer to a piece of new memory of the given size in bytes
   allocated from a named pool. 
*/
void *MallocPool(int size, int pool)
{
    char *m, **prev;
    int bsize;

    /* Round size up to a multiple of wordsize.  The following expression
	only works for WORDSIZE that is a power of 2, by masking last bits of
	incremented size to zero.
	*/
	size = (size + WORDSIZE - 1) & ~(WORDSIZE - 1);

	/* Check whether new block must be allocated.  Note that first word of
	block is reserved for pointer to previous block. */
	if (size > PoolRemain[pool]) {
		bsize = (size + sizeof(char **) > BLOCKSIZE) ?
			size + sizeof(char **) : BLOCKSIZE;
		m = (char*) malloc(bsize);

		if (! m) {
			fprintf(stderr, "ERROR: Ran out of memory.\n");
			abort();
		}
		PoolRemain[pool] = bsize - sizeof(void *);
		/* Fill first word of new block with pointer to previous block. */
		prev = (char **) m;
		prev[0] = PoolBase[pool];
		PoolBase[pool] = m;
	}
	/* Allocate new storage from end of the block. */
	PoolRemain[pool] -= size;
	return (PoolBase[pool] + sizeof(char **) + PoolRemain[pool]);
}


/* Free all storage that was previously allocated with MallocPool from
   a particular named pool. 
*/
void FreeStoragePool(int pool)
{
    char *prev;
   // printf(" Free ");


    while (PoolBase[pool] != NULL) {
	prev = *((char **) PoolBase[pool]);  /* Get pointer to prev block. */
	free(PoolBase[pool]);
	PoolBase[pool] = prev;
    }
    PoolRemain[pool] = 0;
}




/*----------------- 2D matrix and image allocation ------------------------*/

/* Allocate memory for a 2D float matrix of size [row,col].  This returns
     a vector of pointers to the rows of the matrix, so that routines
     can operate on this without knowing the dimensions.
*/
float **AllocMatrix(int rows, int cols, int pool)
{
    int i;
    float **m, *v;

    m = (float **) MallocPool(rows * sizeof(float *), pool);
    v = (float *) MallocPool(rows * cols * sizeof(float), pool);
    for (i = 0; i < rows; i++) {
	m[i] = v;
	v += cols;
    }
    return (m);
}


/* Create a new image with uninitialized pixel values.
*/
LoweImage CreateLoweImage(int rows, int cols, int pool)
{
    LoweImage im;

    im = NEWOBJ(struct LoweImageSt, pool);
    im->rows = rows;
    im->cols = cols;
    im->pixels = AllocMatrix(rows, cols, pool);
    return im;
}




/* Return a new copy of the image.
*/
LoweImage CopyLoweImage(LoweImage image, int pool)
{
    int r, c;
    LoweImage newIma;

    newIma = CreateLoweImage(image->rows, image->cols, pool);

    for (r = 0; r < image->rows; r++)
      for (c = 0; c < image->cols; c++)
	newIma->pixels[r][c] = image->pixels[r][c];
    return newIma;
}


/*----------------------- LoweImage utility routines ----------------------*/

/* Double image size. Use linear interpolation between closest pixels.
   Size is two rows and columns short of double to simplify interpolation.
*/
LoweImage DoubleSize(LoweImage image)
{
   int rows, cols, nrows, ncols, r, c, r2, c2;
   float **im, **newIma;
   LoweImage newimage;
   
   rows = image->rows;
   cols = image->cols;
   nrows = 2 * rows - 2;
   ncols = 2 * cols - 2;
   newimage = CreateLoweImage(nrows, ncols, IMAGE_POOL);
   im = image->pixels;
   newIma = newimage->pixels;
   
   for (r = 0; r < rows - 1; r++)
      for (c = 0; c < cols - 1; c++) {
         r2 = 2 * r;
         c2 = 2 * c;
         newIma[r2][c2] = im[r][c];
         newIma[r2+1][c2] = (float)(0.5 * (im[r][c] + im[r+1][c]));
         newIma[r2][c2+1] = (float)(0.5 * (im[r][c] + im[r][c+1]));
         newIma[r2+1][c2+1] = (float)(0.25 * (im[r][c] + im[r+1][c] + im[r][c+1] + im[r+1][c+1]));
      }
   return newimage;
}


/* Reduce the size of the image by half by selecting alternate pixels on
   every row and column.  We assume image has already been blurred
   enough to avoid aliasing.
*/
LoweImage HalfLoweImageSize(LoweImage image)
{
   int rows, cols, nrows, ncols, r, c, ri, ci;
   float **im, **newIma;
   LoweImage newimage;
   
   rows = image->rows;
   cols = image->cols;
   nrows = rows / 2;
   ncols = cols / 2;
   newimage = CreateLoweImage(nrows, ncols, IMAGE_POOL);
   im = image->pixels;
   newIma = newimage->pixels;
   
   for (r = 0, ri = 0; r < nrows; r++, ri += 2)
      for (c = 0, ci = 0; c < ncols; c++, ci += 2)
	newIma[r][c] = im[ri][ci];

   return newimage;
}


/* Subtract image im2 from im1 and leave result in im1.
*/
void SubtractLoweImage(LoweImage im1, LoweImage im2)
{
   float **pix1, **pix2;
   int r, c;
   
   pix1 = im1->pixels;
   pix2 = im2->pixels;
   
   for (r = 0; r < im1->rows; r++)
      for (c = 0; c < im1->cols; c++)
	pix1[r][c] -= pix2[r][c];
}


/* Given a smoothed image, im, return image gradient and orientation
     at each pixel in grad and ori.  Note that since we will eventually
     be normalizing gradients, we only need to compute their relative
     magnitude within a scale image (so no need to worry about change in
     pixel sampling relative to sigma or other scaling issues).
*/
void GradOriLoweImages(LoweImage im, LoweImage grad, LoweImage ori)
{
    float xgrad, ygrad, **pix;
    int rows, cols, r, c;
   
    rows = im->rows;
    cols = im->cols;
    pix = im->pixels;

    for (r = 0; r < rows; r++)
      for (c = 0; c < cols; c++) {
        if (c == 0)
          xgrad = 2.0f * (pix[r][c+1] - pix[r][c]);
        else if (c == cols-1)
          xgrad = 2.0f * (pix[r][c] - pix[r][c-1]);
        else
          xgrad = pix[r][c+1] - pix[r][c-1];
        if (r == 0)
          ygrad = 2.0f * (pix[r][c] - pix[r+1][c]);
        else if (r == rows-1)
          ygrad = 2.0f * (pix[r-1][c] - pix[r][c]);
        else
          ygrad = pix[r-1][c] - pix[r+1][c];
	grad->pixels[r][c] = (float)sqrt(xgrad * xgrad + ygrad * ygrad);
	ori->pixels[r][c] = (float)atan2 (ygrad, xgrad);
      }
}
   

/* --------------------------- Blur image --------------------------- */

/* Convolve image with a Gaussian of width sigma and store result back
     in image.   This routine creates the Gaussian kernel, and then applies
     it sequentially in horizontal and vertical directions.
*/
void GaussianBlur(LoweImage image, float sigma)
{
    float x, kernel[100], sum = 0.0;
    int ksize, i;

    /* The Gaussian kernel is truncated at GaussTruncate sigmas from
       center.  The kernel size should be odd.
    */
    ksize = (int)(2.0 * GaussTruncate * sigma + 1.0);
    ksize = MAX(3, ksize);    /* Kernel must be at least 3. */
    if (ksize % 2 == 0)       /* Make kernel size odd. */
      ksize++;
    assert(ksize < 100);

    /* Fill in kernel values. */
    for (i = 0; i <= ksize; i++) {
      x = (float)(i - ksize / 2);
      kernel[i] = (float)(exp(- x * x / (2.0 * sigma * sigma)));
      sum += kernel[i];
    }
    /* Normalize kernel values to sum to 1.0. */
    for (i = 0; i < ksize; i++)
      kernel[i] /= sum;

    ConvHorizontal(image, kernel, ksize);
    ConvVertical(image, kernel, ksize);
}


/* Convolve image with the 1-D kernel vector along image rows.  This
   is designed to be as efficient as possible.  Pixels outside the
   image are set to the value of the closest image pixel.
*/
void ConvHorizontal(LoweImage image, float *kernel, int ksize)
{
    int rows, cols, r, c, i, halfsize;
    float **pixels;//, buffer[40001];

    rows = image->rows;
    cols = image->cols;
    halfsize = ksize / 2;
    pixels = image->pixels;
    //assert(cols + ksize < 40001);
	float *buffer = new float[cols+ksize];

    for (r = 0; r < rows; r++) {
	/* Copy the row into buffer with pixels at ends replicated for
	   half the mask size.  This avoids need to check for ends
	   within inner loop. */
	for (i = 0; i < halfsize; i++)
	    buffer[i] = pixels[r][0];
	for (i = 0; i < cols; i++)
	    buffer[halfsize + i] = pixels[r][i];
	for (i = 0; i < halfsize; i++)
	    buffer[halfsize + cols + i] = pixels[r][cols - 1];

	ConvBufferFast(buffer, kernel, cols, ksize);
	for (c = 0; c < cols; c++)
	  pixels[r][c] = buffer[c];
    }

	delete [] buffer;
}


/* Same as ConvHorizontal, but apply to vertical columns of image.
*/
void ConvVertical(LoweImage image, float *kernel, int ksize)
{
    int rows, cols, r, c, i, halfsize;
    float **pixels;//, buffer[40001];

    rows = image->rows;
    cols = image->cols;
    halfsize = ksize / 2;
    pixels = image->pixels;
    //assert(rows + ksize < 40001);
	float* buffer = new float[rows+ksize];

    for (c = 0; c < cols; c++) {
	for (i = 0; i < halfsize; i++)
	    buffer[i] = pixels[0][c];
	for (i = 0; i < rows; i++)
	    buffer[halfsize + i] = pixels[i][c];
	for (i = 0; i < halfsize; i++)
	    buffer[halfsize + rows + i] = pixels[rows - 1][c];

	ConvBufferFast(buffer, kernel, rows, ksize);
	for (r = 0; r < rows; r++)
	  pixels[r][c] = buffer[r];
    }
	delete [] buffer;
}


/* Perform convolution of the kernel with the buffer, returning the
   result at the beginning of the buffer.  rsize is the size
   of the result, which is buffer size - ksize + 1.
*/
void ConvBuffer(float *buffer, float *kernel, int rsize, int ksize)
{
    int i, j;
    float sum, *bp, *kp;

    for (i = 0; i < rsize; i++) {
      sum = 0.0;
      bp = &buffer[i];
      kp = &kernel[0];

      /* Make this inner loop super-efficient. */
      for (j = 0; j < ksize; j++)
	sum += *bp++ * *kp++;
	
      buffer[i] = sum;
    }
}


/* Same as ConvBuffer, but implemented with loop unrolling for increased
   speed.  This is the most time intensive routine in keypoint detection,
   so deserves careful attention to efficiency.  Loop unrolling simply
   sums 5 multiplications at a time to allow the compiler to schedule
   operations better and avoid loop overhead.  This almost triples
   speed of previous version on a Pentium with gcc.
*/
void ConvBufferFast(float *buffer, float *kernel, int rsize, int ksize)
{
    int i;
    float sum, *bp, *kp, *endkp;

    for (i = 0; i < rsize; i++) {
      sum = 0.0;
      bp = &buffer[i];
      kp = &kernel[0];
      endkp = &kernel[ksize];

      /* Loop unrolling: do 5 multiplications at a time. */
      while (kp + 4 < endkp) {
	sum += bp[0] * kp[0] +  bp[1] * kp[1] + bp[2] * kp[2] +
	       bp[3] * kp[3] +  bp[4] * kp[4];
	bp += 5;
	kp += 5;
      }
      /* Do 2 multiplications at a time on remaining items. */
      while (kp + 1 < endkp) {
	sum += bp[0] * kp[0] +  bp[1] * kp[1];
	bp += 2;
	kp += 2;
      }
      /* Finish last one if needed. */
      if (kp < endkp)
	sum += *bp * *kp;
	
      buffer[i] = sum;
    }
}


/*--------------------- Least-squares solutions ---------------------------*/

/* Give a least-squares solution to the system of linear equations given in
   the jacobian and errvec arrays.  Return result in solution.
   This uses the method of solving the normal equations.
*/
void SolveLeastSquares(float *solution, int rows, int cols, float **jacobian,
		       float *errvec, float **sqarray)
{
    int r, c, i;
    float sum;

    assert(rows >= cols);
    /* Multiply Jacobian transpose by Jacobian, and put result in sqarray. */
    for (r = 0; r < cols; r++)
	for (c = 0; c < cols; c++) {
	    sum = 0.0;
	    for (i = 0; i < rows; i++)
		sum += jacobian[i][r] * jacobian[i][c];
	    sqarray[r][c] = sum;
	}
    /* Multiply transpose of Jacobian by errvec, and put result in solution. */
    for (c = 0; c < cols; c++) {
	sum = 0.0;
	for (i = 0; i < rows; i++)
	    sum += jacobian[i][c] * errvec[i];
	solution[c] = sum;
    }
    /* Now, solve square system of equations. */
    SolveLinearSystem(solution, sqarray, cols);
}
  

/* Solve the square system of linear equations, Ax=b, where A is given
   in matrix "sq" and b in the vector "solution".  Result is given in
   solution.  Uses Gaussian elimination with pivoting.
*/
void SolveLinearSystem(float *solution, float **sq, int size)
{
    int row, col, c, pivot = 0, i;
    float maxc, coef, temp, mult, val;

    /* Triangularize the matrix. */
    for (col = 0; col < size - 1; col++) {
	/* Pivot row with largest coefficient to top. */
	maxc = -1.0;
	for (row = col; row < size; row++) {
	    coef = sq[row][col];
	    coef = (coef < 0.0 ? - coef : coef);
	    if (coef > maxc) {
		maxc = coef;
		pivot = row;
	    }
	}
	if (pivot != col) {
	    /* Exchange "pivot" with "col" row (this is no less efficient
	       than having to perform all array accesses indirectly). */
	    for (i = 0; i < size; i++) {
		temp = sq[pivot][i];
		sq[pivot][i] = sq[col][i];
		sq[col][i] = temp;
	    }
	    temp = solution[pivot];
	    solution[pivot] = solution[col];
	    solution[col] = temp;
	}
	/* Do reduction for this column. */
	for (row = col + 1; row < size; row++) {
	    mult = sq[row][col] / sq[col][col];
	    for (c = col; c < size; c++)	/* Could start with c=col+1. */
		sq[row][c] -= mult * sq[col][c];
	    solution[row] -= mult * solution[col];
	}
    }

    /* Do back substitution.  Pivoting does not affect solution order. */
    for (row = size - 1; row >= 0; row--) {
	val = solution[row];
	for (col = size - 1; col > row; col--)
	    val -= solution[col] * sq[row][col];
	solution[row] = val / sq[row][row];
    }
}


/* Return dot product of two vectors with given length.
*/
float DotProd(float *v1, float *v2, int len)
{
    int i;
    float sum = 0.0;

    for (i = 0; i < len; i++)
      sum += v1[i] * v2[i];
    return sum;
}





/*------------------------ Error reporting ----------------------------*/

/* This function prints an error message and exits.  It takes a variable
 *    number of arguments that function just like those in printf.
 *    */
void FatalError(char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
    	fprintf(stderr, "Error: ");
	vfprintf(stderr, fmt, args);
    	fprintf(stderr,"\n");
	va_end(args);
	exit(1);
}



/*----------------- Read and write PGM files ------------------------*/


/* This reads a PGM file from a given filename and returns the image.
 * */
LoweImage ReadPGMFile(char *filename)
{
	FILE *file;

        /* The "b" option is for binary input, which is needed if this is
	 *        compiled under Windows.  It has no effect in Linux.
	 *            */
	file = fopen (filename, "rb");
	if (! file)
	    	FatalError("Could not open file: %s", filename);

	return ReadPGM(file);
}


/* Read a PGM file from the given file pointer and return it as a
 * float LoweImage structure with pixels in the range [0,1].  If the file
 * contains more than one image, then the images will be returned
 * linked by the "next" field of the LoweImage data structure.  
 * See "man pgm" for details on PGM file format.  This handles only
 * the usual 8-bit "raw" PGM format.  Use xv or the PNM tools (such as
 * pnmdepth) to convert from other formats.
 **/
LoweImage ReadPGM(FILE *fp)
{
	int char1, char2, width, height, max, c1, c2, c3, r;
	LoweImage image, nextimage;
	char1 = fgetc(fp);
	char2 = fgetc(fp);
	SkipComments(fp);
	c1 = fscanf(fp, "%d", &width);
	SkipComments(fp);
	c2 = fscanf(fp, "%d", &height);
	SkipComments(fp);
	c3 = fscanf(fp, "%d", &max);
	
	if (char1 != 'P' || char2 != '5' || c1 != 1 || c2 != 1 || c3 != 1 ||
		max > 255)
		FatalError("Input is not a standard raw 8-bit PGM file.\n"
				   "Use xv or pnmdepth to convert file to 8-bit PGM format.\n");
		  
	fgetc(fp);  /* Discard exactly one byte after header. */
	
	/* Create floating point image with pixels in range [0,1]. */
	image = CreateLoweImage(height, width,IMAGE_POOL);
	for (r = 0; r < height; r++)
		
		/* Check if there is another image in this file, as the latest PGM
		*      standard allows for multiple images. */
		SkipComments(fp);
	if (getc(fp) == 'P') {
		ungetc('P', fp);
		nextimage = ReadPGM(fp);
		image->next = nextimage;
	}
	return image;
}


/* PGM files allow a comment starting with '#' to end-of-line.  Skip
 *    white space including any comments.
 *    */
void SkipComments(FILE *fp)
{
	int ch;
	
	fscanf(fp," ");      /* Skip white space. */
	while ((ch = fgetc(fp)) == '#') {
		while ((ch = fgetc(fp)) != '\n'  &&  ch != EOF)
		   ;
	   fscanf(fp," ");
	}
	ungetc(ch, fp);      /* Replace last character read. */
}


/* Write an image to the file fp in PGM format.
 * */
void WritePGM(FILE *fp, LoweImage image)
{
	int r, c, val;
	
	fprintf(fp, "P5\n%d %d\n255\n", image->cols, image->rows);
	
	for (r = 0; r < image->rows; r++)
		for (c = 0; c < image->cols; c++) {
			val = (int) (255.0 * image->pixels[r][c]);
			fputc(MAX(0, MIN(255, val)), fp);
		}
}


/* Draw a white line from (r1,c1) to (r2,c2) on the image.  Both points
 *    must lie within the image.
 *    */
void DrawLine(LoweImage image, int r1, int c1, int r2, int c2)
{
	int i, dr, dc, temp;
	
	if (r1 == r2 && c1 == c2)  /* Line of zero length. */
		return;
	
	/* Is line more horizontal than vertical? */
	if (ABS(r2 - r1) < ABS(c2 - c1)) {
		
		/* Put points in increasing order by column. */
		if (c1 > c2) {
			temp = r1; r1 = r2; r2 = temp;
			temp = c1; c1 = c2; c2 = temp;
		}
		dr = r2 - r1;
		dc = c2 - c1;
		for (i = c1; i <= c2; i++)
							    	image->pixels[r1 + (i - c1) * dr / dc][i] = 1.0;
		
	} else {
		
		if (r1 > r2) {
			temp = r1; r1 = r2; r2 = temp;
			temp = c1; c1 = c2; c2 = temp;
		}
		dr = r2 - r1;
		dc = c2 - c1;
		for (i = r1; i <= r2; i++)
			image->pixels[i][c1 + (i - r1) * dc / dr] = 1.0;
	}
}


/*---------------------- Read keypoint file ---------------------------*/


/* This reads a keypoint file from a given filename and returns the list
 *    of keypoints.
 *    */
LoweKeypoint ReadKeyFile(char *filename)
{
	FILE *file;
	
	file = fopen (filename, "r");
	if (! file) FatalError("Could not open file: %s", filename);
	
	return ReadKeys(file);
}


/* Read keypoints from the given file pointer and return the list of
 *    keypoints.  The file format starts with 2 integers giving the total
 *       number of keypoints and the size of descriptor vector for each
 *          keypoint (currently assumed to be 128). Then each keypoint is
 *             specified by 4 floating point numbers giving subpixel row and
 *                column location, scale, and orientation (in radians from -PI to
 *                   PI).  Then the descriptor vector for each keypoint is given as a
 *                      list of integers in range [0,255].
 *
 *                      */
LoweKeypoint ReadKeys(FILE *fp)
{
	int i, j, num, len, val;
	LoweKeypoint k, keys = NULL;
	
	if (fscanf(fp, "%d %d", &num, &len) != 2)
			    	FatalError("Invalid keypoint file beginning.");
	
	if (len != 128)
		FatalError("LoweKeypoint descriptor length invalid (should be 128).");
	
	for (i = 0; i < num; i++) 
	{
	/* Allocate memory for the keypoint. */
	k = (LoweKeypoint) malloc(sizeof(struct LoweKeypointSt));
	k->next = keys;
	keys = k;
	k->ivec = (unsigned char *)malloc(len);
	
	if (fscanf(fp, "%f %f %f %f", &(k->row), &(k->col), &(k->scale),
				   &(k->ori)) != 4)
	  	FatalError("Invalid keypoint file format.");
		
	for (j = 0; j < len; j++) {
		if (fscanf(fp, "%d", &val) != 1 || val < 0 || val > 255)
			FatalError("Invalid keypoint file value.");
		k->ivec[j] = (unsigned char) val;
		}
	}
	return keys;
}



