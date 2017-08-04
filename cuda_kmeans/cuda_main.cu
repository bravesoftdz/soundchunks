/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         seq_main.c   (an sequential version)                      */
/*   Description:  This program shows an example on how to call a subroutine */
/*                 that implements a simple k-means clustering algorithm     */
/*                 based on Euclid distance.                                 */
/*   Input file format:                                                      */
/*                 ascii  file: each line contains 1 data object             */
/*                 binary file: first 4-byte integer is the number of data   */
/*                 objects and 2nd integer is the no. of features (or        */
/*                 coordinates) of each object                               */
/*                                                                           */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department Northwestern University                         */
/*            email: wkliao@ece.northwestern.edu                             */
/*   Copyright, 2005, Wei-keng Liao                                          */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Copyright (c) 2005 Wei-keng Liao
// Copyright (c) 2011 Serban Giuroiu
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// -----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>     /* strtok() */
#include <sys/types.h>  /* open() */
#include <sys/stat.h>
#include <fcntl.h>
#include "unistd.h"     /* getopt() */
#include "kmeans.h"

int      _debug;

#define malloc2D(name, xDim, yDim, type) do {               \
    name = (type **)malloc(xDim * sizeof(type *));          \
    assert(name != NULL);                                   \
    name[0] = (type *)malloc(xDim * yDim * sizeof(type));   \
    assert(name[0] != NULL);                                \
    for (size_t i = 1; i < xDim; i++)                       \
        name[i] = name[i-1] + yDim;                         \
} while (0)

float** file_read(int, char*, int*, int*);
int     file_write(char*, int, int, int, float**, int*);

/*---< usage() >------------------------------------------------------------*/
static void usage(char *argv0, float threshold) {
    const char *help =
        "Usage: %s [switches] -i filename -n num_clusters\n"
        "       -i filename    : file containing data to be clustered\n"
        "       -c filename    : file containing init centroids\n"
        "       -b             : input file is in binary format (default no)\n"
        "       -n num_clusters: number of clusters (K must > 1)\n"
        "       -t threshold   : threshold value (default %.4f)\n"
        "       -d             : enable debug mode\n";
    fprintf(stderr, help, argv0, threshold);
    exit(-1);
}

typedef float real_t;

/*---< main() >-------------------------------------------------------------*/
int main(int argc, char **argv) {
           int     opt;
    extern char   *optarg;
    extern int     optind;
           int     isBinaryFile;

           int     numClusters, numCoords, numObjs;
           int    *membership;    /* [numObjs] */
           char   *filename, *centFname;
           float **objects;       /* [numObjs][numCoords] data objects */
           float **clusters;      /* [numClusters][numCoords] cluster center */
           float   threshold;

    /* some default values */
    _debug           = 0;
    threshold        = 0.001;
    numClusters      = 0;
    isBinaryFile     = 0;
    filename         = NULL;
    centFname        = NULL;

    setbuf(stdout, NULL);
    setbuf(stderr, NULL);
    
    while ( (opt=getopt(argc,argv,"p:i:c:n:t:abdo"))!= EOF) {
        switch (opt) {
            case 'i': filename=optarg;
                      break;
            case 'c': centFname=optarg;
                      break;
            case 'b': isBinaryFile = 1;
                      break;
            case 't': threshold=atof(optarg);
                      break;
            case 'n': numClusters = atoi(optarg);
                      break;
            case 'd': _debug = 1;
                      break;
            case '?': usage(argv[0], threshold);
                      break;
            default: usage(argv[0], threshold);
                      break;
        }
    }

    if (filename == 0 || numClusters <= 1) usage(argv[0], threshold);

    /* read data points from file ------------------------------------------*/
    objects = file_read(isBinaryFile, filename, &numObjs, &numCoords);
    if (objects == NULL) exit(1);

    if (centFname != 0)
    {
      int yc, xc;
      clusters = file_read(isBinaryFile, centFname, &yc, &xc);
      if (yc != numClusters || xc != numCoords)
      {
        printf("Cendroids mismatch: numCoords %d->%d numClusters %d->%d\n", numCoords, xc, numClusters, yc);
      }
    }
    else
    {
      clusters = NULL;
    }

    /* start the timer for the core computation -----------------------------*/
    /* membership: the cluster id for each data object */
    membership = (int*) malloc(numObjs * sizeof(int));
    assert(membership != NULL);

    for (int i = 0; i < numObjs; i++) {
        membership[i] = i;
    }

    thrust::device_vector<real_t> *data;
    thrust::device_vector<int> *labels;
    thrust::device_vector<real_t> *centroids;
    thrust::device_vector<real_t> *distances;
    
    cudaSetDevice(0);
    data = new thrust::device_vector<real_t>(numObjs*numCoords);
    labels = new thrust::device_vector<int>(numObjs);
    centroids = new thrust::device_vector<real_t>(numClusters * numCoords);
    distances = new thrust::device_vector<real_t>(numObjs);
    
    thrust::copy(&(objects[0][0]), &(objects[numObjs-1][numCoords-1]), data->begin());

    if (clusters != 0)
    {
      thrust::copy(&(clusters[0][0]), &(clusters[numClusters-1][numCoords-1]), centroids->begin());
    }

    thrust::copy(&(membership[0]), &(membership[numObjs-1]), labels->begin());
    
    printf("numObjs %d numCoords %d numClusters %d\n", numObjs, numCoords, numClusters);
      
    kmeans::kmeans<real_t>(numObjs, numCoords, numClusters, &data, &labels, &centroids, &distances, 1, INT_MAX, clusters == 0, threshold);

    malloc2D(clusters, numClusters, numCoords, float);
    thrust::copy(centroids->begin(), centroids->end(), &(clusters[0][0]));
    
    thrust::copy(labels->begin(), labels->end(), &(membership[0]));
    
    free(objects[0]);
    free(objects);

    /* output: the coordinates of the cluster centres ----------------------*/
    file_write(filename, numClusters, numObjs, numCoords, clusters,
               membership);

    free(membership);
    free(clusters[0]);
    free(clusters);

    return(0);
}

