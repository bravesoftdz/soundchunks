#include <thrust/device_vector.h>
#include "timer.h"
#include <iostream>
#include "cuda.h"
#include <cstdlib>
#include "kmeans.h"

template<typename T>
void fill_array(T& array, int m, int n) {
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < n; j++) {
      array[i * n + j] = (i % 2)*3 + j;
    }
  }
}

template<typename T>
void random_data(thrust::device_vector<T>& array, int n, int d, int k) {
  thrust::host_vector<T> host_array(n*d);
  for(int i = 0; i < n; i++) {
  for(int j = 0; j < d; j++) {
    //    host_array[i] = (T)rand()/(T)RAND_MAX;
    host_array[i*d+j] = i%k;
    //    host_array[j*n+i] = i%k;
  }
  }
  array = host_array;
}

void random_labels(thrust::device_vector<int>& labels, int n, int k) {
  thrust::host_vector<int> host_labels(n);
  for(int i = 0; i < n; i++) {
    host_labels[i] = rand() % k;
  }
  labels = host_labels;
}

typedef float real_t;

int main() {
  int max_iterations = 10000;
  int n = 260753;
  //  int d = 298;
  //int k = 100;
  int d = 3;
  int k = 10;
  double thresh = 1e-3;

  int n_gpu;
  cudaGetDeviceCount(&n_gpu);
  n_gpu=1;

  std::cout << n_gpu << " gpus." << std::endl;

  thrust::device_vector<real_t> *data[16];
  thrust::device_vector<int> *labels[16];
  thrust::device_vector<real_t> *centroids[16];
  thrust::device_vector<real_t> *distances[16];
  for (int q = 0; q < n_gpu; q++) {
    cudaSetDevice(q);
    data[q] = new thrust::device_vector<real_t>(n/n_gpu*d);
    labels[q] = new thrust::device_vector<int>(n/n_gpu*d);
    centroids[q] = new thrust::device_vector<real_t>(k * d);
    distances[q] = new thrust::device_vector<real_t>(n);
  }

  std::cout << "Generating random data" << std::endl;
  std::cout << "Number of points: " << n << std::endl;
  std::cout << "Number of dimensions: " << d << std::endl;
  std::cout << "Number of clusters: " << k << std::endl;
  std::cout << "Max. number of iterations: " << max_iterations << std::endl;
  std::cout << "Stopping threshold: " << thresh << std::endl;

  /* Intializes random number generator */
  //srand((unsigned) time(&t));
  srand(777);
  
  for (int q = 0; q < n_gpu; q++) {
    random_data<real_t>(*data[q], n/n_gpu, d, k);
    random_labels(*labels[q], n/n_gpu, k);
  }
  kmeans::timer t;
  t.start();
  kmeans::kmeans<real_t>(n, d, k, data, labels, centroids, distances, n_gpu, max_iterations, true, thresh);
  float time = t.stop();
  std::cout << "  Time: " << time/1000.0 << " s" << std::endl;

  // debug
  int printcenters=1;
  if(printcenters){
    thrust::host_vector<real_t> *ctr = new thrust::host_vector<real_t>(*centroids[0]);
    for(unsigned int ii=0;ii<k;ii++){
      fprintf(stderr,"ii=%d of k=%d ",ii,k);
      for(unsigned int jj=0;jj<d;jj++){
        fprintf(stderr,"%g ",(*ctr)[d*ii+jj]);
      }
      fprintf(stderr,"\n");
      fflush(stderr);
    }
  }
  
  for (int q = 0; q < n_gpu; q++) {
    delete(data[q]);
    delete(labels[q]);
    delete(centroids[q]);
    delete(distances[q]);
  }
}
