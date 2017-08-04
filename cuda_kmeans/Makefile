CUDA_ARCH ?= sm_52
CUB_HOME ?= ./cub

NVCCFLAGS = -O4 --use_fast_math

all: cuda_main

cuda_main: cuda_main.cu labels.o cuda_io.o getopt.o kmeans.h
	nvcc $(NVCCFLAGS) -arch=$(CUDA_ARCH) -Xptxas -v -I$(CUB_HOME) -o cuda_main cuda_main.cu labels.o cuda_io.o getopt.o -lcublas -lineinfo

test: test.cu labels.o timer.o kmeans.h
	nvcc $(NVCCFLAGS) -arch=$(CUDA_ARCH) -Xptxas -v -I$(CUB_HOME) -o test test.cu labels.o timer.o -lcublas -lineinfo

labels.o: labels.cu labels.h
	nvcc $(NVCCFLAGS) -arch=$(CUDA_ARCH) -Xptxas -v -I$(CUB_HOME) -c -o labels.o labels.cu -lineinfo

cuda_io.o: cuda_io.cu
	nvcc $(NVCCFLAGS) -arch=$(CUDA_ARCH) -Xptxas -v -I$(CUB_HOME) -c -o cuda_io.o cuda_io.cu -lineinfo

getopt.o: getopt.cu
	nvcc $(NVCCFLAGS) -arch=$(CUDA_ARCH) -Xptxas -v -I$(CUB_HOME) -c -o getopt.o getopt.cu -lineinfo

timer.o:
	nvcc $(NVCCFLAGS) -arch=$(CUDA_ARCH) -c -o timer.o timer.cu -lineinfo

clean:
	rm -f *.o test
