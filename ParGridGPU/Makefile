NVCC = nvcc
TARGET = main
SRC = Graph/Graph.cpp Graph/MultilayerGraph.cpp Algorithm/ParGridGPU_block.cu Algorithm/ParGridGPU.cu Algorithm/ParGridGPU_k.cu main.cu

all:
	$(NVCC) -O3 -std=c++11 -I./Graph $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)
