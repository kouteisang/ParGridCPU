# Resource-Efficient FirmCore Decomposition on Billion-scale Multilayer Graphs (under review)

## Abstract


Multilayer (ML) graphs offer a convenient paradigm for modeling complex node-to-node interactions, such as social or semantic connections, as layers of a graph. In such graphs, FirmCore decomposition represents a convenient technique to identify cohesive groups of nodes with strong ties across layers. Unfortunately, the fastest FirmCore decomposition method fails to fully harness the resources leading to thread underutilization and idle threads. Our main observation is that FirmCores enjoy a grid structure we call fcgrid, which we exploit to distribute work among threads in our work: we introduce serial and parallel algorithms for multi-core CPUs, as well as the first GPU-based algorithm. Owing to our new structural design, our solutions show greatly improved performance and resource utilization. Our experiments on 12 datasets show 9× speedup on average for our serial version fcgrid compared to existing serial methods. Furthermore, our parallel algorithm achieves an average 100.3× speedup over the state-of-the-art parallel algorithm. For the challenging NP-hard densest subgraph mining problem in ML graphs, our algorithms achieve 15× speedup on average.


## Compile the code

```
cd Build
cmake ..
make
```

## Parameters

dataset(-d):


- homo
- sacchcere
- fao
- internet
- wiki
- obamainisrael
- amazon
- dblp-coauthor
- flickr-growth
- so
- hw1
- hw2


method(-m):

- lpe (SerGrid with LPE traversal order)
- dpe (SerGrid with DPE traversal order)
- ParGrid (ParCore algorithm)
- ParGridPlus (ParGridPlus)
- CoreIndex (Baseline)
- ParCoreIndex (parallel version of CoreIndex Baseline)

num of thread (-num_thread)

Order of layer (-o)

- 0: random order
- 1: order by increasing the number of edge
- 2: order by decreasing the number of edge
- 3: order by increasing the graph density
- 4: order by decreasing the graph density



## example

To run the hw1 dataset with ParGrid algorithm with 2 threads.

```
./firmcore_baseline -d hw1 -m ParGrid -num_thread 2
```

## Dataset

The datasets (H) Homo, (S) Sacchcere, (F) Fao, (I) Internet, (W) Wiki, (OI) ObamaInIsrael, (A) Amazon, (DL) DBLP-Large, (FG) FlickrGrowth, and (SO) StackOverflow are available at:

https://drive.google.com/file/d/1msDDYSIC2sQODqrBjpa55yUIm15D2bda/view?usp=sharing

The datasets (EM) Email-EuAll and (LJ) LiveJournal, can be downloaded from

https://snap.stanford.edu/data/index.html

The datasets (HW1) Hollywood-1 and (HW2) Hollywood-2 can be downloaded from:

https://law.di.unimi.it/datasets.php

The Google+ dataset used in our paper is available at:

https://people.duke.edu/~zg70/gplus.html

https://drive.google.com/file/d/1HF8Q2N_hxsaQ26MarKYxZEQhqI66qAxV/view
