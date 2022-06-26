# sgat
This repo contains the implementation of many sequence-to-graph alignment algorithms. 

## Installation

### Install dependencies
CMake >=3.3.

### Download and install
First get the repo:
```
git clone git@github.com:haowenz/sgat.git
```
Then just run:
```
cd sgat && mkdir build && cd build
cmake -DSGAL_BUILD_TESTING=ON ..
make
```
This will build the tests. Then the tests can be run:
```
ctest
```

## Usage
Sgat is a C++ library. It supports loading graphs in gfa format and sequence files in fasta/q format. You can easily use its API to build your own applications. An example on how to use the library to align sequences to graphs is given. After the installation, you can simply test it by running:
```
./sgat_example graph_file read_file
```
