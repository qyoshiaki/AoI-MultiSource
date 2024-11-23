# Computing the Mean Age of Information (AoI) in the PH+M/GI+GI/1 Queue

## Overview
This repository contains a C++ implementation of the computational algorithm for the mean Age of Information (AoI) developed in the following paper:

> Y. Inoue and M. Mandjes, "Characterizing the Age of Information with Multiple Coexisting Data Streams," arXiv:2404.15623, [https://arxiv.org/abs/2404.15623](https://arxiv.org/abs/2404.15623), 2024.

## Usage

### Building the Program
To build the source code, run the following command:

```bash
make
```
*Note: This step is required only once.*

### Running the Program
You can run the program using the following command:

```bash
./a.out --EG=${EG} --EH=${EH} --EH_bg=${EH_bg} --CvG=${CvG} --CvH=${CvH} --CvH_bg=${CvH_bg} --lmd_bg=${lmd_bg} --mu=${mu} --output=${outfile}
```

#### Parameters
- **EG**: Mean inter-generation time
- **EH**: Mean service time of the tagged stream
- **EH_bg**: Mean service time of the background stream
- **CvG**: Coefficient of variation of inter-generation times
- **CvH**: Coefficient of variation of service times of the tagged stream
- **CvH_bg**: Coefficient of variation of service times of the background stream
- **lmd_bg**: Arrival rate of the background stream
- **mu**: Service rate
- **output**: Output file name

For Linux and Mac users, you can also use `run.sh` to set these parameters conveniently.
