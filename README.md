# Submatrix-Search-Parallel-Implementation

## Overview

Welcome to the Parallel Image Processing project! This software, focuses on optimizing image-related operations through parallel computing. The implementation utilizes MPI (Message Passing Interface) and OpenMP (Open Multi-Processing) to distribute tasks efficiently among processes, enhancing the overall performance of the image processing workflow.

## Project Description

### Execution Instructions

To run the program, follow these steps in the terminal:

```bash
mpicc prog.c -o p -lm -fopenmp
mpiexec -n 2 ./p
```

Feel free to adjust the number of processes, but ensure it doesn't exceed the number of pictures.

### Code Structure

- The master process (Process 0) serves as the coordinator, handling data distribution and result aggregation.
- Image and object data are organized in a three-dimensional array (`int*** picArr`) for streamlined processing.
- The project employs MPI for inter-process communication and OpenMP for parallelization within each process.

### Parallelization Strategy

- The number of processes must be less than or equal to the number of images to ensure optimal workload distribution.
- Process 0 orchestrates the distribution of data among processes using MPI_Scatterv.
- Each process independently executes computations on its assigned images, optimizing performance through OpenMP.

### Algorithm Overview

- The core algorithm involves six nested loops, iterating over images, objects, rows, and columns to calculate matching scores.
- Parallelization is achieved using `#pragma omp parallel for` to distribute the workload among threads.

### Project Instructions

- The master process (Process 0) receives data from the "input" file, distributing the workload among processes.
- Data is stored in three-dimensional arrays for both images and objects.
- The execution workflow involves parallel processing, result collection, and final output generation.

## Conclusion

This project showcases an efficient approach to parallel image processing, leveraging the combined power of MPI and OpenMP. The strategic organization of data and the use of parallelization techniques contribute to the optimization of computation, making the software scalable and capable of handling large datasets.

Feel free to explore the code, experiment with different numbers of processes, and contribute to the project's development. For any issues or inquiries, please don't hesitate to reach out.
