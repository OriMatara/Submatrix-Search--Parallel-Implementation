# Submatrix-Search-Parallel-Implementation

Parallel Image Processing
Overview
This project focuses on parallel image processing using MPI (Message Passing Interface) and OpenMP (Open Multi-Processing). The program is designed to distribute the workload across multiple processes to enhance the efficiency of image-related operations.

Table of Contents
Getting Started
Prerequisites
Installation
Running the Program
Execution Instructions
Contributing
License
Contact
Getting Started
Prerequisites
MPI (Message Passing Interface)
OpenMP (Open Multi-Processing)
Installation
Clone the repository:

bash
Copy code
git clone https://github.com/[YourUsername]/parallel-image-processing.git
Navigate to the project directory:

bash
Copy code
cd parallel-image-processing
Compile the program:

bash
Copy code
mpicc prog.c -o p -lm -fopenmp
Running the Program
To run the program with 2 MPI processes:

bash
Copy code
mpiexec -n 2 ./p
Adjust the number of MPI processes based on your system capabilities, but ensure it doesn't exceed the number of pictures.

Contributing
If you would like to contribute to this project, please follow these guidelines:

Submit bug reports or feature requests through the issue tracker.
Fork the repository, create a branch, and submit a pull request for code contributions.
License
This project is licensed under the MIT License - see the LICENSE file for details.

Contact
John Doe - john.doe@email.com

Project Link: https://github.com/[YourUsername]/parallel-image-processing
