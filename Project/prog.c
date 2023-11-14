#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<mpi.h>
#include<time.h>


// reading all pictures from the file Using a while loop together with the fscanf() function
void readingPictures(int index, int totalNumberOfPics, int picID, int size, int*** picsArr_0, int* picturesSizes_0, FILE* fp)
{
	while (index < totalNumberOfPics)//will run for each pic
	{
		fscanf(fp, "%d",&picID);
		fscanf(fp, "%d", &size);// get pic spesific size 
				
		// every cell in this array represent spesific rows of a spesific picture
		picsArr_0[picID - 1] = (int**)malloc(size * sizeof(int*));
		for (int row = 0; row < size; row++)
		{
			//row is the number of the row in spesific matrix of spesific picture
			// and every line creating malloc of the size number ( for example: 100)
			// every cell in this array represent spesific columns of a spesific picture
			picsArr_0[picID - 1][row] = (int*)malloc(size * sizeof(int));
		}

		// store the sizes of pictures in array
		picturesSizes_0[picID - 1] = size;
				
		//reading pictures matrixs
		for (int row = 0; row < size; row++)
		{
			for (int col = 0; col < size; col++)
			{
				fscanf(fp, "%d", &picsArr_0[picID - 1][row][col]);
			}
		}
		index++;
	}
}

// reading all objects from the file Using a while loop together with the fscanf() function
void readingObjects(int index, int totalNumberOfObjs, int objID, int size, int*** objsArr_0, int* objectsSizes_0, FILE* fp)
{
	while (index < totalNumberOfObjs)//will run for each obj
	{	
		fscanf(fp, "%d", &objID);
		fscanf(fp, "%d", &size);// get obj spesific size
				
		// every cell in this array represent spesific rows of a spesific object
		objsArr_0[objID - 1] = (int**)malloc(size * sizeof(int*));
		for (int row = 0; row < size; row++)
		{
			//row is the number of the row in spesific matrix of spesific object
			// and every line creating malloc of the size number ( for example: 20)
			// every cell in this array represent spesific columns of a spesific picture
			objsArr_0[objID - 1][row] = (int*)malloc(size * sizeof(int));
		}
		// store the sizes of objects in array
		objectsSizes_0[objID - 1] = size;
				
		//reading objects matrixs
		for (int row = 0; row < size; row++)
		{
			for (int col = 0; col < size; col++)
			{
				fscanf(fp, "%d", &objsArr_0[objID - 1][row][col]);
			}
		}
		index++;
	}
}

// free allocated memory for pictures and objects, pictures sizes and objects sizes
void freeAllocations(int*** picsArr_0, int* picturesSizes_0, int totalNumberOfPics, int*** objsArr_0, int* objectsSizes_0, int totalNumberOfObjs)
{
	for (int i = 0; i < totalNumberOfPics; i++)
	{
		for (int row = 0; row < picturesSizes_0[i]; row++)
		{
			free(picsArr_0[i][row]);
		}
		free(picsArr_0[i]);
	}
	free(picsArr_0);
	
	for (int i = 0; i < totalNumberOfObjs; i++)
	{
		for (int row = 0; row < objectsSizes_0[i]; row++)
		{
			free(objsArr_0[i][row]);
		}
		free(objsArr_0[i]);
	}
	free(objsArr_0);

	free(picturesSizes_0);
	free(objectsSizes_0);
}

//check if the number of pictures can split by all processes without remaining pics,
//and calculating number of pictures for each process,
//if there is remaining pics add them to last process
int splitPicsForEachProcess(int totalNumberOfPics, int numprocs, int rank)
{
	int picsPerProcess = 0;
	if (totalNumberOfPics % numprocs != 0)
	{
		picsPerProcess = totalNumberOfPics / numprocs;
		if (rank == numprocs - 1)
			picsPerProcess += totalNumberOfPics % numprocs; //add to last process the remaining pics
	}
	else
	{
		picsPerProcess = totalNumberOfPics / numprocs;
	}
	return picsPerProcess;
}

//allocate memory for pictures and objects for each process to dill with
void picsAndObjsAllocations(int*** picsArr, int picsPerProcess, int* picturesSizesArr, int*** objsArr, int totalNumberOfObjs, int* objectsSizesArr)
{
	//creating pictures array of that sizes for each process
	for (int i = 0; i < picsPerProcess; i++)
	{
		picsArr[i] = (int**)malloc(picturesSizesArr[i] * sizeof(int*));
		for (int j = 0; j < picturesSizesArr[i]; j++)
			picsArr[i][j] = (int*)malloc(picturesSizesArr[i] * sizeof(int));
	}
	
	//creating objects array of that sizes for each process
	for (int i = 0; i < totalNumberOfObjs; i++)
	{
		objsArr[i] = (int**)malloc(objectsSizesArr[i] * sizeof(int*));
		for (int j = 0; j < objectsSizesArr[i]; j++)
			objsArr[i][j] = (int*)malloc(objectsSizesArr[i] * sizeof(int));
	}
}

//matching between picture and object in specific start cell of the picture
int matching(int* objectsSizesArr, int difference, int*** picsArr, int*** objsArr, int x, int objNumber, int startCellRow, int startCellCol)
{
	for (int x2 = 0; x2 < objectsSizesArr[objNumber]; x2++)
	{	
		for (int y2 = 0; y2 < objectsSizesArr[objNumber]; y2++)
		{
			//main calculation statement
			difference = difference + abs(((float)picsArr[x][x2 + startCellRow][y2 + startCellCol] - (float)objsArr[objNumber][x2][y2]) / picsArr[x][x2 + startCellRow][y2 + startCellCol]);
		}
	}
	return difference;
}

int main(int argc, char* argv[])
{
	//measure time
	clock_t start, end;
	double runningTime = 0.0;
	start = clock();
		
	int rank, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	//limit threads to 4
	omp_set_num_threads(4);

	//Local processes variables and pointers for all ranks//
	int*** picsArr;// to store pictures, 3d array
	int*** objsArr; // to store objects, 3d array
	int* picturesSizesArr; // array to store picture sizes
	int* objectsSizesArr;// array to store object sizes
	int check = 0;// checkFlag to dill with the threads and to see if there is data in the file input

	//Rank 0 only- variables and pointers //
	int*** picsArr_0 = NULL;// to store pictures
	int*** objsArr_0 = NULL;// to store objects
	int* picturesSizes_0 = NULL;// array to store picture sizes only for process 0
	int* objectsSizes_0 = NULL;// array to store object sizes only for process 0

	float ending_condition = 0.0;// variable for ending condition (0.100000)
	
	// to store total number of pictures, objects and specific picture and object
	int totalNumberOfPics = 0, totalNumberOfObjs = 0, picID = 0, objID = 0;

	// Read from the text file by process 0//
	if (rank == 0)//only for process 0
	{
		int size;//size of spesific pics and after size of spesific objs
		int index = 0;
		FILE* fp = fopen("input.txt", "r");//create file object
		
		// checks if the file opened succesfully 
		if (fp != NULL)
		{
			//reading values:
			fscanf(fp, "%f", &ending_condition);
			fscanf(fp, "%d", &totalNumberOfPics);
			
			//allocation pics sizes:
			picturesSizes_0 = (int*)malloc(totalNumberOfPics * sizeof(int));
			//allocation pics array:
			// every cell in this array pointer each picture
			picsArr_0 = (int***)malloc(totalNumberOfPics * sizeof(int**));

			//reading all pictures from the file Using a while loop together with the fscanf() function
			readingPictures(index, totalNumberOfPics, picID, size, picsArr_0, picturesSizes_0, fp);
			
			//reading total number of objects:
			fscanf(fp, "%d", &totalNumberOfObjs);
			
			//allocation objs sizes array:
			objectsSizes_0 = (int*)malloc(totalNumberOfObjs * sizeof(int));
			//allocation objs array:
			objsArr_0 = (int***)malloc(totalNumberOfObjs * sizeof(int**));
			
			index = 0;
			
			//reading all objects from the file Using a while loop together with the fscanf() function
			readingObjects(index, totalNumberOfObjs, objID, size, objsArr_0, objectsSizes_0, fp);
			
			//checks if there is data in the file, 0 for empty, 1 for not empty
			check = 1;//not empty
		}
		else {
			check = 0;//empty
		}
		fclose(fp);//close file pointer
	}

	//sending from rank 0 to all processes to exit if there is nothing in file or program failed to open the file
	MPI_Bcast(&check, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (check != 1)
	{
		if (rank == 0)
		{
			printf("error- File not found\n");
			//rank 0 - free allocated memory for pictures and objects, pictures sizes and objects sizes
			freeAllocations(picsArr_0, picturesSizes_0, totalNumberOfPics, objsArr_0, objectsSizes_0, totalNumberOfObjs);
		}
		MPI_Finalize();
		exit(0);//terminate if there is nothing in file
	}

	//rank 0 send to all processes
	MPI_Bcast(&totalNumberOfPics, 1, MPI_INT, 0, MPI_COMM_WORLD);//send count of pictures to all processses
	MPI_Bcast(&totalNumberOfObjs, 1, MPI_INT, 0, MPI_COMM_WORLD);//send count of objects to all processes
	MPI_Bcast(&ending_condition, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);//seding ending condition to all process

	//program exit if the number of processes is higher than the number of pictures
	if (totalNumberOfPics < numprocs)
	{
		if (rank == 0)
		{
			printf("the number of processes must be less or equal to the number of pictures\n");
			//rank 0 - free allocated memory for pictures and objects, pictures sizes and objects sizes
			freeAllocations(picsArr_0, picturesSizes_0, totalNumberOfPics, objsArr_0, objectsSizes_0, totalNumberOfObjs);
		}
		MPI_Finalize();
		exit(0);
	}

	//if all goes well,
	//than working on sending the DATA and Matching!
	int picsPerProcess;
	int *sendcounts = NULL;//array describing how many elements to send to each process
	int *displs = NULL;// array describing the position where each segment begins
	MPI_Status *status = (MPI_Status*)malloc(sizeof(MPI_Status));//for the MPI_Recv
	
	//check if the number of pictures can split by all processes without remaining pics,
	//and calculating number of pictures for each process,
	//if there is remaining pics add them to last process
	picsPerProcess = splitPicsForEachProcess(totalNumberOfPics, numprocs, rank);
	
	if (rank == 0)
	{
		//sendcounts is a variable for the chunk sizes that each process receive
		sendcounts = (int*)malloc(numprocs * sizeof(int));
		sendcounts[0] = picsPerProcess;
		for (int i = 1; i < numprocs; i++)
			MPI_Recv(&sendcounts[i], 1, MPI_INT, i, 1, MPI_COMM_WORLD, status);//rank 0 receiving from each process how many pictures have to send to each process

		displs = (int*)malloc(numprocs * sizeof(int));
		displs[0] = 0;
		for (int i = 1; i < numprocs; i++)
			displs[i] = sendcounts[i - 1] + displs[i - 1]; //calculation of the position where each segment begins for each process
	}
	else
	{
		//sending from each process to rank 0 how many pictures rank 0 have to send each process
		MPI_Send(&picsPerProcess, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}

	//creating pictures sizes array and objects sizes for each process 
	picturesSizesArr = (int*)malloc(picsPerProcess * sizeof(int));
	objectsSizesArr = (int*)malloc(totalNumberOfObjs * sizeof(int));
	
	//rank 0 sending pictures sizes to all process by the chank of each process to manage
	MPI_Scatterv(picturesSizes_0, sendcounts, displs, MPI_INT, picturesSizesArr, picsPerProcess, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (rank == 0)
	{
		//sending to each processes all of the objects sizes
		for (int i = 1; i < numprocs; i++)
		{
			MPI_Send(objectsSizes_0, totalNumberOfObjs, MPI_INT, i, 1, MPI_COMM_WORLD);
		}
		objectsSizesArr = objectsSizes_0;
	}
	else
	{
		//receiving to each processes all of the objects sizes from rank 0
		MPI_Recv(objectsSizesArr, totalNumberOfObjs, MPI_INT, 0, 1, MPI_COMM_WORLD,status);
	}
	
	//allocate memory for pictures and objects for each process to dill with
	picsArr = (int***)malloc(picsPerProcess * sizeof(int**));
	objsArr = (int***)malloc(totalNumberOfObjs * sizeof(int**));
	picsAndObjsAllocations(picsArr, picsPerProcess, picturesSizesArr, objsArr, totalNumberOfObjs, objectsSizesArr);
	

	//sending to each precess the spesific chunk of pictures from process 0 
	if (rank == 0)
	{
		for (int i = 1; i < numprocs; i++)
		{
			int position = displs[i];
			for (int j = 0; j < sendcounts[i]; j++)
			{
				for (int row = 0; row < picturesSizes_0[position]; row++)
				{
					for (int col = 0; col < picturesSizes_0[position]; col++)
					{
						MPI_Send(&picsArr_0[position][row][col], 1, MPI_INT, i, 1, MPI_COMM_WORLD);
					}
				}
				position++;
			}
		}
		
		//sending to each precess all objects from process 0 
		for (int i = 1; i < numprocs; i++)
		{
			for (int j = 0; j < totalNumberOfObjs; j++)
			{
				for (int row = 0; row < objectsSizesArr[j]; row++)
				{
					for (int col = 0; col < objectsSizesArr[j]; col++)
					{
						MPI_Send(&objsArr_0[j][row][col], 1, MPI_INT, i, 2, MPI_COMM_WORLD);
					}
				}
			}
		}
		
		//initialize rank 0 with the proper chank of pictures that he need to dill with
		for (int i = 0; i < picsPerProcess; i++)
		{
			for (int row = 0; row < picturesSizesArr[i]; row++)
			{
				for (int col = 0; col < picturesSizesArr[i]; col++)
				{
					picsArr[i][row][col] = picsArr_0[i][row][col];
				}
			}
		}
		
		//initialize rank 0 with all objects
		for (int i = 0; i < totalNumberOfObjs; i++)
		{
			for (int row = 0; row < objectsSizesArr[i]; row++)
			{
				for (int col = 0; col < objectsSizesArr[i]; col++)
				{
					objsArr[i][row][col] = objsArr_0[i][row][col];
				}
			}
		}
	}
	//each process receiving pictures and objects from process 0
	else
	{
		//each process receiving the chunk of pictures he need to dill with from process 0
		for (int i = 0; i < picsPerProcess; i++)
		{
			for (int row = 0; row < picturesSizesArr[i]; row++)
			{
				for (int col = 0; col < picturesSizesArr[i]; col++)
				{
					MPI_Recv(&picsArr[i][row][col], 1, MPI_INT, 0, 1, MPI_COMM_WORLD, status);
				}
			}
		}
		
		//each process receiving all objects
		for (int i = 0; i < totalNumberOfObjs; i++)
		{
			for (int row = 0; row < objectsSizesArr[i]; row++)
			{
				for (int col = 0; col < objectsSizesArr[i]; col++)
				{
					MPI_Recv(&objsArr[i][row][col], 1, MPI_INT, 0, 2, MPI_COMM_WORLD, status);
				}
			}
		}
	}
	
	//synchronizing processes by stopping each one when they get to the Barrier
	MPI_Barrier(MPI_COMM_WORLD);
	
	float difference = 0;//to store diffrence value between picture and object in a specific cell
	float sumDifference = 0.0;//to store all diffrences value between picture and object in all cells of specific position
	float matchingValue = 0.0;//to store all diffrences value after deviding by size of all object matrix in a specific position
	int check3ObjInPic; // check if there is 3 objects for each picture
	int foundObjFlag = 0; ////any object found in picture
	int exactSize;
	int** objsNumInPicsPerProcessArr = NULL;
	int** objRowIndexsInPicsPerProcessArr = NULL;
	int** objColIndexsInPicsPerProcessArr = NULL;
	int* numberObjsFoundInPicsPerProcessArr = NULL;
	
	//arrays to store processing output 
	objsNumInPicsPerProcessArr = (int**)malloc(picsPerProcess * sizeof(int*));
	objRowIndexsInPicsPerProcessArr = (int**)malloc(picsPerProcess * sizeof(int*));
	objColIndexsInPicsPerProcessArr = (int**)malloc(picsPerProcess * sizeof(int*));
	numberObjsFoundInPicsPerProcessArr = (int*)malloc(picsPerProcess * sizeof(int));
	for (int i = 0; i < picsPerProcess; i++)
	{
		objsNumInPicsPerProcessArr[i] = (int*)malloc(3 * sizeof(int));
		objRowIndexsInPicsPerProcessArr[i] = (int*)malloc(3 * sizeof(int));
		objColIndexsInPicsPerProcessArr[i] = (int*)malloc(3 * sizeof(int));
	}

	//initializing arrays to save processing data per picture
	int matchingObjectsNumberArr[3] = { -1,-1,-1 };
	int rowIndexsStartCellsArr[3] = { -1,-1,-1 };
	int colIndexsStartCellsArr[3] = { -1,-1,-1 };
	
	int x, i, objNumber, startCellRow, startCellCol, j;
	
	//main loop- each process pictures will be parallelaized with threads,
	//thread for each picture
	#pragma omp parallel for private(x, i, objNumber, foundObjFlag, check3ObjInPic, startCellRow, exactSize, sumDifference, matchingValue)
	//will run for the specific number of pictures for each process to dill with
	for ( x = 0; x < picsPerProcess; x++)
	{
		foundObjFlag = 0;
		check3ObjInPic = 0;
		
		//initializing arrays again to save processing data per picture
		for ( i = 0; i < 3; i++)
		{
			matchingObjectsNumberArr[i] = -1;
			rowIndexsStartCellsArr[i] = -1;
			colIndexsStartCellsArr[i] = -1;
		}
		
		//will run for each picture on total number of objects
		for ( objNumber = 0; objNumber < totalNumberOfObjs; objNumber++) 
		{
			foundObjFlag = 0;
			
			//"exactSize" variable runs on the appropriate indexes to avoid overflow,
			//that the object wouldn't go outside of the picture
			//when doing the calculation on the matching
			exactSize = picturesSizesArr[x] - objectsSizesArr[objNumber] + 1;
			for ( startCellRow = 0; startCellRow < exactSize; startCellRow++)
			{
				for ( startCellCol = 0; startCellCol < exactSize; startCellCol++)
				{
					//to store diffrence value between picture and object in a specific cell
					difference = 0.0;
					//to store all diffrences value between picture and object in all cells of specific position
					sumDifference = 0.0;
					//to store all diffrences value after deviding by size of all object matrix in a specific position
					matchingValue = 0.0;
					
					sumDifference = matching(objectsSizesArr, difference, picsArr, objsArr, x, objNumber, startCellRow, startCellCol);
					
					//divide with object size
					matchingValue = (float)sumDifference / (float)(objectsSizesArr[objNumber] * objectsSizesArr[objNumber]);
					
					//checking ending condition (if lower than 0.100000)
					if (matchingValue < ending_condition)
					{
						//save the object number that is lower than ending condition
						matchingObjectsNumberArr[check3ObjInPic] = objNumber;
						//save the row indexs of starting cells 
						rowIndexsStartCellsArr[check3ObjInPic] = startCellRow;
						//save the col indexs of starting cells 
						colIndexsStartCellsArr[check3ObjInPic] = startCellCol;
						
						//to check there must be three objects
						check3ObjInPic++;
						//any object found
						foundObjFlag = 1;
					}

					//check if there is three objects in picture,
					//if yes end the search for this picture
					if (check3ObjInPic == 3)
					{
						objNumber = totalNumberOfObjs;
					}
					//break, to check next object
					if (foundObjFlag == 1)
						startCellCol = exactSize;
				}
				//break, to check next object
				if (foundObjFlag == 1)
					startCellRow = exactSize;
			}
		}
		
		//temporal printing to terminal to see the results:
		//means there is no three objects
		if (check3ObjInPic != 3)
		{
			printf("Rank %d Picture %d: No three different objects were found\n",rank,x+1);
		}
		//there are three objects founded in picture
		else if (check3ObjInPic == 3) 
		{
			//printing the objects and position
			printf("Rank %d Picture %d: Found Objects:",rank,x+1);
			for (int x = 0; x < 3; x++) {
				printf(" %d Position (%d,%d);",matchingObjectsNumberArr[x]+1,rowIndexsStartCellsArr[x],colIndexsStartCellsArr[x]);

			}printf("\n");
		}
		
		for ( j = 0; j < 3; j++)
		{
			//saving all objects numbers for each picture in array of all pictures in specific process
			objsNumInPicsPerProcessArr[x][j] = matchingObjectsNumberArr[j];
			//saving all row indexs of starting cells for each picture in array of all pictures in specific process
			objRowIndexsInPicsPerProcessArr[x][j] = rowIndexsStartCellsArr[j];
			//saving all col indexs of starting cells for each picture in array of all pictures in specific process
			objColIndexsInPicsPerProcessArr[x][j] = colIndexsStartCellsArr[j];
		}
		//saving the number of objects found in picture in array of all pictures in specific process
		numberObjsFoundInPicsPerProcessArr[x] = check3ObjInPic;
		
	}
	
	//only for process 0, to receive the results from all processes,
	//and writing into the "output" file
	if (rank == 0)
	{
		//saving all objects numbers for each picture in array of all pictures in all processes
		int*** objsNumInPicsAllProcessesArr = NULL;
		//saving all row indexs of starting cells for each picture in array of all pictures in all processes
		int*** objRowIndexsInPicsAllProcessArr = NULL;
		//saving all col indexs of starting cells for each picture in array of all pictures in all processes
		int*** objColIndexsInPicsAllProcessArr = NULL;
		
		//saving the number of objects found in picture in array of all pictures in all processes
		int** numberObjsFoundInPicsAllProcessArr = NULL;

		//allocating memory of arrays that receive the data
		objsNumInPicsAllProcessesArr = (int***)malloc(numprocs * sizeof(int**));
		objRowIndexsInPicsAllProcessArr = (int***)malloc(numprocs * sizeof(int**));
		objColIndexsInPicsAllProcessArr = (int***)malloc(numprocs * sizeof(int**));
		numberObjsFoundInPicsAllProcessArr = (int**)malloc(numprocs * sizeof(int*));
		for (int i = 0; i < numprocs; i++)
		{
			objsNumInPicsAllProcessesArr[i] = (int**)malloc(sendcounts[i] * sizeof(int*));
			objRowIndexsInPicsAllProcessArr[i] = (int**)malloc(sendcounts[i] * sizeof(int*));
			objColIndexsInPicsAllProcessArr[i] = (int**)malloc(sendcounts[i] * sizeof(int*));
			numberObjsFoundInPicsAllProcessArr[i] = (int*)malloc(sendcounts[i] * sizeof(int));
			for (int j = 0; j < sendcounts[i]; j++)
			{
				objsNumInPicsAllProcessesArr[i][j] = (int*)malloc(3 * sizeof(int));
				objRowIndexsInPicsAllProcessArr[i][j] = (int*)malloc(3 * sizeof(int));
				objColIndexsInPicsAllProcessArr[i][j] = (int*)malloc(3 * sizeof(int));
			}
		}

		//initializing the arrays with rank 0 values
		numberObjsFoundInPicsAllProcessArr[0] = numberObjsFoundInPicsPerProcessArr;
		for (int i = 0; i < picsPerProcess; i++)
		{
			objsNumInPicsAllProcessesArr[0][i] = objsNumInPicsPerProcessArr[i];
			objRowIndexsInPicsAllProcessArr[0][i] = objRowIndexsInPicsPerProcessArr[i];
			objColIndexsInPicsAllProcessArr[0][i] = objColIndexsInPicsPerProcessArr[i];	
		}

		//initializing the arrays with all other processes values
		//by recieving the information with "MPI_Recv" from all processes
		for (int i = 1; i < numprocs; i++)
		{
			MPI_Recv(numberObjsFoundInPicsAllProcessArr[i], sendcounts[i], MPI_INT, i, 1, MPI_COMM_WORLD, status);
			
		}
		for (int i = 1; i < numprocs; i++)
		{
			for (int j = 0; j < sendcounts[i]; j++)
			{
				MPI_Recv(objsNumInPicsAllProcessesArr[i][j], 3, MPI_INT, i, 1, MPI_COMM_WORLD, status);//receiving data of object found in picture
				MPI_Recv(objRowIndexsInPicsAllProcessArr[i][j], 3, MPI_INT, i, 1, MPI_COMM_WORLD, status);//receiving data of row indexs
				MPI_Recv(objColIndexsInPicsAllProcessArr[i][j], 3, MPI_INT, i, 1, MPI_COMM_WORLD, status);//receiving data of col indexs
			}
		}

		//writing the results in the "output" file
		FILE* fp = fopen("output.txt", "w");
		int pictureNum = 1;
		//writing result in file
		for (int i = 0; i < numprocs; i++)
		{
			for (int j = 0; j < sendcounts[i]; j++)
			{
				if (numberObjsFoundInPicsAllProcessArr[i][j] != 3)
				{
					fprintf(fp, "Picture %d: No three different Objects were found\n", pictureNum);
				}
				else
				{
					fprintf(fp, "Picture %d: found Objects: ",pictureNum);
					for (int k = 0; k < 3; k++)
					{
						fprintf(fp," %d Position (%d,%d); ", objsNumInPicsAllProcessesArr[i][j][k]+1, objRowIndexsInPicsAllProcessArr[i][j][k],objColIndexsInPicsAllProcessArr[i][j][k]);
					}
					fprintf(fp, "\n");
				}
				pictureNum++;
			}
		}
		fclose(fp);
		
		//free all allocated memory
		for (int i = 0; i < numprocs; i++)
		{
			for (int j = 0; j < sendcounts[i]; j++)
			{
				free(objsNumInPicsAllProcessesArr[i][j]);
				free(objRowIndexsInPicsAllProcessArr[i][j]);
				free(objColIndexsInPicsAllProcessArr[i][j]);
			}	
			free(objsNumInPicsAllProcessesArr[i]);
			free(objRowIndexsInPicsAllProcessArr[i]);
			free(objColIndexsInPicsAllProcessArr[i]);
			free(numberObjsFoundInPicsAllProcessArr[i]);
		}
		free(objsNumInPicsAllProcessesArr);
		free(objRowIndexsInPicsAllProcessArr);
		free(objColIndexsInPicsAllProcessArr);
		free(numberObjsFoundInPicsAllProcessArr);
	}

	//all other processes are sending the findings information to process 0
	else
	{
		MPI_Send(numberObjsFoundInPicsPerProcessArr, picsPerProcess, MPI_INT, 0, 1, MPI_COMM_WORLD);

		for (int i = 0; i < picsPerProcess; i++)
		{
			MPI_Send(objsNumInPicsPerProcessArr[i], 3, MPI_INT, 0, 1, MPI_COMM_WORLD);
			MPI_Send(objRowIndexsInPicsPerProcessArr[i], 3, MPI_INT, 0, 1, MPI_COMM_WORLD);
			MPI_Send(objColIndexsInPicsPerProcessArr[i], 3, MPI_INT, 0, 1, MPI_COMM_WORLD);
		}
	}
	
	
	//all processes - free allocated memory for pictures and objects, pictures sizes and objects sizes
	freeAllocations(picsArr, picturesSizesArr, picsPerProcess, objsArr, objectsSizesArr, totalNumberOfObjs);
	
	//free all allocated memory
	if (!rank == 0)
	{
		for (int i = 0; i < picsPerProcess; i++)
		{
			free(objsNumInPicsPerProcessArr[i]);
			free(objRowIndexsInPicsPerProcessArr[i]);
			free(objColIndexsInPicsPerProcessArr[i]);
		}
		free(objsNumInPicsPerProcessArr);
		free(objRowIndexsInPicsPerProcessArr);
		free(objColIndexsInPicsPerProcessArr);
		free(numberObjsFoundInPicsPerProcessArr);
	}
	
	free(sendcounts);
	free(displs);
	
	MPI_Finalize();
	
	//for measure the time
	end = clock();
	runningTime = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("---------------------\n");
	printf("running time: %f\n", runningTime);
	printf("---------------------\n");
	
	return 0;
}
