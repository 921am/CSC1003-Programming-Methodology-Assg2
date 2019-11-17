#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define GNUPLOT "/usr/local/Cellar/gnuplot/5.2.7_1/bin/gnuplot -persist" // change this based on where your gnuplot executable file is located. To locate it, type `brew ls gnuplot` in the terminal
#define N 10000 // The amount of values in the dataset

//test path: /Users/magdalene/Desktop/SIT-UofG/programMeth/progMeth_assg2/Group9_15.txt

char DATASET_FILEPATH[1000];
float X_VALUES[N], Y_VALUES[N], N_VALUES[N];

void readDataFromFile();
void plotGraph(float slope, float yIntercept);

//swaps 2 elements
void swap(float* a, float* b) 
{ 
    float t = *a; 
    *a = *b; 
    *b = t; 
} 
  
//uses last element as the pivot
//elements < pivot will be placed on the left of pivot
//elements > pivot will be placed on the right of pivot
int partition (float sortN[], int low, int high) 
{ 
    float pivot = sortN[high];    // pivot 
    int i = (low - 1);  // Index of smaller element 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        //if current element is smaller than the pivot
        if (sortN[j] < pivot) 
        { 
            i++;    //increment index of smaller element by one
            swap(&sortN[i], &sortN[j]); 
        } 
    } 
    swap(&sortN[i + 1], &sortN[high]); 
    return (i + 1); 
} 
  
//sort N_VALUES using quicksort
//low represents starting index and high is the ending index
void quickSort(float sortN[], int low, int high) 
{ 
    if (low < high) 
    { 
        //puts sortN[] in the right position
        int partitioningIndex = partition(sortN, low, high); 
  
        //Sort the elements before and after partition seperatelty
        quickSort(sortN, low, partitioningIndex - 1); 
        quickSort(sortN, partitioningIndex + 1, high); 
    } 
} 
  
// Function to print an array
void printArray(float sortN[], int size) 
{ 
    int i; 
    for (i=0; i < size; i++) 
        printf("%f ", sortN[i]); 
    printf("n"); 
}

int main()
{
    clock_t begin = clock();

    readDataFromFile();
    
    // variables needed for calculating equation for the regression line
    float sumX = 0.0, sumY = 0.0, xMean, yMean, slopeNumer = 0.0, slopeDenom = 0.0, b0, b1;

    // variables needed for calculating correlation coefficient & coefficient of determination
    float sumXY, sumXSq, sumYSq, sqSumX, sqSumY, correlationCoefficient, coefficientOfDetermination;
    
    // variables needed for calculating standard error
    float stdError, stdErrorNumer;

    for (int i = 0; i < N; i++)
    {
        float x = X_VALUES[i];
        float y = Y_VALUES[i];

        // calculate all summations needed
        sumX += x;
        sumY += y;
        sumXY += x*y;
        sumXSq += pow(x, 2);
        sumYSq += pow(y, 2);
    }

    // calculate square of the summation of x
    sqSumX += pow(sumX, 2);
    // calculate square of the summation of y
    sqSumY += pow(sumY, 2);

    //calculate the mean of x and y
    xMean = sumX / N;
    yMean = sumY / N;

    // Calculate correlation coefficient (r) & coefficient of determination (r squared)
    float rNumer = (N*sumXY)-(sumX*sumY);
    float rDenom = sqrt(((N*sumXSq)-(sqSumX))*((N*sumYSq)-sqSumY));
    correlationCoefficient = rNumer/rDenom;
    coefficientOfDetermination = pow(correlationCoefficient, 2) * 100;

    for (int i = 0; i < N; i++)
    {
        float xMinusXMean = X_VALUES[i] - xMean;
        float yMinusYMean = Y_VALUES[i] - yMean;

        // calculate the numerator & denominator of regression slope equation
        slopeNumer += xMinusXMean * yMinusYMean;
        slopeDenom += pow(xMinusXMean, 2);
    }

    // calculate standard error
    stdErrorNumer = slopeDenom; // the numerator for the standard error equation is the same as the denominator of regression slope equation
    stdError = sqrt(stdErrorNumer/(N-1));

    //calculate b0 (y-intercept) and b1 (slope) of regression line
    b1 = slopeNumer / slopeDenom;
    b0 = yMean - (b1*xMean);

    //calculate the n variable
    for (int i = 0; i < N; i++)
    {
        N_VALUES[i] = Y_VALUES[i] - (b1*X_VALUES[i]) -  b0;
        //printf("The n varianle is: %0.2f \n", N_VALUES[i]);
    }

    quickSort(N_VALUES, 0, N-1); 
    printf("Sorted array: \n"); 
    printArray(N_VALUES, N); 

    // printf("y = %0.2f + %0.2fx\n", b0, b1);
    // printf("sum of x and y: %0.2f and %0.2f\n", sumX, sumY);
    // printf("Mean of x and y: %0.2f and %0.2f\n", xMean, yMean);
    // printf("The correlation coefficient is %0.2f\n", correlationCoefficient);
    // printf("The coefficient of determination is %0.2f%%\n", coefficientOfDetermination);
    // printf("The standard error is %f\n", stdError);

    //plotGraph(b1, b0);
    
    //calculate the time taken for the the program to run.
    clock_t end = clock();
    double progTime = (double) (end-begin)/CLOCKS_PER_SEC;
    printf("The program took: %lf seconds \n", progTime);

    return 0;
}

void readDataFromFile()
{
    char x[1000], *y, countChar;
    FILE *fptr;
    int count = 0;

    //declare the delim
    const char delim[2] = ",";

    printf("\nPlease enter the file path for your dataset: ");
    scanf("%s", &DATASET_FILEPATH);
    fptr = fopen(DATASET_FILEPATH, "r");

    // file path of Mag's desktop (file path for testing)- /Users/magdalene/Desktop/SIT-UofG/programMeth/progMeth_assg1/Group9_15.txt
    if (fptr == NULL)
    {
        printf("ERROR! Please try again. Unable to read the following file from path: %s", DATASET_FILEPATH);
        // Program exits if file pointer returns NULL.
        exit(1);
    }

    //loop while countChar is not at the End of File
    while (countChar != EOF)
    {
        //Count whenever new line is encountered
        if (countChar == '\n')
        {
            count = count + 1;
        }

        //get the x coordinates and store in x
        fscanf(fptr,"%10000s[^\n]", x);

        //get the first token (getting the y coordinates)
        y = strtok(x, delim);
        
        // walk through other tokens while the token is null
        while( y != NULL ) {
            //printf( " %s\n", token );
            Y_VALUES[count] = (float)atof(y); //store token to the Y_VALUES for every count
            y = strtok(NULL, delim); //make the token null after storing
        }

        //store the x coorindates in array X_VALUES
        X_VALUES[count] = (float)atof(x);

        //take next character from file.
        countChar = getc(fptr);
    }

    fclose(fptr); //close file.
    //end of pulling coordinates from file
}

void plotGraph(float slope, float yIntercept)
{
    // Plot graph
    FILE *gp;
    gp = popen(GNUPLOT, "w"); // pipe to gnuplot program
    if (gp == NULL) {
        printf("Error opening pipe to GNU plot.\n"
            "Install with 'sudo apt-get install gnuplot' or 'brew install gnuplot'.\n");
        exit(0);
    }

    fprintf(gp, "set datafile separator comma\n");
    fprintf(gp, "f(x) = m*x + b\n");
    fprintf(gp, "set fit quiet\n"); // disables automatic output values from GNUPlot
    fprintf(gp, "fit f(x) '%s' using 1:2 via m, b\n", DATASET_FILEPATH);
    fprintf(gp, "plot '%s', f(x) title 'Regression Line y=%0.2fx+%0.2f'\n", DATASET_FILEPATH, slope, yIntercept);
    fclose(gp);

}