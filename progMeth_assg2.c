#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define GNUPLOT "/usr/local/Cellar/gnuplot/5.2.7_1/bin/gnuplot -persist" // change this based on where your gnuexecutable file is located. To locate it, type `brew ls gnu in the terminal
#define N 10000 // The amount of values in the dataset

char DATASET_FILEPATH[1000]; // "/Users/magdalene/Desktop/SIT-UofG/programMeth/progMeth_assg1/Group9_15.txt";
char DATAWRITE_FILEPATH [1000] = "N_VALUES.txt";
float X_VALUES[N], Y_VALUES[N], Y2_VALUES[N], N_VALUES[N];

// variables needed for calculating correlation coefficient & coefficient of determination
float sumXY, sumXSq, sumYSq, sqSumX, sqSumY, correlationCoefficient, coefficientOfDetermination;

void readDataFromFile(); 
void storeDataToFile ();

void plotGraph(float slope, float yIntercept, float n_mean, float standardErr);

void swap(float* a, float* b);
void quickSort(float sortValues[], int low, int high);
void printArray(float sortValues[], int size);

struct timer
{
    clock_t begin, end;
    double timeTaken;
};

int main()
{
    struct timer timeProg;
    timeProg.begin = clock();

    readDataFromFile();
    
    // variables needed for calculating equation for the regression line
    float sumX = 0.0, sumY = 0.0, xMean, yMean, slopeNumer = 0.0, slopeDenom = 0.0, b0, b1;
    
    // variables needed for calculating standard error
    float stdError, stdErrorNumer;

    // variables needed for calculating Y2_VALUES
    float sumN, nMean;

    float deviationNumer, deviationDenom, varN, standardDeviation;

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

    //calculate b0 (y-intercept) and b1 (slope) of regression line
    b1 = slopeNumer / slopeDenom;
    b0 = yMean - (b1*xMean);

    for (int i = 0; i < N; i++)
    {
        //calculation of standard error
        Y2_VALUES[i] = (b1*X_VALUES[i]) + b0;

        float YCorruptMinusYActual = Y_VALUES[i] - Y2_VALUES[i];
        stdErrorNumer += pow(YCorruptMinusYActual, 2); 

        //calculation of N_VALUES
        float n;
        N_VALUES[i] = Y_VALUES[i] - (b1*X_VALUES[i]) -  b0;
        n = N_VALUES[i];
        sumN += n;
    }

    // calculate standard error
    stdError = sqrt(stdErrorNumer/(N-2));

    nMean = sumN / N; // mean of all the N_VALUES

    //calculation of standard deviation
    for (int i = 0; i < N; i++)
    {
        float nMinusNmean = N_VALUES[i] - nMean;
        deviationNumer += pow(nMinusNmean, 2); 
    }

    varN = deviationNumer / (N-1);
    standardDeviation = sqrt(varN);

    //sort N_VALUES, X_VALUES and Y_VALUES in ascending order using quicksort
    quickSort(N_VALUES, 0, N-1);
    storeDataToFile();

    printf("y = %fx + %f\n", b1, b0);
    printf("sum of x and y: %0.2f and %0.2f\n", sumX, sumY);
    printf("Mean of x and y: %0.2f and %0.2f\n", xMean, yMean);
    printf("The correlation coefficient is %0.2f\n", correlationCoefficient);
    printf("The coefficient of determination is %0.2f%%\n", coefficientOfDetermination);
    printf("The standard error is %f\n", stdError);
    printf("The standard deviation based on n values is: %f \n", standardDeviation);
    printf("The mean of all the n values is %f\n", nMean);
    //printArray(N_VALUES, N);

    plotGraph(b1, b0, nMean, stdError);

    timeProg.end = clock();
    timeProg.timeTaken = (double) (timeProg.end-timeProg.begin) / CLOCKS_PER_SEC;
    printf("The program took: %lf seconds \n", timeProg.timeTaken);

    return 0;
}

void readDataFromFile()
{
    char x[1000], *y, countChar;
    FILE *fptr;
    int count = 0;

    //declare the delim
    const char delim[2] = ",";

    printf("\nPlease enter the file path to your dataset: ");
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

void storeDataToFile ()
{
    FILE *fptr;
    fptr = fopen(DATAWRITE_FILEPATH, "w");

    // file path of Mag's desktop (file path for testing)- /Users/magdalene/Desktop/SIT-UofG/programMeth/progMeth_assg1/Group9_15.txt
    if (fptr == NULL)
    {
        printf("ERROR! Please try again. Unable to read the following file from path: %s", DATAWRITE_FILEPATH);
        // Program exits if file pointer returns NULL.
        exit(1);
    }

    for (int i = 0; i < N; i++)
    {
        fprintf(fptr, "%f\n", N_VALUES[i]);
    }

    fclose(fptr);
}

void plotGraph(float slope, float yIntercept, float n_mean, float standardErr)
{
    // graph
   FILE *fptr;
    fptr = popen(GNUPLOT, "w"); // pipe to gnuplot program
    if (fptr == NULL) {
        printf("Error opening pipe to GNU \n"
            "Install with 'sudo apt-get install gnu or 'brew install gnu.\n");
        exit(0);
    }

    fprintf(fptr, "set multiplot layout 2,1 columnsfirst\n"); // set GNUPLOT to display 2 graphs at once, in one column
    
    //plot all 10 0000 points and regression line
    fprintf(fptr, "set datafile separator comma\n");
    fprintf(fptr, "f(x) = m*x + b\n");
    fprintf(fptr, "set fit quiet\n"); // disables automatic output values from GNUPlot
    fprintf(fptr, "fit f(x) '%s' using 1:2 via m, b\n", DATASET_FILEPATH);
    fprintf(fptr, "set key title 'Linear Regression Graph'\n");
    fprintf(fptr, "plot '%s', f(x) title 'Regression Line y=%0.2fx+%0.2f'\n", DATASET_FILEPATH, slope, yIntercept);
    
    //plot histogram
    fprintf(fptr, "binwidth=0.5\n");
    fprintf(fptr, "set boxwidth binwidth\n");
    fprintf(fptr, "bin(x,width)=width*floor(x/width) + binwidth/2.0\n");

    //plot curve
    fprintf(fptr, "set arrow from '%f', graph 0 to '%f', graph 1 nohead lw 2 lc rgb 'red'\n", n_mean, n_mean);
    fprintf(fptr, "set arrow from '%f', graph 0 to '%f', graph 1 nohead lw 2 lc rgb 'red'\n", (n_mean-standardErr), (n_mean-standardErr));
    fprintf(fptr, "set arrow from '%f', graph 0 to '%f', graph 1 nohead lw 2 lc rgb 'red'\n", (n_mean+standardErr), (n_mean+standardErr));
    fprintf(fptr, "Gauss(x,mu,sigma) = 1./(sigma*sqrt(2*pi)) * exp( -(x-mu)**2 / (2*sigma**2) )\n");
    fprintf(fptr, "d1(x) = Gauss(x,'%f','%f')*binwidth*'%d'\n", n_mean, standardErr,N); // multiply by N and binwidth to scale curve to histogram.
    fprintf(fptr, "set key title 'Histogram'\n");
    fprintf(fptr, "plot d1(x), '%s' using (bin($1,binwidth)):(1.0) smooth freq with boxes lt rgb '#1b9646', '%s' using (bin($1,binwidth)):(1.0) smooth freq with line lt rgb '#000080'\n", DATAWRITE_FILEPATH, DATAWRITE_FILEPATH);
    fclose(fptr);
}

// functions for to quicksort values
//swaps 2 elements
void swap(float* a, float* b) 
{ 
    float t = *a; 
    *a = *b; 
    *b = t; 
} 

//uses last element as the pivot
//elements < pivot value will be placed on the left of pivot
//elements > pivot value will be placed on the right of pivot
int partition (float sortValues[], int low, int high) 
{ 
    float pivot = sortValues[high];    // pivot 
    int smallerIndex = (low - 1);  // Index of smaller element 
  
    for (int j = low; j <= high- 1; j++) 
    { 
        //if current element is smaller than the pivot
        if (sortValues[j] < pivot) 
        { 
            smallerIndex++;    //increment index of smaller element by one
            swap(&sortValues[smallerIndex], &sortValues[j]); 
        } 
    } 
    swap(&sortValues[smallerIndex + 1], &sortValues[high]); 
    return (smallerIndex + 1); 
} 
  
//sort N_VALUES using quicksort
//low represents starting index and high is the ending index
void quickSort(float sortValues[], int low, int high) 
{ 
    if (low < high) 
    { 
        //puts sortValues[] in the right position
        int partitioningIndex = partition(sortValues, low, high); 
  
        // Sort the elements before and after partition, seperatelty
        quickSort(sortValues, low, partitioningIndex - 1); 
        quickSort(sortValues, partitioningIndex + 1, high); 
    } 
}
