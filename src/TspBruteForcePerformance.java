import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.io.*;
import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Random;
import java.util.function.Function;

public class TspBruteForcePerformance {

    static ThreadMXBean bean = ManagementFactory.getThreadMXBean();
    static int permCount = 0; // permutation counter
    static int maxPerms = 0; // maximum number of permutations
    static double shortestLength = 2147483647; // max value of int
    static int[] shortestPath;
    static double[][] costMatrix;
    static double[] xCoord;
    static double[] yCoord;
    static double[] expectedPathLength;

    /* define constants */
    static long MAXVALUE = 2000000000;
    static long MINVALUE = -2000000000;
    static int numberOfTrials = 1;                   // adjust numberOfTrials and MAXINPUTSIZE based on available
    static int MAXINPUTSIZE = 14;  // time, processor speed, and available memory
    static int MININPUTSIZE = 1;

    static String ResultsFolderPath = "/home/steve/Results/"; // pathname to results folder
    static FileWriter resultsFile;
    static PrintWriter resultsWriter;

    public static void main(String[] args) {
        // run the whole experiment at least twice, and expect to throw away the data from the earlier runs, before java has fully optimized
        System.out.println("Running first full experiment...");
        runFullExperiment("TspBruteForce-Exp1-ThrowAway.txt");
        System.out.println("Running second full experiment...");
        runFullExperiment("TspBruteForce-Exp2.txt");
        System.out.println("Running third full experiment...");
        runFullExperiment("TspBruteForce-Exp3.txt");

        int numVertices = 5;
        int[] verticesArray = new int[numVertices];
        for (int i = 0; i < numVertices; i++)
            verticesArray[i] = i;
        costMatrix = generateRandomCircularGraphCostMatrix(numVertices, 100);
        maxPerms = factorial(numVertices-1);
        TspBruteForcePerformance permutation = new TspBruteForcePerformance();
        permutation.permute(verticesArray, 0, numVertices-1);

        // append 0 to shortestPath
        int appendedPath[] = new int[shortestPath.length + 1];
        for (int i = 0; i < numVertices; i++)
            appendedPath[i] = shortestPath[i];
        appendedPath[numVertices] = 0;

        System.out.println("\nCost Matrix:");
        print2dDoubleArrayPretty(costMatrix);
        System.out.println("\nTesting using circular matrix:");
        System.out.println("The shortest path is: " + Arrays.toString(appendedPath));
        System.out.printf("The shortest length is: %.2f", shortestLength);
        System.out.println("\n");
        System.out.print("Vertex Coords:  ");
        printIntArrayPretty(verticesArray);
        System.out.print("The x coords are: ");
        printDoubleArrayPretty(xCoord);
        System.out.print("The y coords are: ");
        printDoubleArrayPretty(yCoord);
    }

    static void runFullExperiment(String resultsFileName) {
        try {
            resultsFile = new FileWriter(ResultsFolderPath + resultsFileName);
            resultsWriter = new PrintWriter(resultsFile);
        } catch (Exception e) {
            System.out.println("*****!!!!!  Had a problem opening the results file " + ResultsFolderPath + resultsFileName);
            return; // not very foolproof... but we do expect to be able to create/open the file...
        }

        ThreadCpuStopWatch BatchStopwatch = new ThreadCpuStopWatch(); // for timing an entire set of trials
        ThreadCpuStopWatch TrialStopwatch = new ThreadCpuStopWatch(); // for timing an individual trial

        resultsWriter.println("#InputSize    AverageTime"); // # marks a comment in gnuplot data
        resultsWriter.flush();

        /* for each size of input we want to test: in this case starting small and doubling the size each time */
        for (int inputSize = MININPUTSIZE; inputSize <= MAXINPUTSIZE; inputSize++) {
            // progress message...
            System.out.println("Running test for input size " + inputSize + " ... ");

            /* repeat for desired number of trials (for a specific size of input)... */
            long batchElapsedTime = 0;
            // generate a list of random integers in random order to use as test input
            // In this case we're generating one list to use for the entire set of trials (of a given input size)
            //System.out.print("    Generating test data...");
            //long[] testList = createRandomIntegerList(inputSize);
            //System.out.println("...done.");
            //System.out.print("    Running trial batch...");

            /* force garbage collection before each batch of trials run so it is not included in the time */
            System.gc();

            // instead of timing each individual trial, we will time the entire set of trials (for a given input size)
            // and divide by the number of trials -- this reduces the impact of the amount of time it takes to call the
            // stopWatch methods themselves
            //BatchStopwatch.start(); // comment this line if timing trials individually

            // run the trials
            for (long trial = 0; trial < numberOfTrials; trial++) {
                // generate a random list of integers each trial
                double[][] testMatrix = generateRandomEuclideanCostMatrix(inputSize, 100, 100);
                int[] testArray = new int[testMatrix[0].length];

                // generate a random key to search in the range of a the min/max numbers in the list
                // long testSearchKey = (long) (0 + Math.random() * (testList[testList.length - 1]));
                /* force garbage collection before each trial run so it is not included in the time */
                // System.gc();

                TrialStopwatch.start(); // *** uncomment this line if timing trials individually
                /* run the function we're testing on the trial input */
                TspBruteForcePerformance expPermutation = new TspBruteForcePerformance();
                expPermutation.permute(testArray, 0, inputSize-1);
                batchElapsedTime = batchElapsedTime + TrialStopwatch.elapsedTime(); // *** uncomment this line if timing trials individually
            }
            //batchElapsedTime = BatchStopwatch.elapsedTime(); // *** comment this line if timing trials individually
            double averageTimePerTrialInBatch = (double) batchElapsedTime / (double) numberOfTrials; // calculate the average time per trial in this batch

            /* print data for this size of input */
            resultsWriter.printf("%12d  %15.2f \n", inputSize, averageTimePerTrialInBatch); // might as well make the columns look nice
            resultsWriter.flush();
            System.out.println(" ....done.");
        }
    }

    // Gets all permutations of an array and calculates best path via brute force algorithm.
    // I modified and adapted the code from:
    // https://www.geeksforgeeks.org/write-a-c-program-to-print-all-permutations-of-a-given-string/
    private void permute(int[] arr, int startIndex, int endIndex)
    {
        if ((startIndex == endIndex) && (permCount < maxPerms)) {
            double pathLength = costMatrix[arr[arr.length-1]][arr[0]];
            for (int i = 0; i < arr.length-1; i++){
                pathLength = pathLength + costMatrix[arr[i]][arr[i+1]];
            }
            //System.out.println(pathLength);
            if (pathLength < shortestLength){
                shortestLength = pathLength;
                shortestPath = arr.clone();
            }
            permCount++;
        }
        else
        {
            for (int i = startIndex; i <= endIndex; i++)
            {
                arr = swap(arr,startIndex,i);
                permute(arr, startIndex+1, endIndex);
                arr = swap(arr,startIndex,i);
            }
        }
    }

    // swap i and j
    public int[] swap(int[] arr, int i, int j){
        int temp = arr[i] ;
        arr[i] = arr[j];
        arr[j] = temp;
        return arr;
    }

    // Method to find factorial of given number
    static int factorial(int n){
        int res = 1, i;
        for (i=2; i<=n; i++)
            res *= i;
        return res;
    }

    public static double[][] generateRandomCostMatrix(int numVertices, double maxEdgeLength) {
        double[][] costMatrix = new double[numVertices][numVertices];
        for (int i = 0; i < numVertices; i++) {
            costMatrix[i][i] = 0;
            for (int j = i+1; j < numVertices; j++) {
                costMatrix[i][j] = (Math.random() * ((maxEdgeLength - 1) + 1)) + 1;
                costMatrix[j][i] = costMatrix[i][j];
            }
        }
        return costMatrix;
    }

    public static double[][] generateRandomEuclideanCostMatrix(int numVertices, int xMax, int yMax){
        double [][] costMatrix = new double[numVertices][numVertices];
        int[] x = new int[numVertices];
        int[] y = new int[numVertices];
        // generate random coordinates
        for (int i = 0; i < numVertices; i++){
            x[i] = (int) (Math.random() * ((xMax - 1) + 1)) + 1;
            y[i] = (int) (Math.random() * ((yMax - 1) + 1)) + 1;
        }
        // build costMatrix
        for (int i = 0; i < numVertices; i++){
            costMatrix[i][i] = 0;
            for (int j = i+1; j < numVertices; j++){
                costMatrix[i][j] = Math.sqrt((Math.pow((x[i]-x[j]),2) + Math.pow((y[i]-y[j]),2)));
                costMatrix[j][i] = costMatrix[i][j];
            }
        }
        return costMatrix;
    }

    public static double[][] generateRandomCircularGraphCostMatrix(int numVertices, int radius){
        double [][] costMatrix = new double[numVertices][numVertices];
        double[] x = new double[numVertices];
        double[] y = new double[numVertices];
        double stepAngle = 2*Math.PI/numVertices;

        // generate circular coordinates
        for (int i = 0; i < numVertices; i++){
            x[i] = radius * Math.sin(i * stepAngle);
            y[i] = radius * Math.cos(i * stepAngle);
        }
        shuffleArray(x, y);

        // build costMatrix
        for (int i = 0; i < numVertices; i++){
            costMatrix[i][i] = 0;
            for (int j = i+1; j < numVertices; j++){
                costMatrix[i][j] = Math.sqrt((Math.pow((x[i]-x[j]),2) + Math.pow((y[i]-y[j]),2)));
                costMatrix[j][i] = costMatrix[i][j];
            }
        }

        // for verification testing
        xCoord = x.clone();
        yCoord = y.clone();
        //expectedPathLength =

        return costMatrix;
    }

    // shuffle array
    // code modified from: https://www.vogella.com/tutorials/JavaAlgorithmsShuffle/article.html
    public static void shuffleArray(double[] x, double[] y) {
        int n = x.length;
        Random random = new Random();
        random.nextInt();
        for (int i = 0; i < n; i++) {
            int change = i + random.nextInt(n - i);
            swapDoubles(x, i, change);
            swapDoubles(y, i, change);
        }
    }
    // swap i and j
    public static double[] swapDoubles(double[] arr, int i, int j)
    {
        double temp = arr[i] ;
        arr[i] = arr[j];
        arr[j] = temp;
        return arr;
    }

    public static void print2dIntArray(int[][] mat){
        for (int[] row : mat)
            System.out.println(Arrays.toString(row));
    }

    public static void print2dDoubleArray(double[][] arr){
        for (double[] row : arr)
            System.out.println(Arrays.toString(row));
    }

    public static void print2dDoubleArrayPretty(double[][] arr){
        for (double[] row : arr) {
            for (double x : row)
                System.out.printf("%.2f\t", x);
            System.out.println();
        }
    }

    public static void printDoubleArrayPretty(double[] arr){
        for (int i = 0; i < arr.length; i++)
                System.out.printf("%8.2f", arr[i]);
        System.out.println();
    }

    public static void printIntArrayPretty(int[] arr){
        for (int i = 0; i < arr.length; i++)
            System.out.printf("%8d", arr[i]);
        System.out.println();
    }
}