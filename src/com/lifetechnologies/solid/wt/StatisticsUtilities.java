package com.lifetechnologies.solid.wt;

import java.util.*;

/**
 * User: btuch
 * Date: Aug 19, 2004
 * Revision: $Rev$
 */
public class StatisticsUtilities {

    public static HashMap<Integer, Double> factorialLookupTable = new HashMap<Integer, Double>();

    public static double calculateMedian(Double[] sortedDoubles) {
        int middle = sortedDoubles.length/2;
        if (sortedDoubles.length%2 == 1) {
            // Odd number of elements -- return the middle one.
            return sortedDoubles[middle];
        } else {
           // Even number -- return average of middle two
           return (sortedDoubles[middle-1] + sortedDoubles[middle]) / 2.0;
        }
    }
    
    public static double calculateOverlapPValue(String countString1, String countString2, String countStringDelimiter,
                                                int minNumberOfOccurences1, int minNumberOfOccurences2, boolean scrambleCounts) {

        String string1[] = countString1.split(countStringDelimiter);
        String string2[] = countString2.split(countStringDelimiter);
        if (scrambleCounts)
            string1 = scrambleStringArray(string1);

        double countsArray1[] = new double[string1.length];
        double countsArray2[] = new double[string2.length];
        for (int i = 0; i < countsArray1.length; i++) {
            countsArray1[i] = Double.parseDouble(string1[i] + "");
            countsArray2[i] = Double.parseDouble(string2[i] + "");
        }

        double pValue = calculateOverlapPValue(countsArray1, countsArray2, minNumberOfOccurences1, minNumberOfOccurences2);

        return pValue;
    }


    public static double calculateOverlapPValue(double[] countString1, double[] countString2,
                                                int minNumberOfOccurences1, int minNumberOfOccurences2) {
        double pValue = 0.0;
        int numberOfOccurrencesAboveThresholdInString1 = 0;
        int numberOfOccurrencesAboveThresholdInString2 = 0;
        int overlap = 0;
        int populationSize = countString1.length;
        for (int i = 0; i < countString1.length; i++) {
            if (countString1[i] >= minNumberOfOccurences1)
                numberOfOccurrencesAboveThresholdInString1++;
            if (countString2[i] >= minNumberOfOccurences2)
                numberOfOccurrencesAboveThresholdInString2++;
            if (countString1[i] >= minNumberOfOccurences1 &&
                    countString2[i] >= minNumberOfOccurences2)
                overlap++;
        }

        pValue = calculateOverlapPValue(overlap, numberOfOccurrencesAboveThresholdInString1, numberOfOccurrencesAboveThresholdInString2, populationSize);

        return pValue;
    }

    public static double calculateOverlapPValue(int overlap, int numberOfOccurrences1, int numberOfOccurrences2, int populationSize) {
        double pValue = 0.0;
        for (int i = overlap; i <= Math.min(numberOfOccurrences1, numberOfOccurrences2); i++)
            pValue += Math.pow(Math.E, calculateLnHypergeometricProbability(populationSize, numberOfOccurrences1, numberOfOccurrences2, i));
        return pValue;
    }

    
    private static double calculateLnHypergeometricProbability(int sampleSize, int numberOfType1, int numberOfType2, int overlap) {
        return lnNChooseM(numberOfType1, overlap) + lnNChooseM(sampleSize - numberOfType1, numberOfType2 - overlap) - lnNChooseM(sampleSize, numberOfType2);

    }

    private static double lnNChooseM(int N, int M) {
        return lnFactorial(N) - lnFactorial(M) - lnFactorial(N - M);
    }

//    private static double lnFactorial2(int n) {
//        if (n < 2) {
//            return 0.0;
//        } else if (factorialLookupTable.containsKey(new Integer(n))) {
//            return ((Double)factorialLookupTable.get(new Integer(n))).doubleValue();
//        } else {
//            double lnFactorial = Math.log(n) + lnFactorial(n - 1);
//            factorialLookupTable.put(new Integer(n), new Double(lnFactorial));
//            return lnFactorial;
//        }
//    }

    private static double lnFactorial(int n) {

        int currentN = n;
        double lnFactorial = 0.0;
        while (currentN > 0 ) {
            if (factorialLookupTable.containsKey(new Integer(currentN))) {
                lnFactorial += ((Double)factorialLookupTable.get(new Integer(currentN))).doubleValue();
                currentN = 0;
            } else
                lnFactorial += Math.log(currentN--);
        }

        factorialLookupTable.put(new Integer(n), new Double(lnFactorial));

        return lnFactorial;
    }

    public static String scrambleString(String string) {

        String scrambledSequence = "";
        while (string.length() > 0) {
            int randomIndex = (int) (Math.random() * string.length());
            scrambledSequence += string.charAt(randomIndex);
            string = string.substring(0, randomIndex) + string.substring(randomIndex + 1);
        }

        return scrambledSequence;
    }

    public static Double[] scrambleDoubleArray(Double[] array) {

        Double scrambledArray[] = new Double[array.length];

        Vector<Double> arrayAsVector = new Vector<Double>();
        for (int i = 0; i < array.length; i++)
            arrayAsVector.add(new Double(array[i]));

        for (int i = 0; i < scrambledArray.length; i++) {
            int randomIndex = (int) (Math.random() * arrayAsVector.size());
            scrambledArray[i] = ((Double) arrayAsVector.remove(randomIndex)).doubleValue();
        }

        return scrambledArray;
    }

//    private static int[] scrambleIntegerArray(int[] array) {
//        int scrambledArray[] = new int[array.length];
//
//        Vector arrayAsVector = new Vector();
//        for (int i = 0; i < array.length; i++)
//            arrayAsVector.add(new Integer(array[i]));
//
//        for (int i = 0; i < scrambledArray.length; i++) {
//            int randomIndex = (int) (Math.random() * arrayAsVector.size());
//            scrambledArray[i] = ((Integer) arrayAsVector.remove(randomIndex)).intValue();
//        }
//
//        return scrambledArray;
//    }

    public static String[] scrambleStringArray(String[] array) {

        String scrambledArray[] = new String[array.length];

        Vector<String> arrayAsVector = new Vector<String>();
        for (int i = 0; i < array.length; i++)
            arrayAsVector.add(new String(array[i]));

        for (int i = 0; i < scrambledArray.length; i++) {
            int randomIndex = (int) (Math.random() * arrayAsVector.size());
            scrambledArray[i] = ((String) arrayAsVector.remove(randomIndex));
        }

        return scrambledArray;
    }

    public static double calculateMean(Double[] data) {
        double mean = 0.0;
        int dataPoints = 0;
        for (int i = 0; i < data.length; i++)
            if (!Double.isNaN(data[i])) {
                mean += data[i];
                dataPoints++;
            }

        mean /= dataPoints;

        return mean;
    }

    public static double calculateVariance(Double[] data) {
        double variance = 0.0;
        double mean = calculateMean(data);
        for (int i = 0; i < data.length; i++)
            if (!Double.isNaN(data[i])) {
                variance += (data[i] - mean) * (data[i] - mean);
            }

        variance /= data.length - 1;
        return variance;
    }

    public static double calculateVarianceOfPopulation(Double[] data) {
        double variance = 0.0;
        double mean = calculateMean(data);
        for (int i = 0; i < data.length; i++)
            if (!Double.isNaN(data[i])) {
                variance += (data[i] - mean) * (data[i] - mean);
            }

        variance /= data.length;
        return variance;
    }



    // the following subroutine returns
    // Log_10 of the Sum of Poission distribution
    // Log_10[ Sum(m >= n) x^n/n!*exp(-x)
    // taken from Hao's motif_contrast.pl script
    public static double calculateLogSumPoisson(int count, double expectation) {

        if (count == 0) return 0.0;     

        // calculate the log of first term in the summation
        double first = count * Math.log(expectation) - expectation;
        for (int i = 1; i <= count; i++)
            first -= Math.log(1.0 * i);

        double factor = 1.0;
        double addition = 1.0;

        for (int i = count + 1; addition / factor > 1.0e-15; i++) {
            addition *= (expectation / i);
            factor += addition;
        }

        return (first + Math.log(factor)) / Math.log(10.0);
    }


    public static double calculatePearsonCorrelation(Double vector1[], Double vector2[]) {

        double correlation = 0.0;

        double counts1Mean = calculateMean(vector1);
        double counts2Mean = calculateMean(vector2);
        double counts1StdDev = Math.sqrt(calculateVariance(vector1));
        double counts2StdDev = Math.sqrt(calculateVariance(vector2));

        for (int i = 0; i < vector1.length; i++) {
            correlation += (vector1[i] - counts1Mean) * (vector2[i] - counts2Mean) / counts1StdDev / counts2StdDev;
        }

        return correlation / vector1.length;
    }

    public static double calculatePearsonCorrelation(Vector<Double> dataVector1, Vector<Double> dataVector2) {

        Double vector1[] = new Double[dataVector1.size()];
        Double vector2[] = new Double[dataVector2.size()];
        for (int i = 0; i < vector1.length; i++) {
            vector1[i] = (Double)dataVector1.get(i);
            vector2[i] = (Double)dataVector2.get(i);
        }

        return calculatePearsonCorrelation(vector1, vector2);
    }

    public static double calculatePearsonCorrelationSignificance(Double[] vector1, Double[] vector2) {

        double significance = 0.0;
        double pearsonCorrelation = calculatePearsonCorrelation(vector1, vector2);
        double countGreaterThanObserved = 0;
        for (int i = 0; i < 10000; i++) {
            Double scrambledVector1[] = scrambleDoubleArray(vector1);
            if (calculatePearsonCorrelation(scrambledVector1, vector2) >= pearsonCorrelation)
                countGreaterThanObserved++;
        }

        significance = countGreaterThanObserved / 10000;
        return significance;
    }

    public static double findMax(double[] data) {
        double max = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < data.length; i++)
            if (data[i] > max)  max = data[i];

        return max;
    }

    public static void shuffleCounts(int[] countsA, int[] countsB, double probabilityToA, int[] countsShuffledA, int[] countsShuffledB) throws Exception {

        if (countsA.length != countsB.length || countsB.length != countsShuffledA.length || countsShuffledA.length != countsShuffledB.length)
            throw new Exception("Counts arrrays must all have the same size.");

        for (int indexArray = 0; indexArray < countsA.length; indexArray++) {
            int sum = countsA[indexArray] + countsB[indexArray];
            for (int indexCount = 0; indexCount < sum; indexCount++) {
                if (Math.random() < probabilityToA)
                    countsShuffledA[indexArray]++;
                else
                    countsShuffledB[indexArray]++;
            }
        }

    }
}