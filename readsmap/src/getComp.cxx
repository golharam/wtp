#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "getComp.h"

getComp::getComp()
{
  int i;
  for (i = 0; i < 128; i++) {
	 addColors[0][i] =  addColors[1][i] = addColors[2][i] = addColors[3][i] = addColors[4][i] =  4;
  }
  addColors[0][0] = 0;
  addColors[1][0] = 1;
  addColors[0][1] = 1;
  addColors[1][1] = 0;  //correction!
  addColors[2][0] = 2;
  addColors[0][2] = 2;
  addColors[2][1] = 3;
  addColors[1][2] = 3;
  addColors[2][2] = 0;
  addColors[3][0] = 3;
  addColors[0][3] = 3;
  addColors[3][1] = 2;
  addColors[1][3] = 2;
  addColors[3][2] = 1;
  addColors[2][3] = 1;
  addColors[3][3] = 0;
  char C[4][4] = {"0Aa", "1Cc", "2Gg", "3Tt"}; 
  for (i = 0; i < 4; i++) {
    int j, k;
    for (j = 0; j < 4; j++) {
	for (k = 0; k < 3; k++) {
	    addColors[i][C[j][k]] = addColors[i][j];
	}
    }
  } 
}

/*Given an array of integers for two sequences and the length of them
  Returns 1 if the color sum of the two sequences are the same,
  0 otherwise.
*/
int getComp::is_compatible(char *seq1, int seq1_length, char *seq2, int seq2_length)
{
    return (get_sum(seq1, seq1_length) == get_sum(seq2, seq2_length));
}

/*Given an array of integers for a sequence and the length of it,
  returns the color sum.
  To see if an indel is valid, take the sequence of the deleted region
  and see if it sums to 0.  This assumes that there are no mismatches
  at the ends.
  seq can be using values 0-3 or ascii number 0-3 or letter ACGT(acgt) to code 0-3
*/
int getComp::get_sum(char *seq, int seq_length)
{
  int i,j;

  int seqSum = addColors[0][seq[0]];

  if (seq_length>1) {
    for(i=1;i<seq_length;i++) {
      //printf("addColors[%d][%d]=", seqSum, seq[i]);
      seqSum=addColors[seqSum][seq[i]];
      //printf("%d\n", seqSum);
    }
  }
  
  //printf("seqSum=%d\n", seqSum);
  return seqSum;
}


