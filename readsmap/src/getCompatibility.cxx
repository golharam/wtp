#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "util.h"
#include "getComp.cxx"

void set_addColors(int addColors[4][4])
{
  addColors[0][0] = 0;
  addColors[1][0] = 1;
  addColors[0][1] = 1;
  addColors[1][1] = 1;
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
}

/*Given an array of integers for two sequences and the length of them
  Returns 1 if the color sum of the two sequences are the same,
  0 otherwise.
*/
int is_compatible(int *seq1, int seq1_length, int *seq2, int seq2_length)
{
  int addColors[4][4];
  set_addColors(&addColors[0]);
  int i,j;

  /*
  for(i=0; i<4; i++) {
    for(j=0; j<4; j++) {
      printf("%d %d %d\n", i, j, addColors[i][j]);
    }
  }
  */

  /*
  const int seq1_length=8;
  const int seq2_length=2;
  int seq1[seq1_length]={0,2,2,1,3,2,1,2};
  int seq2[seq2_length]={3,2};
  */

  int seq1sum = seq1[0];
  int seq2sum = seq2[0];

  if (seq1_length>1) {
    for(i=1;i<seq1_length;i++) {
      printf("addColors[%d][%d]=", seq1sum, seq1[i]);
      seq1sum=addColors[seq1sum][seq1[i]];
      printf("%d\n", seq1sum);
    }
  }
  if (seq2_length>1) {
    for(j=1;j<seq2_length;j++) {
      seq2sum=addColors[seq2sum][seq2[j]];
    }
  }
  
  printf("seq1sum=%d, seq2sum=%d\n", seq1sum, seq2sum);
  if(seq1sum==seq2sum)
    return TRUE;
  else
    return FALSE;
}

/*Given an array of integers for a sequence and the length of it,
  returns the color sum.
  To see if an indel is valid, take the sequence of the deleted region
  and see if it sums to 0.  This assumes that there are no mismatches
  at the ends.
*/
int get_sum(int *seq, int seq_length) //, int *seq2, int seq2_length)
{
  int addColors[4][4];
  set_addColors(&addColors[0]);
  int i,j;

  /*
  const int seq_length=8;
  int seq[seq_length]={0,2,2,1,3,2,1,2};
  */

  int seqSum = seq[0];

  if (seq_length>1) {
    for(i=1;i<seq_length;i++) {
      printf("addColors[%d][%d]=", seqSum, seq[i]);
      seqSum=addColors[seqSum][seq[i]];
      printf("%d\n", seqSum);
    }
  }
  
  printf("seqSum=%d\n", seqSum);
  return seqSum;
}


void main_old(int argc, char *argv[])
{
  {
    //test code 1
    printf("test1:\n");
    int i,j;
    const int seq1_length=8;
    const int seq2_length=3;
    int seq1[seq1_length]={1,2,2,1,3,2,1,2};
    int seq2[seq2_length]={1,1,2};
    
    printf("seq1= (");
    for(i=0;i<seq1_length;i++) {
      printf("%d ", seq1[i]);
    }
    printf(")\nseq2= (");
    for(j=0;j<seq2_length;j++) {
      printf("%d ", seq2[j]);
    }
    printf(")\n");
    printf("is_compatible = %d\n",
	   is_compatible(seq1,seq1_length,seq2,seq2_length)
	   );
  }
  {
    //test code 2
    printf("\ntest2:\n");
    int i,j;
    const int seq1_length=8;
    const int seq2_length=2;
    int seq1[seq1_length]={1,3,0,2,3,0,1,2};
    int seq2[seq2_length]={1,2};
    
    printf("seq1= (");
    for(i=0;i<seq1_length;i++) {
      printf("%d ", seq1[i]);
    }
    printf(")\nseq2= (");
    for(j=0;j<seq2_length;j++) {
      printf("%d ", seq2[j]);
    }
    printf(")\n");
    printf("is_compatible = %d\n",
	   is_compatible(seq1,seq1_length,seq2,seq2_length)
	   );
  }

  {
    //test code 3
    printf("\ntest3:  (doing sume only)\n");
    int i,j;
    const int seq1_length=8;
    const int seq2_length=4;

    int seq1[seq1_length]={1,2,2,1,3,2,1,2};
    int seq2[seq2_length]={3,0,2,1};
    
    printf("seq1= (");
    for(i=0;i<seq1_length;i++) {
      printf("%d ", seq1[i]);
    }
    printf(")\n");
    printf("sum=%d\n", get_sum(seq1,seq1_length));

    printf("seq2= (");
    for(j=0;j<seq2_length;j++) {
      printf("%d ", seq2[j]);
    }
    printf(")\n");
    printf("sum=%d\n", get_sum(seq2,seq2_length));
  }


}

int main(int argc, char *argv[]) {

  int i,j;
  getComp comparison;
  char* seq1 = "13023032";
  int seq1_length=8;

  char* seq2 = "13";
  int seq2_length=2;

  printf ("test1\n");
  printf("seq1= %s\tsum= %d\n", seq1, comparison.get_sum(seq1,seq1_length));
  printf("seq2= %s\tsum= %d\n", seq2, comparison.get_sum(seq2,seq2_length));
  printf("is_compatible = %d\n",  
	 comparison.is_compatible(seq1,seq1_length,seq2,seq2_length));

  printf ("test2\n");
  seq1="30232112";
  seq2="21";
  printf("seq1= %s\tsum= %d\n", seq1, comparison.get_sum(seq1,seq1_length));
  printf("seq2= %s\tsum= %d\n", seq2, comparison.get_sum(seq2,seq2_length));
  printf("is_compatible = %d\n",  
	 comparison.is_compatible(seq1,seq1_length,seq2,seq2_length));
}
