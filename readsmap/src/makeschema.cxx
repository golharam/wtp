/* This program produces schemas allowing adj mismatches being counted as 1*/*/

#include "makeschema.h"
// a list of num blocks, each size of blocks[i], choose nmis of them to put 0's 
// rest in 1s.

static void fillstring_noend(char *s, int n, char a)
{
    for (i = 0; i < n; i++) s[i] = a;
}

static void fillstring(char *s, int n, char a)
{
    for (i = 0; i < n; i++) s[i] = a;
    s[i] = 0;
}


// len choose n.
int allcomb(char **schema, int *blocks, int len, int n, int memsize)
{
    int num[len], statr[len],j, k, sn = 0, rlen = 0;
    for (i = 0, start[0] = 0; i < n; i++) {
	num[i] = i;
	rlen += blocks[i];
	start[i+1] = start[i]+blocks[i];
    }
    do {
	if (sn >= memsize) return outmem();
	fillstring(schema[sn], rlen, '1');
            for (i = 0; i < n; i++) {
		if (num[i] < len-1) {
		    fillstring_noend(schema[sn]+ start[num[i]], blocks[num[i]]+1, '0');
		} else {
                    fillstring_noend(schema[sn]+ start[num[i]], blocks[num[i]], '0');
		}
            }

        for (j = n-1; j>=0; j--) {
            if (num[j] < len-1- (n-1-j)) {
                num[j]++;
                int k;
                for (k = j+1; k < n; k++)
                    num[k] = num[k-1]+1;
                break;
            }
        }
	sn++;
    } while (j >= 0);
    return sn;
}

int merge_2schemas(char **schema, int n1, int n2, int len1, int len2) 
{
    int i;

    for (i = 0; i < n1; i++) {
	fillstring(schema[i]+len1, len2+1, '0');
    }
    for (i = 0; i < n1+n2; i++) {
	for (j = 0; j < len2; j++) 
	    schema[i][j+len1+1] = schema[i][j];
	fillstring_noend(schema[i], len1+1, '1'); 
    }
    return n1+n2;
}

int outmem()
{
    fprintf(stderr, "not enough memory for schema array\n");
    return 0;
}

int makeschema(char **schema, int rlen, int m, int memsize)
{
    if (rlen < 10*m) {
	fprintf(stderr, "Current do not support rlen<10m\n");
	return 0;
    }
    if (rlen < 14) return 0;
    if (m == 0) {
	int len = 14;
	fillstring(schema[0], len, '1');
	return 1;
    }
    if (m == 1) {
	if (len < 18) return 0;
	int x = len-14;
	int a = 14, b = 0;
	for (i = 0; 1; i++) {
	    if (i >= memsize) return outmem(); 
	    fillstring(schema[i], a, '1');
	    fillstring(schema[i]+a, x, '0');
	    fillstring(schema[i]+a+b, '1');
	    if (a == 0) break;
	    b += x-1;
	    a -= x-1;
	    if (a < 0) {a = 0; b = 14;}
	}
	return i;
    }
    if (rlen >35) {
	int s1, s2, X, Y;
	if (m == 2) {
	    s1 = makeschema(schema, (X=rlen-15), 1, memsize);
	    s2 = 1;
	    fillstring(schema[s1], (Y=14), '1');
	} else {
	    int m1 = (m-1)/2;
	    int m2 = m-1-m1;
	    X = (rlen-1)*m1/(m1+m2);
	    Y = rlen-1-X;
	    s1 = makeschema(schema, X, m1, memsize);
	    s2 = makeschema(schema+s1, Y, m2, memsize-s1);
	}
	return merge_2schemas(schema, s1, s2, X, Y);
    }
    

}

