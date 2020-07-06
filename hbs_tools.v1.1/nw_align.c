/* Global sequence alignment by dynamic programming.
 *
 * The program is ANSI C and should compile on any machine that has
 * a C compiler, with a command line like:
 *     gcc -o hbs_align hbs_align.c
 * To run the program:
 *     ./hbs_align
 *
*/

/* This program is Copyright (C) 2015-16, Ming-an Sun (mingansun@gmail.com)

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

*/

#define MATCH      3
#define MISMATCH  -3
#define INDEL     -5

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[])
{   
  
    if(argc != 2){
        printf("\nDescription:\n");
        printf("  This script aligns the read1 and read2 from HB-Seq using needleman-wunsch algorithm.\n");
        printf("  Mismatches probably due to bisulfite conversion (T-C, G-A) are considered as match.\n");
        printf("\nUsage:\n");
        printf("  %s <file with name, seq1 and seq2>\n\n", argv[0]);
        printf("Version:\n");
        printf("  Ver 1.1, last modified on August 25, 2015\n\n");
        exit(1);
    }
    
	char buf[1000];
    char *seqfile = argv[1];
    FILE *fp;
    fp = fopen(seqfile, "r");
	while(fgets(buf, sizeof(buf), fp)){
    	char        *x,*y;          /* the two sequences, x_1..M and y_1..N */
    	int         M,N;            /* lengths of sequences x,y             */
    	int         i,j;            /* residue coords in x and y            */
    	int         **S;            /* dynamic programming matrix           */
    	char        *ax,*ay;        /* result: the two aligned sequences    */
    	int         K;              /* length of pairwise alignment         */
    	int         k;              /* index in pairwise alignment, 0..K-1  */
    	int         sc;             /* tmp variable used to find best path  */
    	char        tmp;            /* tmp variable used to swap chars      */
		char *str = strtok(buf, "	 \n");
		char * NAME = str;
		str = strtok(NULL, "	 \n");
		char *SEQX = str;
		str = strtok(NULL, "	 \n");
		char *SEQY = str;

        /*****************************************************************
        * Initial
        *****************************************************************/
        /* Initialize the sequences x,y, using SEQX and SEQY */
        M = strlen(SEQX);
        N = strlen(SEQY);
        x = malloc(sizeof(char) * (M+2)); /* +2: leading blank, and trailing \0 */
        y = malloc(sizeof(char) * (N+2));
        strcpy(x+1, SEQX);            /* x_1..x_M now defined; x_0 undefined */
        strcpy(y+1, SEQY);            /* y_1..y_N now defined; y_0 undefined */

        /* Allocate the dynamic programming matrix, 0..M rows by 0..N columns */
        S = malloc(sizeof(int *) * (M+1));
        for (i = 0; i <= M; i++)
            S[i] = malloc(sizeof(int) * (N+1));


        /*****************************************************************
        * Recursion
        *****************************************************************/
        /* Initialization */
        S[0][0] = 0;
        for (i = 1; i <= M; i++)  S[i][0] = i * INDEL;
        for (j = 1; j <= N; j++)  S[0][j] = j * INDEL;

        /* The dynamic programming, global alignment (Needleman/Wunsch) recursion. */
        for (i = 1; i <= M; i++){
            for (j = 1; j <= N; j++){
                /* case #1: i,j are aligned */
                if (x[i] == y[j]){
                    S[i][j] = S[i-1][j-1] + MATCH;
                /* case #2: mismatch due to bisulfite conversion */
                }else if ((x[i]=='T' && y[j]=='C') || (x[i]=='G' && y[j]=='A')){
                    S[i][j] = S[i-1][j-1] + MATCH;
                /* case #3: other mismatch */
                }else{
                    S[i][j] = S[i-1][j-1] + MISMATCH;
                }
                
                sc = S[i-1][j] + INDEL;               /* case #2: i aligned to -  */
                if (sc > S[i][j]) S[i][j] = sc;

                sc = S[i][j-1] + INDEL;               /* case #3: j aligned to -  */
                if (sc > S[i][j]) S[i][j] = sc;
            }
        }    


        /*****************************************************************
        * Traceback
        *****************************************************************/

        /* Allocate for our aligned strings (which will be 0..K-1) */
        ax = malloc(sizeof(char) * (M+N+1));
        ay = malloc(sizeof(char) * (M+N+1));

        /* Start at M,N; work backwards */
        i = M;
        j = N;
        k = 0;
        while (i > 0 || j > 0){
            if(i*4>=M+N && j*4>=M+N){
                
                /* Case 2: best was x_i aligned to -?  */
                if (i != 0) {  /* watch out for boundaries */
                    if (S[i-1][j] + INDEL == S[i][j]) {
                        ax[k] = x[i]; ay[k] = '-'; k++;     /* record the alignment       */
                        i--;                                /* move to next cell in trace */
                        continue; 
                    }    
                }

                /* Case 3: best was - aligned to y_j? */
                if (j != 0) { /* watch out for boundaries */
                    if (S[i][j-1] + INDEL == S[i][j]) { 
                        ax[k] = '-'; ay[k] = y[j]; k++; /* record the alignment       */
                        j--;                            /* move to next cell in trace */
                        continue;
                    }
                }
                /* Case 1: best was x_i aligned to x_j? */
                if (i > 0 && j > 0) { /* watch out for boundaries */
                    sc = S[i-1][j-1];
                    if (x[i] == y[j] || ((x[i]=='T' && y[j]=='C') || (x[i]=='G' && y[j]=='A'))){
                        sc += MATCH;
                    }else{
                        sc += MISMATCH;
                    }
                    if (sc == S[i][j]) {                    /* residue i aligns to j:     */
                        ax[k] = x[i]; ay[k] = y[j]; k++;    /* record the alignment       */
                        i--; j--;                           /* move to next cell in trace */
                        continue; 
                    } 
                }
            }else{
                
                /* Case 1: best was x_i aligned to x_j? */
                if (i > 0 && j > 0) { /* watch out for boundaries */
                    sc = S[i-1][j-1];
                    if (x[i] == y[j] || ((x[i]=='T' && y[j]=='C') || (x[i]=='G' && y[j]=='A'))){
                        sc += MATCH;
                    }else{
                        sc += MISMATCH;
                    }
                    if (sc == S[i][j]) {                    /* residue i aligns to j:     */
                        ax[k] = x[i]; ay[k] = y[j]; k++;    /* record the alignment       */
                        i--; j--;                           /* move to next cell in trace */
                        continue; 
                    } 
                }
                /* Case 2: best was x_i aligned to -?  */
                if (i != 0) {  /* watch out for boundaries */
                    if (S[i-1][j] + INDEL == S[i][j]) {
                        ax[k] = x[i]; ay[k] = '-'; k++;     /* record the alignment       */
                        i--;                                /* move to next cell in trace */
                        continue; 
                    }    
                }

                /* Case 3: best was - aligned to y_j? */
                if (j != 0) { /* watch out for boundaries */
                    if (S[i][j-1] + INDEL == S[i][j]) { 
                        ax[k] = '-'; ay[k] = y[j]; k++; /* record the alignment       */
                        j--;                            /* move to next cell in trace */
                        continue;
                    }
                }

            }

           
        }
        /* Done with the traceback.
        * Now ax[0..k-1] and ay[0..k-1] are aligned strings, but backwards;
        * and k is our alignment length.
        * Reverse the strings and we're done. 
        */
        K = k;            /* K = remember the alignment length */
        for (k = 0; k < K/2; k++){    /* a sly, efficient in-place reversal by swapping pairs of chars: */
            tmp = ax[K-k-1]; ax[K-k-1] = ax[k]; ax[k] = tmp; 
            tmp = ay[K-k-1]; ay[K-k-1] = ay[k]; ay[k] = tmp; 
        }


        /*****************************************************************
        * Output
        *****************************************************************/
      	printf("%s\t", NAME);
		for(i=0;i<K;i++){
			printf("%c", ax[i]);
		}
		printf("\t");
		for(i=0;i<K;i++){
			printf("%c", ay[i]);
		}
		printf("\n");
        /*printf("%s\n%s\n%s\n", NAME, ax, ay);*/
      
        /* Free the memory we allocated; then exit. */
        free(x);  free(y);
        free(ax); free(ay);
        for (i = 0; i <= M; i++) free(S[i]);
        free(S);
    }
    fclose(fp);
    exit(0);
}

