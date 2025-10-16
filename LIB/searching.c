/* Searching utilities for GECKO-A */
#include "searching.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Binary tree search for string */
int srh5(const char *aseek, char **alist, int nlist) {
    int jhi, jlo, jold, j;
    
    /* initialize */
    jold = 0;
    jlo = 0;  /* C arrays are 0-indexed, so start from 0 */
    jhi = nlist;
    
    /* search loop */
    while (1) {
        j = (jhi + jlo) / 2;
        if (j == jold) break;
        
        jold = j;
        
        int cmp = strcmp(aseek, alist[j]);
        if (cmp > 0) {
            jlo = j;
            continue;
        }
        if (cmp == 0) {
            return j + 1;  /* Return 1-indexed position for compatibility */
        }
        jhi = j;
    }
    
    /* string not found */
    return -(j + 1);  /* Return negative 1-indexed position */
}

/* Binary tree search for formula in dictionary */
int srch(int nrec, const char *chem, char **dict) {
    int jhi, jlo, jold, j;
    char dict_substr[121];  /* For characters 10:129 (0-indexed: 9:129) */
    
    /* initialize */
    jold = 0;
    jlo = 0;  /* C arrays are 0-indexed */
    jhi = nrec;
    
    /* search loop */
    while (1) {
        j = (jhi + jlo) / 2;
        if (j == jold) break;
        
        jold = j;
        
        /* Extract substring from position 10:129 (Fortran) = 9:129 (C, 0-indexed) */
        /* dict(j)(10:129) in Fortran means characters 10 through 129 */
        if (strlen(dict[j]) >= 129) {
            strncpy(dict_substr, dict[j] + 9, 120);
            dict_substr[120] = '\0';
        } else if (strlen(dict[j]) > 9) {
            strcpy(dict_substr, dict[j] + 9);
        } else {
            dict_substr[0] = '\0';
        }
        
        int cmp = strcmp(chem, dict_substr);
        if (cmp > 0) {
            jlo = j;
            continue;
        }
        if (cmp == 0) {
            return j + 1;  /* Return 1-indexed position for compatibility */
        }
        jhi = j;
    }
    
    /* string not found */
    return -(j + 1);  /* Return negative 1-indexed position */
}

/* Return the rank (position) of an integer in a sorted list */
int search_ipos(int iseek, const int *ilist, int list_size) {
    int i;
    
    if (iseek < 1) {
        fprintf(stderr, "--error-- in search_ipos, requested search using <1 int.\n");
        exit(EXIT_FAILURE);
    }
    
    for (i = 0; i < list_size; i++) {
        if (iseek > ilist[i]) {
            return i + 1;  /* Return 1-indexed position for compatibility */
        }
    }
    
    fprintf(stderr, "--error-- in search_ipos, no slot found\n");
    exit(EXIT_FAILURE);
}
