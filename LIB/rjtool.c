/* Ring-joining tool module for GECKO-A */
#include "rjtool.h"
#include "keyparameter.h"
#include <string.h>
#include <stdio.h>

/* Add ring-joining numerical characters to group strings */
void rjgadd(int nring, char **group, int **rjg, int group_size, int rjg_size) {
    int n, i, ii, j;
    char temp[256];
    
    /* loop in reverse order => chars in numerical order if >1 exist at any node */
    for (n = nring - 1; n >= 0; n--) {
        for (ii = 0; ii < 2; ii++) {
            i = rjg[n][ii] - 1;  /* Convert to 0-indexed */
            if (i >= 0) {
                if ((strlen(group[i]) >= 2 && strncmp(group[i], "-O", 2) == 0) ||
                    (strlen(group[i]) >= 2 && group[i][1] == 'd')) {
                    j = 2;  /* Position 3 in Fortran = index 2 in C */
                } else {
                    j = 1;  /* Position 2 in Fortran = index 1 in C */
                }
                
                /* Insert digit at position j */
                strcpy(temp, group[i]);
                temp[j] = gecko_digit[n];
                strcpy(temp + j + 1, group[i] + j);
                strcpy(group[i], temp);
            }
        }
    }
}

/* Strip ring-joining numerical characters from group strings */
void rjgrm(int nring, char **group, int **rjg, int group_size, int rjg_size) {
    int n, i, j;
    
    /* Initialize rjg to 0 */
    for (n = 0; n < rjg_size; n++) {
        rjg[n][0] = 0;
        rjg[n][1] = 0;
    }
    
    for (n = 0; n < nring; n++) {
        for (i = 0; i < group_size; i++) {
            if ((strlen(group[i]) >= 2 && strncmp(group[i], "-O", 2) == 0) ||
                (strlen(group[i]) >= 2 && group[i][1] == 'd')) {
                j = 2;  /* Position 3 in Fortran = index 2 in C */
            } else {
                j = 1;  /* Position 2 in Fortran = index 1 in C */
            }
            
            if (group[i][j] == gecko_digit[n]) {
                /* Remove the digit by shifting the string */
                memmove(group[i] + j, group[i] + j + 1, strlen(group[i] + j + 1) + 1);
                
                if (rjg[n][0] == 0) {
                    rjg[n][0] = i + 1;  /* Store as 1-indexed */
                } else {
                    rjg[n][1] = i + 1;  /* Store as 1-indexed */
                }
            }
        }
    }
}

/* Add ring-joining numerical characters to chem string */
void rjsadd(int nring, char *chem, int **rjs, int rjs_size) {
    int n, i, j;
    char temp[256];
    
    /* loop in reverse order => chars in numerical order if >1 exist at any node */
    for (n = nring - 1; n >= 0; n--) {
        for (j = 1; j >= 0; j--) {
            i = rjs[n][j] - 1;  /* Convert to 0-indexed */
            if (i >= 0 && i < (int)strlen(chem)) {
                /* Insert digit at position i */
                strcpy(temp, chem);
                temp[i] = gecko_digit[n];
                strcpy(temp + i + 1, chem + i);
                strcpy(chem, temp);
            }
        }
    }
}

/* Strip ring-joining numerical characters from chem string */
void rjsrm(int nring, char *chem, int **rjs, int rjs_size) {
    int n, i, j, ptr, lenchem;
    
    /* Initialize rjs to 0 */
    for (n = 0; n < rjs_size; n++) {
        rjs[n][0] = 0;
        rjs[n][1] = 0;
    }
    
    /* Find length of chem (up to first space or end) */
    lenchem = strcspn(chem, " ");
    if (lenchem == 0) lenchem = strlen(chem);
    
    for (n = 0; n < nring; n++) {
        for (j = 0; j < 2; j++) {
            ptr = 0;
            for (i = 0; i < lenchem; i++) {
                if (chem[i] == 'C' || chem[i] == 'c' || chem[i] == '-') {
                    ptr = i + 1;
                    if (ptr < lenchem && chem[ptr] == 'd') ptr++;
                    if (ptr >= 2 && strncmp(chem + ptr - 2, "-O", 2) == 0) ptr++;
                    
                    if (ptr < lenchem && chem[ptr] == gecko_digit[n]) {
                        rjs[n][j] = ptr + 1;  /* Store as 1-indexed */
                        /* Remove the digit */
                        memmove(chem + ptr, chem + ptr + 1, strlen(chem + ptr + 1) + 1);
                        lenchem--;
                        break;
                    }
                }
            }
        }
    }
}
