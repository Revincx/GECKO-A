/* Prime tool module for GECKO-A */
#include "primetool.h"
#include "keyparameter.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/* Rank numerical priorities */
void sortrank(int nca, int **prank, int *rank, int prank_rows, int prank_cols) {
    int **pr;
    int store[3];
    int i, ii, j, ctr;
    
    /* Just to be sure */
    if (prank_cols != 3) {
        fprintf(stderr, "in sortrank, size array issue\n");
        exit(EXIT_FAILURE);
    }
    
    /* Allocate and copy prank to pr */
    pr = (int **)malloc(prank_rows * sizeof(int *));
    for (i = 0; i < prank_rows; i++) {
        pr[i] = (int *)malloc(prank_cols * sizeof(int));
        memcpy(pr[i], prank[i], prank_cols * sizeof(int));
    }
    
    /* First pass - bubble sort by column 2 */
    bool swapped = true;
    while (swapped) {
        swapped = false;
        for (i = 0; i < nca - 1; i++) {
            ii = i + 1;
            if (pr[i][1] > pr[ii][1]) {  /* Column 2 is index 1 (0-indexed) */
                for (j = 0; j < 3; j++) {
                    store[j] = pr[ii][j];
                    pr[ii][j] = pr[i][j];
                    pr[i][j] = store[j];
                }
                swapped = true;
            }
        }
    }
    
    /* Second pass - bubble sort by column 3 within same column 2 */
    swapped = true;
    while (swapped) {
        swapped = false;
        for (i = 0; i < nca - 1; i++) {
            ii = i + 1;
            if (pr[i][1] == pr[ii][1] && pr[i][2] > pr[ii][2]) {
                for (j = 0; j < 3; j++) {
                    store[j] = pr[ii][j];
                    pr[ii][j] = pr[i][j];
                    pr[i][j] = store[j];
                }
                swapped = true;
            }
        }
    }
    
    /* Find new ranks */
    ctr = 1;
    rank[pr[0][0] - 1] = ctr;  /* Convert to 0-indexed */
    
    for (ii = 1; ii < nca; ii++) {
        i = ii - 1;
        if (pr[ii][1] == pr[i][1] && pr[ii][2] == pr[i][2]) {
            /* ctr stays same */
        } else {
            ctr++;
        }
        rank[pr[ii][0] - 1] = ctr;  /* Convert to 0-indexed */
    }
    
    /* Free allocated memory */
    for (i = 0; i < prank_rows; i++) {
        free(pr[i]);
    }
    free(pr);
}

/* Update rank of nodes */
void primes(int nca, int **bond, int *rank, int bond_size) {
    int *cprim;
    int *pp;
    int **prank;
    int i, j, maxr;
    bool changed;
    
    cprim = (int *)malloc(bond_size * sizeof(int));
    pp = (int *)malloc(bond_size * sizeof(int));
    prank = (int **)malloc(bond_size * sizeof(int *));
    for (i = 0; i < bond_size; i++) {
        prank[i] = (int *)malloc(3 * sizeof(int));
    }
    
    /* Re-entry point if iteration required */
    do {
        /* Assign primes */
        for (i = 0; i < nca; i++) {
            cprim[i] = prim[rank[i] - 1];  /* rank is 1-indexed */
            pp[i] = 1;
        }
        
        /* Find product of connected primes for nodes */
        for (i = 0; i < nca; i++) {
            for (j = 0; j < nca; j++) {
                if (i != j && bond[i][j] > 0) {
                    pp[i] *= cprim[j];
                }
            }
        }
        
        /* Construct rank vectors for nodes */
        maxr = 0;
        for (i = 0; i < nca; i++) {
            prank[i][0] = i + 1;  /* Store as 1-indexed */
            prank[i][1] = rank[i];
            prank[i][2] = pp[i];
            if (rank[i] > maxr) maxr = rank[i];
        }
        
        /* Save old ranks to check for changes */
        int *old_rank = (int *)malloc(nca * sizeof(int));
        memcpy(old_rank, rank, nca * sizeof(int));
        
        sortrank(nca, prank, rank, bond_size, 3);
        
        /* Check if rank changed */
        changed = false;
        for (i = 0; i < nca; i++) {
            if (rank[i] != old_rank[i]) {
                changed = true;
                break;
            }
        }
        
        free(old_rank);
        
    } while (changed);
    
    /* Free allocated memory */
    free(cprim);
    free(pp);
    for (i = 0; i < bond_size; i++) {
        free(prank[i]);
    }
    free(prank);
}
