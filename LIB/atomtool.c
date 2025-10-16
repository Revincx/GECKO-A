/* Atom counting and manipulation tools for GECKO-A */
#include "atomtool.h"
#include "keyparameter.h"
#include "rjtool.h"
#include "toolbox.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/* Return the # of carbon atoms in a molecule */
int cnum(const char *chem) {
    int ic = 0;
    int nc = strcspn(chem, " ");
    if (nc == 0) nc = strlen(chem);
    
    for (int i = 0; i < nc; i++) {
        if (chem[i] == 'C' && !(i + 1 < nc && chem[i+1] == 'l')) {
            ic++;
        }
        if (chem[i] == 'c') {
            ic++;
        }
    }
    
    return ic;
}

/* Return the # of C-O-C (ether) in a molecule */
int onum(const char *chem) {
    int io = 0;
    int nc = strcspn(chem, " ");
    if (nc == 0) nc = strlen(chem);
    if (nc > 0) nc--;
    
    for (int i = 0; i < nc; i++) {
        if (i + 1 < (int)strlen(chem) && strncmp(chem + i, "-O", 2) == 0) {
            io++;
        }
    }
    
    return io;
}

/* Compute the molecular weight of the species */
void molweight(const char *chem, double *weight) {
    int ica, iha, ina, ioa, ira, isa, ifl, ibr, icl;
    
    getatoms(chem, &ica, &iha, &ina, &ioa, &ira, &isa, &ifl, &ibr, &icl);
    
    *weight = ica * 12.01 + iha * 1.01 + ina * 14.01 + ioa * 16.0 + 
              isa * 32.07 + ifl * 19.0 + ibr * 79.9 + icl * 35.45;
}

/* Find the number of each possible atom in the species CHEM */
void getatoms(const char *chem, int *ica, int *ih, int *init, int *io, 
              int *ir, int *isul, int *ifl, int *ibr, int *icl) {
    int nring, i, n, nc;
    char *tchem;
    int **rjs;
    bool loring;
    
    /* Initialize */
    nc = strcspn(chem, " ");
    if (nc == 0) nc = strlen(chem);
    
    tchem = (char *)malloc(strlen(chem) + 1);
    strcpy(tchem, chem);
    
    *ifl = 0; *ibr = 0; *icl = 0;
    *ica = 0; *ih = 0; *init = 0; *io = 0;
    *isul = 0; *ir = 0;
    n = 1;
    
    /* Check for rings */
    loring = false;
    if (strstr(tchem, "C1") != NULL) loring = true;
    if (strstr(tchem, "Cd1") != NULL) loring = true;
    if (strstr(tchem, "c1") != NULL) loring = true;
    if (strstr(tchem, "-O1") != NULL) loring = true;
    
    if (loring) {
        nring = 2;  /* set max # of rings */
        rjs = (int **)malloc(nring * sizeof(int *));
        for (i = 0; i < nring; i++) {
            rjs[i] = (int *)malloc(2 * sizeof(int));
        }
        rjsrm(nring, tchem, rjs, nring);
        
        for (i = 0; i < nring; i++) {
            free(rjs[i]);
        }
        free(rjs);
    }
    
    /* Start counting from the end of the string */
    i = nc;
    while (i > 0) {
        i--;
        
        if (tchem[i] == '2') n = 2;
        if (tchem[i] == '3') n = 3;
        if (tchem[i] == '4') n = 4;
        
        if (i > 0 && strncmp(tchem + i - 1, "Cl", 2) == 0) {
            *icl += n;
            n = 1;
            i--;
            continue;
        }
        
        if (tchem[i] == 'F') {
            *ifl += n;
            n = 1;
        } else if (i + 1 < nc && strncmp(tchem + i, "Br", 2) == 0) {
            *ibr += n;
            n = 1;
        } else if (tchem[i] == 'H') {
            *ih += n;
            n = 1;
        } else if (tchem[i] == 'N') {
            *init += n;
            n = 1;
        } else if (tchem[i] == 'O') {
            *io += n;
            n = 1;
        } else if (tchem[i] == 'S') {
            *isul += 1;
            n = 1;
        } else if (tchem[i] == '.') {
            *ir += n;
            n = 1;
        } else if ((tchem[i] == 'C' || tchem[i] == 'c') && 
                   !(i + 1 < nc && tchem[i+1] == 'l')) {
            *ica += 1;
            n = 1;
        }
    }
    
    free(tchem);
}

/* Calculate mean oxidation number of carbon - stub for now */
void oxnum(const char *chem, double *numoxC) {
    /* This function requires stdgrbond which hasn't been converted yet */
    /* For now, return 0 as a placeholder */
    (void)chem;  /* Unused */
    *numoxC = 0.0;
    
    /* TODO: Implement when stdgrbond is converted */
}
