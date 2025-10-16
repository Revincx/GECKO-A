/* Toolbox utility functions for GECKO-A */
#include "toolbox.h"
#include "keyparameter.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Add a reference to the reference list */
void addref(const char *progname, const char *reference, int *nref, 
            char **reflist, int maxrefs, const char *chem) {
    char species[MXLFO+1];
    int i;
    
    /* Check if chem is provided */
    if (chem != NULL) {
        strncpy(species, chem, MXLFO);
        species[MXLFO] = '\0';
    } else {
        species[0] = '\0';
    }
    
    /* Check if reference already in list */
    for (i = 0; i < *nref; i++) {
        if (strcmp(reflist[i], reference) == 0) {
            return;  /* ref already in list */
        }
    }
    
    (*nref)++;
    if (*nref > maxrefs) {
        fprintf(stderr, "--error-- too many references added for the list\n");
        fprintf(stderr, "          in subroutine: %s\n", progname);
        fprintf(stderr, "size of reflist is: %d\n", maxrefs);
        if (species[0] != '\0') {
            fprintf(stderr, "          for the species: %s\n", species);
        }
        fprintf(stderr, "          reference was: %s\n", reference);
        exit(EXIT_FAILURE);
    }
    strcpy(reflist[*nref - 1], reference);
}

/* Increment number of reactions */
void addrx(const char *progname, const char *chem, int *nr, 
           int *flag, int maxflag) {
    char mesg[200];
    
    (*nr)++;
    if (*nr > maxflag) {
        snprintf(mesg, sizeof(mesg), 
                "-error- too many reactions created for species: %s in addrx", chem);
        stoperr(progname, mesg, chem);
    }
    flag[*nr - 1] = 1;
}

/* Increment number of products in a reaction */
void add1tonp(const char *progname, const char *chem, int *np) {
    (*np)++;
    if (*np > MXPD) {
        fprintf(stderr, "--error-- too many products created, np > mxpd\n");
        fprintf(stderr, "          for species: %s\n", chem);
        fprintf(stderr, "          in subroutine: %s\n", progname);
        exit(EXIT_FAILURE);
    }
}

/* Generate error message output */
void stoperr(const char *prog, const char *mesg, const char *chem) {
    const bool kill = true;  /* Set to false to continue on error */
    
    fprintf(stderr, "\n");
    fprintf(stderr, "--error-- in: %s\n", prog);
    fprintf(stderr, "%s\n", mesg);
    fprintf(stderr, "for species: %s\n", chem);
    
    if (kill) {
        exit(EXIT_FAILURE);
    }
}

/* Set bond value (symmetric) */
void setbond(int **bondtb, int x, int y, int bondvalue) {
    bondtb[x][y] = bondvalue;
    bondtb[y][x] = bondvalue;
}

/* Remove empty nodes in group and bond */
void erase_blank(int **bond, char **group, int size1, int size2) {
    int i, j, k, last, n, nca, icheck;
    int **tbond;
    char **tgroup;
    
    /* Allocate temporary arrays */
    tbond = (int **)malloc(size1 * sizeof(int *));
    for (i = 0; i < size1; i++) {
        tbond[i] = (int *)malloc(size2 * sizeof(int));
        memcpy(tbond[i], bond[i], size2 * sizeof(int));
    }
    
    tgroup = (char **)malloc(size1 * sizeof(char *));
    for (i = 0; i < size1; i++) {
        tgroup[i] = (char *)malloc(strlen(group[i]) + 1);
        strcpy(tgroup[i], group[i]);
    }
    
    /* Count non-empty groups */
    nca = 0;
    last = 0;
    for (n = 0; n < size1; n++) {
        if (tgroup[n][0] != ' ' && tgroup[n][0] != '\0') {
            nca++;
            last = n;
        }
    }
    
    if (nca < 2 || last == nca - 1) {
        /* Clean up and return */
        for (i = 0; i < size1; i++) {
            free(tbond[i]);
            free(tgroup[i]);
        }
        free(tbond);
        free(tgroup);
        return;
    }
    
    /* Erase blank lines in bond and group */
    icheck = 0;
    i = 0;
    
    while (i < nca) {
        if (tgroup[i][0] == ' ' || tgroup[i][0] == '\0') {
            /* Shift groups */
            for (j = i; j < last; j++) {
                strcpy(tgroup[j], tgroup[j+1]);
            }
            tgroup[last][0] = '\0';
            
            /* Shift bond rows */
            for (j = 0; j <= last; j++) {
                for (k = i; k < last; k++) {
                    tbond[k][j] = tbond[k+1][j];
                }
                tbond[last][j] = 0;
            }
            
            /* Shift bond columns */
            for (j = 0; j <= last; j++) {
                for (k = i; k < last; k++) {
                    tbond[j][k] = tbond[j][k+1];
                }
                tbond[j][last] = 0;
            }
        } else {
            i++;
        }
        
        icheck++;
        if (icheck > size1) {
            fprintf(stderr, "--error-- in erase_blank\n");
            fprintf(stderr, "infinite loop when erasing blank lines\n");
            exit(EXIT_FAILURE);
        }
    }
    
    /* Copy back results */
    for (i = 0; i < size1; i++) {
        memcpy(bond[i], tbond[i], size2 * sizeof(int));
        strcpy(group[i], tgroup[i]);
    }
    
    /* Clean up */
    for (i = 0; i < size1; i++) {
        free(tbond[i]);
        free(tgroup[i]);
    }
    free(tbond);
    free(tgroup);
}

/* Count the number of occurrences of a substring in a string */
int countstring(const char *line, const char *str) {
    int n = 0;
    const char *p;
    size_t len;
    
    len = strlen(str);
    if (len == 0) return 0;
    
    p = line;
    while ((p = strstr(p, str)) != NULL) {
        n++;
        p += len;
    }
    
    return n;
}

/* Return the rate constant at T for the Arrhenius parameters */
double kval(const double *arrh, double T) {
    /* arrh[0] = A, arrh[1] = n, arrh[2] = Ea/R */
    return arrh[0] * pow(T, arrh[1]) * exp(-arrh[2] / T);
}
