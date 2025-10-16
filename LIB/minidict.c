/* Mini dictionary module for GECKO-A */
#include "minidict.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/* Global data definitions */
int nmini = 0;
char nam_minid[MXMINID][MXLCO+1];
char fo_minid[MXMINID][MXLFO+1];

/* Get formula for a given name */
void get_fo(const char *nam, char *formula) {
    int i;
    
    formula[0] = '\0';  /* initialize to empty string */
    
    for (i = 0; i < nmini; i++) {
        if (strcmp(nam_minid[i], nam) == 0) {
            strcpy(formula, fo_minid[i]);
            return;
        }
    }
    
    /* no formula found if that point is reached */
    strcpy(formula, nam);
}

/* Add formula and name to the mini dictionary */
void add_fo(const char *nam, const char *formula) {
    int i;
    bool toadd = true;
    
    for (i = 0; i < nmini; i++) {
        if (strcmp(formula, fo_minid[i]) == 0) {
            toadd = false;
            break;
        }
    }
    
    if (toadd) {
        if (nmini >= MXMINID) {
            fprintf(stderr, "maximum formula in minid reached - see module minidict\n");
            exit(EXIT_FAILURE);
        }
        
        strcpy(fo_minid[nmini], formula);
        strcpy(nam_minid[nmini], nam);
        nmini++;
    }
}

/* Clean/reset the mini dictionary */
void clean_minid(void) {
    int i;
    nmini = 0;
    
    for (i = 0; i < MXMINID; i++) {
        nam_minid[i][0] = '\0';
        fo_minid[i][0] = '\0';
    }
}
