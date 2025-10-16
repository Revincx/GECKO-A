/* Toolbox utility functions for GECKO-A */
#ifndef TOOLBOX_H
#define TOOLBOX_H

#include <stdbool.h>

/* Add a reference to the reference list */
void addref(const char *progname, const char *reference, int *nref, 
            char **reflist, int maxrefs, const char *chem);

/* Increment number of reactions */
void addrx(const char *progname, const char *chem, int *nr, 
           int *flag, int maxflag);

/* Increment number of products in a reaction */
void add1tonp(const char *progname, const char *chem, int *np);

/* Generate error message output */
void stoperr(const char *prog, const char *mesg, const char *chem);

/* Set bond value (symmetric) */
void setbond(int **bondtb, int x, int y, int bondvalue);

/* Remove empty nodes in group and bond */
void erase_blank(int **bond, char **group, int size1, int size2);

/* Count the number of occurrences of a substring in a string */
int countstring(const char *line, const char *str);

/* Return the rate constant at T for the Arrhenius parameters */
double kval(const double *arrh, double T);

#endif /* TOOLBOX_H */
