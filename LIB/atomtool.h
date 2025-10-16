/* Atom counting and manipulation tools for GECKO-A */
#ifndef ATOMTOOL_H
#define ATOMTOOL_H

/* Return the # of carbon atoms in a molecule */
int cnum(const char *chem);

/* Return the # of C-O-C (ether) in a molecule */
int onum(const char *chem);

/* Compute the molecular weight of the species */
void molweight(const char *chem, double *weight);

/* Find the number of each possible atom in the species CHEM */
void getatoms(const char *chem, int *ica, int *ih, int *init, int *io, 
              int *ir, int *isul, int *ifl, int *ibr, int *icl);

/* Calculate mean oxidation number of carbon */
void oxnum(const char *chem, double *numoxC);

#endif /* ATOMTOOL_H */
