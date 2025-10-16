/* Temporary flags module */
#ifndef TEMPFLAG_H
#define TEMPFLAG_H

/* ignore reaction products (e.g. after a given # of generation) */
extern int iflost;   /* flag raised if products are not written */
extern int xxc;      /* # of C atoms in the reactant */
extern int xxh;      /* # of H atoms in the reactant */
extern int xxn;      /* # of N atoms in the reactant */
extern int xxo;      /* # of O atoms in the reactant */
extern int xxr;      /* # of radical dots in the reactant */
extern int xxs;      /* # of S atoms in the reactant */
extern int xxfl;     /* # of F atoms in the reactant */
extern int xxbr;     /* # of Br atoms in the reactant */
extern int xxcl;     /* # of Cl atoms in the reactant */

#endif /* TEMPFLAG_H */
