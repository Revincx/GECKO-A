/* Temporary OCI data module */
#ifndef TEMPOCI_H
#define TEMPOCI_H

extern double primoOH, secoOH, primosci, ycarbo;

#define MXSCI 15
extern int nnsci;
extern double yysci[MXSCI];
extern char chemsci[MXSCI][120];

#define MXCARB 6
extern char focarb[MXCARB][120];
extern int nc_carb[MXCARB];
extern double ycarb[MXCARB];

extern int c1ciflag;

#endif /* TEMPOCI_H */
