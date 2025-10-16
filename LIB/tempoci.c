/* Temporary OCI data module */
#include "tempoci.h"

double primoOH, secoOH, primosci, ycarbo;

int nnsci;
double yysci[MXSCI];
char chemsci[MXSCI][120];

char focarb[MXCARB][120];
int nc_carb[MXCARB];
double ycarb[MXCARB];

int c1ciflag;
