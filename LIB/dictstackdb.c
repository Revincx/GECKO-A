/* Dictionary and stack database module */
#include "dictstackdb.h"

/* dictionaries */
int nrec = 0;
int ninorg = 0;
int nwpspe = 0;
char dict[MXSPE][MXLDI+1];
char namlst[MXSPE][MXLCO+1];
char inorglst[MXC1][MXLDI+1];
double dbrch[MXSPE];

/* to be introduced in a "type dict" */
double dctmw[MXSPE];
int dctatom[MXSPE][9];
double dcthenry[MXSPE];
double dctnan[MXSPE][2];
double dctsim[MXSPE][2];
double dctmyr[MXSPE][2];

/* stack */
int nhldvoc = 0;
char holdvoc[MXSVOC][LENSS+1];
int nhldrad = 0;
char holdrad[MXSRAD][LENSS+1];
int stabl = 0;
int level = 0;
bool lotopstack = false;

/* isomer stack */
int diccri[MXSPE][MXCRI];

/* tetrahydrofuran (CHA) stack */
int ncha = 0;
char chatab[MXCHA][MXLCO+1];
