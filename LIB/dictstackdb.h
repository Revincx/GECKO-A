/* Dictionary and stack database module */
#ifndef DICTSTACKDB_H
#define DICTSTACKDB_H

#include "keyparameter.h"
#include <stdbool.h>

/* dictionaries */
#define MXSPE 7000000     /* max # of species allowed in the mechanism */
#define MXC1 3000         /* max # of inorganic and C1 species allowed in the mechanism */

extern int nrec;                        /* # of species recorded in dict */
extern int ninorg;                      /* # of inorganic species */
extern int nwpspe;                      /* # of species recorded in the condensed (wall,part.) phase */
extern char dict[MXSPE][MXLDI+1];       /* dictionary line (code+formula+fg) */
extern char namlst[MXSPE][MXLCO+1];     /* name (lco=6 characters) of the species used */
extern char inorglst[MXC1][MXLDI+1];    /* list of inorganic species (code+formula+fg) */
extern double dbrch[MXSPE];             /* yield attach to a formula in dict */

/* to be introduced in a "type dict" */
extern double dctmw[MXSPE];             /* molecular weight of species in dict */
extern int dctatom[MXSPE][9];           /* atoms & radical for species in dict */
extern double dcthenry[MXSPE];          /* Henry's law coeff. for species in dict */
extern double dctnan[MXSPE][2];         /* Pvap(:,1) & heat(:,2) (Nannoonal) for species in dict */
extern double dctsim[MXSPE][2];         /* Pvap(:,1) & heat(:,2) (Simpol) for species in dict */
extern double dctmyr[MXSPE][2];         /* Pvap(:,1) & heat(:,2) (Myrdal&Yalkowski) for species in dict */

/* stack */
#define MXSVOC 1500000    /* max # of species in the voc stack */
#define MXSRAD 200        /* max # of species in the radical stack */
#define LENSS 132         /* length of a stack string (code+formula+i3+i3) */

extern int nhldvoc;                     /* # of (stable) VOC in the stack */
extern char holdvoc[MXSVOC][LENSS+1];   /* VOCs in the stack (name[a6]+formula[a120]+stabl[i3]+level[i3]) */
extern int nhldrad;                     /* # of radical in the stack */
extern char holdrad[MXSRAD][LENSS+1];   /* radicals in the stack(name[a6]+formula[a120]+stabl[i3]+level[i3]) */
extern int stabl;                       /* # of generations needed to produce the current species */
extern int level;                       /* # of intermediates (rad + stable) needed to produce the current species */
extern bool lotopstack;                 /* if true, add species on top of VOC stack (default is false) */

/* isomer stack */
#define MXISO 2000        /* max # of isomers for a given molecule */
#define MXCRI 21          /* max # of criteria used to discriminate isomers */
extern int diccri[MXSPE][MXCRI];

/* tetrahydrofuran (CHA) stack (to manage phase partitioning) */
#define MXCHA 10000       /* max # of tetrahydrofuran (CHA) allowed */
extern int ncha;                        /* # of species recorded in dict */
extern char chatab[MXCHA][MXLCO+1];     /* idnam of the CHA */

#endif /* DICTSTACKDB_H */
