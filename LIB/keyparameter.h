/* Set the list of key parameters for gecko */
#ifndef KEYPARAMETER_H
#define KEYPARAMETER_H

/* maximum length and size parameters */
#define MXLCO 6       /* maximum length of species names (code) */
#define MXLFO 120     /* maximum length of a formula */
#define MXLFL 15      /* maximum length of functionality list (string) */
#define MXLDI 146     /* string length in the dictionary */
#define MXNODE 35     /* maximum number of nodes allowed */
#define MXLGR 21      /* maximum length of a string in a group */
#define MXRING 4      /* maximum rings allowed */
#define MXCP 99       /* maximum # of copies of formula allowed & max # of "tracks" */
#define MXHYD 10      /* maximum # of distinct position to add H2O (i.e. # of carbonyls) */
#define MXHISO 5000   /* maximum # of hydrate isomer a molecule can have */
#define MXPD 38       /* maximum # of products per reaction */
#define MXNR 35       /* maximum # of reactions per species */
#define MXCOPD 10     /* maximum # of coproducts per generated species */
#define MXNP 4        /* maximum # of product per reaction in the output mechanism */
#define MXPS 600      /* maximum # of "primary" species that can be given as input */
#define MXTRK 5       /* maximum # of distinct Cd track in a molecule */
#define MXLCD 4       /* maximum length of a Cd track (current max is 4 for C=C-C=C) */
#define MXLEST 6      /* maximum length of an ester track (current max is 6: CO-O-CO-O-CO-O) */

/* file unit numbers */
#define DCTU 7        /* dictionary file unit */
#define PRMU 14       /* list of primary species used (for findname) */
#define LOGU 15       /* log file unit */
#define REFU 16       /* mechanism with reference file unit */
#define MECU 17       /* mechanism file unit */
#define SCRU 18       /* file unit to redirect info from screen to file */
#define KOHU 19       /* kivoci file unit (record k_OH for the various species) */
#define KNO3U 59      /* kjvocj file unit (record k_NO3 for the various species) */
#define WARU 20       /* warning file unit */
#define GASU 21       /* gas phase species file unit */
#define PRTU 22       /* particle phase species file unit */
#define WALU 23       /* wall phase species file unit */
#define OHU 26        /* rate constant OH file unit (e.g. SAR assessment) */
#define O3U 27        /* rate constant O3 file unit (e.g. SAR assessment) */
#define NO3U 28       /* rate constant NO3 file unit (e.g. SAR assessment) */
#define DHFU 29       /* dummy info file for dhf reactions */
#define SARU 30       /* SAR info about group ... */
#define TFU1 31       /* temporary file unit 1 */
#define TFU2 32       /* temporary file unit 2 */
#define TFU3 33       /* temporary file unit 3 */
#define TFU4 34       /* temporary file unit 4 */
#define MCOU 37       /* mechanism file unit with comment! */
#define PFU1 41       /* peroxy file unit 1 */
#define PFU2 42       /* peroxy file unit 2 */
#define PFU3 43       /* peroxy file unit 3 */
#define PFU4 44       /* peroxy file unit 4 */
#define PFU5 45       /* peroxy file unit 5 */
#define PFU6 46       /* peroxy file unit 6 */
#define PFU7 47       /* peroxy file unit 7 */
#define PFU8 48       /* peroxy file unit 8 */
#define PFU9 49       /* peroxy file unit 9 */
#define NMLU 70       /* namelist file unit */
#define RSZU 789      /* record stack size unit (information file) */

/* directory for the input/output files */
extern char dirgecko[512];
extern char dirout[512];

/* digit (to rate groups and ring joining char) */
extern const char digit[4];

/* priorities to rate/sort functionalities in groups */
#define PRI "CHOCHC CH(CH2 CH3O.)OO. OOH ONOF  Br2Br)Cl2Cl)NO2NO)OH)"

/* prime numbers */
extern const int prim[35];

#endif /* KEYPARAMETER_H */
