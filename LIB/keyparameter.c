/* Set the list of key parameters for gecko */
#include "keyparameter.h"

/* directory for the input/output files */
char dirgecko[512] = "";
char dirout[512] = "";

/* digit (to rate groups and ring joining char) */
const char digit[4] = {'1', '2', '3', '4'};

/* prime numbers */
const int prim[35] = {
    2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
   31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
   73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149
};
