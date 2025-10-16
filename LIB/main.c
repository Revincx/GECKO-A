/* GECKO-A Main Program - C Version
 * 
 * PURPOSE: Generate the oxidation mechanism for organic compounds under
 * tropospheric conditions.
 */

#include <stdio.h>
#include <stdlib.h>
#include "keyparameter.h"
#include "keyflag.h"
#include "minidict.h"
#include "tempoci.h"
#include "references.h"
#include "tempflag.h"

int main(int argc, char *argv[]) {
    printf("===========================================\n");
    printf("GECKO-A - C Version\n");
    printf("===========================================\n");
    printf("Generate the oxidation mechanism for organic\n");
    printf("compounds under tropospheric conditions.\n");
    printf("===========================================\n\n");
    
    /* Initialize default parameters */
    printf("Initializing default parameters...\n");
    define_defaults();
    
    printf("Directory for GECKO: %s\n", dirgecko);
    printf("Output directory: %s\n", dirout);
    printf("Maximum generations: %d\n", maxgen);
    printf("Temperature: %.1f K\n", TK);
    printf("Critical vapor pressure: %.1f\n", critvp);
    printf("\n");
    
    /* Test minidict functionality */
    printf("Testing minidict module...\n");
    clean_minid();
    add_fo("CH4", "CH4");
    add_fo("C2H6", "CC");
    
    char formula[MXLFO+1];
    get_fo("CH4", formula);
    printf("Formula for CH4: %s\n", formula);
    get_fo("C2H6", formula);
    printf("Formula for C2H6: %s\n", formula);
    printf("Number of entries in minidict: %d\n", nmini);
    printf("\n");
    
    /* Display some constants */
    printf("Key parameters:\n");
    printf("  MXLCO (max species name length): %d\n", MXLCO);
    printf("  MXLFO (max formula length): %d\n", MXLFO);
    printf("  MXNODE (max nodes): %d\n", MXNODE);
    printf("  MXPS (max primary species): %d\n", MXPS);
    printf("\n");
    
    printf("Prime numbers array (first 10):\n  ");
    for (int i = 0; i < 10; i++) {
        printf("%d ", prim[i]);
    }
    printf("\n\n");
    
    printf("GECKO-A C version initialized successfully!\n");
    printf("Note: This is a minimal demonstration.\n");
    printf("Full mechanism generation not yet implemented.\n");
    
    return 0;
}
