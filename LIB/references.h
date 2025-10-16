/* References module */
#ifndef REFERENCES_H
#define REFERENCES_H

#define MXREAC_INFO 1000  /* max # of comment/references */
#define MXLCOD 10         /* max length of comment's code in database */
#define MXLREAC_INFO 400  /* max length of a comment/reference */

extern int nreac_info;                          /* # of available "reac_info" (references) */
extern char code[MXREAC_INFO][MXLCOD+1];        /* Code for each "reac_info" */
extern char reac_info[MXREAC_INFO][MXLREAC_INFO+1]; /* Full text corresponding to the reac_info */

#endif /* REFERENCES_H */
