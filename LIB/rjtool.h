/* Ring-joining tool module for GECKO-A */
#ifndef RJTOOL_H
#define RJTOOL_H

/* Add ring-joining numerical characters to group strings at nodes given by rjg */
void rjgadd(int nring, char **group, int **rjg, int group_size, int rjg_size);

/* Strip ring-joining numerical characters from group strings */
void rjgrm(int nring, char **group, int **rjg, int group_size, int rjg_size);

/* Add ring-joining numerical characters to group strings at nodes given by rjs */
void rjsadd(int nring, char *chem, int **rjs, int rjs_size);

/* Strip ring-joining numerical characters from chem strings */
void rjsrm(int nring, char *chem, int **rjs, int rjs_size);

#endif /* RJTOOL_H */
