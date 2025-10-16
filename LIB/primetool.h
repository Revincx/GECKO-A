/* Prime tool module for GECKO-A */
#ifndef PRIMETOOL_H
#define PRIMETOOL_H

/* Update rank of nodes in chem. formula by finding and ranking 
 * the product of adjacent primes of input rank
 */
void primes(int nca, int **bond, int *rank, int bond_size);

/* Rank numerical priorities */
void sortrank(int nca, int **prank, int *rank, int prank_rows, int prank_cols);

#endif /* PRIMETOOL_H */
