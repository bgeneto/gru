#ifndef _STDIO_H
#include <stdio.h>
#endif

#ifndef _MYTOOLS_C
#define _MYTOOLS_C
#endif

void print_program_opts(FILE *fd, int argc, char *argv[]);
void print_lint_vec(FILE *fd, long int *s, long int l);

// Print the program call to file fd
void print_program_opts(FILE *fd, int argc, char *argv[])
{
  int i;
  for ( i = 0; i < argc ; i++ )
    fprintf(fd,"%s ",argv[i]);  
}

// Print a list of long intergers inside file fd
void print_lint_vec(FILE *fd, long int *s, long int l)
{
  long int i;
  for ( i=0; i<l; i++ ) fprintf(fd,"%ld ",s[i]);
  fprintf(fd,"\n");
}


/* FROM http://rosettacode.org/wiki/Knuth_shuffle#C
 * RANDOMIZES A SET OF NUMBERS...
#include <stdlib.h>
#include <string.h>
 
int rrand(int m)
{
  return (int)((double)m * ( rand() / (RAND_MAX+1.0) ));
}
 
#define BYTE(X) ((unsigned char *)(X)) 
void shuffle(void *obj, size_t nmemb, size_t size)
{
  void *temp = malloc(size);
  size_t n = nmemb;
  while ( n > 1 ) {
    size_t k = rrand(n--);
    memcpy(temp, BYTE(obj) + n*size, size);
    memcpy(BYTE(obj) + n*size, BYTE(obj) + k*size, size);
    memcpy(BYTE(obj) + k*size, temp, size);
  }
  free(temp);
} */

// int le_configuracao_inicial( long int *s, char *nome_arq_inicial);  // fazer esta funcao - read_integer_seq
// {
//   return 0;
// }

// int main (int argc, char *argv[])
// {
//   print_program_opts ( stdout, argc, argv );  
// }