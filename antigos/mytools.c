#ifndef _STDIO_H
#include <stdio.h>
#endif

// Print the program call to file fd
void print_program_opts(FILE *fd, int argc, char *argv[])
{
  int i;
  for ( i = 0; i < argc ; i++ )
    fprintf(fd,"%s ",argv[i]);  
}

// int main (int argc, char *argv[])
// {
//   print_program_opts ( stdout, argc, argv );  
// }