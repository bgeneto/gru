#include <stdio.h>
#include <stdlib.h>

#include "mc.h"

int main()
{
    int i;
    
    srand(time(NULL));
    
    for (i=0;i<1000;i++) printf("%d %d\n", i, sorteia_salto(2) );
    printf("\n\n");
    for (i=0;i<1000;i++) printf("%d %d\n", i, sorteia_numero_salto(2) );
}
