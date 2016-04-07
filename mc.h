#ifndef _MC_H
#define _MC_H
#endif

#ifndef _MATH_H
#include <math.h>
#endif


#define RAND (((double) rand() )/(RAND_MAX+1.0))

#define ALCANCE_INT 2 // Note que estas funções são intermediárias... quando a interação for generalizada elas desaparecem!!!
#define ALCANCE_SALTO 2
#define INTERVALO_METROPOLIS ALCANCE_SALTO + ALCANCE_INT
#define SRAND ( (float) (2.0*ALCANCE_SALTO+1.0)*RAND -ALCANCE_SALTO) // // SRAND é o salto antigo

// TODO Transferir as funções de salto de manipulação de números aleatórios para "myrand.c"
// TODO Eliminar as variávels ALCANCE_INT INTERVALO_METROPOLIS ALCANCE_SALTO

// Definição básica de um sistema 
typedef struct
{
  double e1, e2, t, energia; // Parametros do gas de rede //TODO: ESTUDAR DOXYGEN
  long int *r, 			// ponteiro para uma rede unidimensional
       L, 			// tamanho da rede
      *p,  			// ponteiro com as posicoes de cada particula na rede
      *d,  			// deslocamento individual de cada particula
       N; 			// numero de particulas
} t_sistema;

// Estrutura para coletar os dados de uma simulação
// Incluir configurações (?) - Aparentemente não é necessário!!!!!
typedef struct 
{
    long int np,                   // número de pontos
             *tempo;               // aloca np valores de tempo
    double *energia,               // para alocar np valores de energia
           *dqm;                   // para alocar np valores de deslocamento quadrático médio.
} t_dados;


// Definição básica de uma simulação
typedef struct
{
    unsigned int semente;     // semente aleatória
    
    long int pmc,       // passos de monte carlo
             imp;       // passos entre impressões
            
    
    FILE  *arq_simulacao; //*arq_inicial

    int relatorio,               // Habilita/desabilita a impressão de relatório de simulação
        imprime_configuracao;    // Habilita/desabilita a impressão da configuração do sistema
    
    t_sistema *sis;
    t_dados *dados;
    
} t_simulacao;

// Criar uma função para inicializar vetor saltos arbitrários...
long int salto[4] = {-2,-1,1,2}; 

inline long int sorteia_salto(long int alcance);
inline long int sorteia_numero_salto (long int alcance);
inline long int sorteia_particula (long int N);

inline long int sorteia_salto ( long int alcance )
{
    return salto[ sorteia_numero_salto( alcance ) ];
}

inline long int sorteia_numero_salto ( long int alcance )
{
    return (long int) floor ( (double) 2.0*alcance*RAND );
}

inline long int sorteia_particula ( long int N)
{
    return floor( (double) N*RAND ) + 1;
}