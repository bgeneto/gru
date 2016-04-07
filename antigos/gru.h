#define RAND (((double) rand())/(RAND_MAX+1.0))

#define ALCANCE_INT 2
#define ALCANCE_SALTO 2
#define INTERVALO_METROPOLIS ALCANCE_SALTO + ALCANCE_INT
// Neste programa SRAND é o salto
#define SRAND ( (float) (2.0*ALCANCE_SALTO+1.0)*RAND -ALCANCE_SALTO)

typedef struct
{
  double e1, e2, t, energia; // Parametros do gas de rede
  long int *r, // ponteiro para uma rede unidimensional
       L, // tamanho da rede
      *p,  // ponteiro com as posicoes de cada particula na rede
      *d,  // memória do deslocamento total de cada particula
       N; // numero de particulas

} t_sistema;

typedef struct
{
  long int PMC,      // passos de monte carlo
           IMP,      // passos para impressao
           PGRAVA;   // passos para gravacao
  unsigned int semente;     // semente aleatória

  long int pmc,      // variaveis da simulacao
           imp;

  FILE *arq_inicial, *arq_simulacao;  
  int argc, *argv;
  
  t_sistema *sis;
  
} t_simulacao;
