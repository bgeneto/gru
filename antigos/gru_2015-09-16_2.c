/* GRU (Gás de rede unidimensional)
 * Autor: MARCO AURÉLIO A. BARBOSA
 * e-mail: aureliobarbosa@gmail.com
 * 
 * FUNCIONALIDADE
 * O programa realiza simulações de Monte Carlo usando o critério de Metrópolis em uma rede unidimensional 
 * com partículas com caroço duro e interações de primeiros e segundos vizinhos.
 *   
 * 
 * DETALHES DO CÓDIGO
 * - Note que as posições na rede estão indexadas de 0 a L-1, mas que as partículas estão indexadas de 1 a N
 * 
 * TAREFAS PENDENTES
 * - Relatorios de simulação não estão funcionando!!!!
 * - Mudar a função energia_particula() para que ela receba um número de partícula, e não uma posição na rede.
 * - Generalizar o critério de Metrópolis (Para a colaboração com o Miguel Rubi!)
 * - Rodar um número arbitrário de simulações com condições iniciais diferentes e extrair a média.
 *  
 * ERROS
 * 1) > ./gru -l 100 -N 72 -t 10000 -s 10 -T 0.4 -a gru-l100N72T0.4 
 *    > Falha de segmentação (imagem do núcleo gravada) 
 *    Este erro ocorre apenas quando é solicitado que a simulação seja gravada em um arquivo!
 *    O arquivo .rel é criado e o .dat nem mesmo chega a ser criado.
 *    # CORRIGIDO HOJE!
 * 
 * IDEIAS PARA FUTURAS IMPLEMENTAÇÕES
 * - Ler a lista de argumentos de forma mais inteligente e/ou padronizada (getopt_long)
 * - Gravar os arquivos utilizando uma nomenclatura padronizada com os dados da simulação!
 * - leitura da configuração inicial a partir de um arquivo
 * - número de n-ésimos vizinhos
 * - Generalizar o tipo de salto
 * - Generalizar o tipo de interação
 * - criar uma tabela de interaçoes
 * 
 * HISTÓRICO
 * O código vem sendo desenvolvido desde o final do meu doutorado, em 2008, e foi utilizado para rodar simulações 
 * e calcular a constante de auto-difusão de modelos com propriedades de água em uma dimensão. Os artigos foram publicados em 2011 e 2015 na revista 
 * Journal of Chemical Physics.
 * 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

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
  
  // DADOS DO SISTEMA NO INSTANTE 
  long int *nv; //  
} t_sistema;

typedef struct 
{
  long int PMC,      // passos de monte carlo
           IMP,      // passos para impressao
           PGRAVA;   // passos para gravacao
  unsigned int seed;     // semente aleatória
  
  long int pmc,      // variaveis da simulacao
           imp;	   
  double ti,    // temperatura inicial
         tf,    // temperatura final
         dt;    // variacao de temp
  char *arq_inicial, *arq_simulacao, *arq_relatorio;  // nomes dos arquivos 
 
} t_simulacao;

void parametros(int argc, char *argv[]);
void configuracao_aleatoria();
int le_configuracao_inicial( long int *s, char *arq_inicial);  // fazer esta funcao
void imprime_configuracao(long int *s, long int l, FILE *arq);

void  calcula_dados_sistema(); //


double energia_particula( long int x );
double energia_local( long int x );

t_sistema sis;
t_simulacao sim;

int main (int argc, char *argv[])
{
  FILE *f_inicial=NULL, *f_simulacao=NULL, *f_relatorio=NULL;

  double enei, enef, dene, x2d;
  long int i, part, x, xv, eta, etav, salto, x2;
 
  parametros( argc, argv );
 
  f_simulacao  = ( sim.arq_simulacao == NULL ) ? stdout : fopen(sim.arq_simulacao, "w");
  f_relatorio  = ( sim.arq_relatorio == NULL ) ? stdout : fopen(sim.arq_relatorio, "w");
 
  fprintf(f_relatorio,"  ");
  
  if ( sim.arq_inicial == NULL )
    configuracao_aleatoria( );
  else
    if ( !le_configuracao_inicial(sis.r, sim.arq_inicial) )
  {
    printf("O arquivo com a configuracao inicial: %s, esta com algum problema!\n", sim.arq_inicial);
    exit(0);
  }

  // inicio da simulacao
  // calculo dos dados iniciais
  calcula_dados_sistema();
  
  //* printf("\n");  
  // inicio da corrida
  for (sim.pmc=0;sim.pmc<sim.PMC;sim.pmc++)
  { 
      //* printf("\n* %ld\n", sim.pmc); //
    for (sim.imp=0;sim.imp<sis.N*sim.IMP;sim.imp++)
    {
      // sorteio da particula 
      part = (long int) floor( (double) sis.N*RAND )+1; 
      x = sis.p[part];
      
      // sorteio da movimentacao da rede
      do 
      { 
	    salto = (long int) floor(SRAND); 
      } while ( !salto ); 
      
      //*printf("%ld  %ld  %+ld ", sim.imp, part, salto ); // *
      
      xv =  (x + salto + sis.L )%sis.L;
      
      if ( !sis.r[xv] ) // NÃO COLIDIU
      {	        
	  // criterio de Metropolis 
	  enei = energia_local(x);
	  sis.r[xv] = sis.r[x];
	  sis.r[x]  = 0;
	  enef = energia_local(x);
	  dene = enef - enei;
	  
	  //* printf("%+2.2f %+2.2f & ", enei, enef); //*
    
	  if ( RAND > exp(-dene/sis.t) )  	// Não passou por Metrópolis
	  {
	      sis.r[x] = sis.r[xv];  		
	      sis.r[xv] = 0;	
	  } else 				// Passou por Metrópolis  -- > Atualizacao dos valores      
	  {
	      sis.p[part] = xv;                 
	      sis.d[part] += salto;
	      sis.energia += dene;     
	  }	  
      }  //* else  printf("              ",dene); //*     
      
      
	  //*imprime_configuracao(sis.r,sis.L,stdout); // ver impressao de arquivo      // *
    } // fim do laco de impressao
    
    x2 = 0; // Deslocamento quadrático médio
    for ( i=1;i<sis.N+1;i++) 
    {
	x2 += sis.d[i]*sis.d[i];
    }
    x2d = (double) x2/sis.N;
    fprintf(f_simulacao,"%ld %.6f %.6f \n", sim.pmc*sim.IMP, sis.energia, x2d);
       
    if ( sim.PGRAVA && !( sim.pmc % (sim.PGRAVA+1) ) ) imprime_configuracao( sis.r, sis.L, f_simulacao); // ver impressao de arquivo

  } // fim do laco da corrida
}

void parametros(int argc, char *argv[])
{
  int i;

  sis.e1 = -1.00000000000;
  sis.t  = 1.00000000000;
  sis.e2 = 0.00000000000;
  sis.t  = 1.00000000000; 
  sis.L = 100;
  sis.N = 50;
  
  sim.seed = time(NULL);
  sim.arq_inicial = sim.arq_simulacao = sim.arq_relatorio = NULL;
  sim.PGRAVA = 0;
  sim.IMP = 1;
    
  for ( i = 1; i < argc ; i += 2)
    switch(argv[i][1])
    {
      case 'l':  // tamanho da rede
        sis.L = atoi( argv[i+1] );
        break;
      case 'N': // numero de particulas na rede
        sis.N = atoi( argv[i+1] );
        break;        
      case 't': // tempo de simulacao
        sim.PMC = atoi( argv[i+1] );
        break;
      case 's': // grava os dados
        sim.IMP = atoi( argv[i+1] );
        break;
      case 'g': // grava as configuracões
        sim.PGRAVA = atoi( argv[i+1] );
        break;	
      case 'T': // temperatura
        sis.t = atof( argv[i+1] );
        break;
      case 'S': // random seed
        sim.seed = atoi( argv[i+1] );
        break;
      case 'i':  // arquivo de inicializaćão
        sim.arq_inicial  = ( char * ) calloc(256,sizeof(char) );
        strcpy( sim.arq_inicial, argv[i+1] ); //        arq_inicial = argv[i+1];
        break;
      case 'a':
        sim.arq_simulacao  = ( char * ) calloc(256,sizeof(char) );
        strcpy( sim.arq_simulacao, argv[i+1] );
        strcat( sim.arq_simulacao, ".dat");
        
        sim.arq_relatorio = ( char * ) calloc(256,sizeof(char) );
        strcpy( sim.arq_relatorio, argv[i+1] );	
        strcat( sim.arq_relatorio, ".rel");
        break;
      case 'P': // interacao com o segundo vizinho
        sis.e2 = atof(argv[i+1]);
        break;
      case 'W': // interacao com o primeiro vizinho
        sis.e1 = atof(argv[i+1]);
        break;
    case 'h':
      printf("\nMonte Carlo - Canonico - gas de rede unidimensional\n"
          "Marco A. A. Barbosa - aureliobarbosa@gmailcom\n"
          "Última atualização: 2015-10\n"
          "Opcoes:\n"
          "  -l 99999   Tamanho da rede.\n"
          "  -N 99999   Numero de particulas na rede.\n"          
          "  -t 99999   Numero total de passos.\n"
          "  -s 999     Passos entre impressoes. \n"
          "  -g 999     Passos entre a gravacao das configuracao. Por padrao e desabilitado.\n"
          "  -a xxxx    Nome que definira os arquivos a serem gravados. No maximo 255 caracteres!\n"
          "  -T 9.999   Temperatura.\n"
          "  -S 999     Semente Aleatoria. Por padrao utiliza o tempo do sistema para obte-la.\n"
          "  -i arq     Arquivo com a conformacao inicial( por ser implementada ). Conformacao aleatoria por padrao.\n"
          "  -P 9.999   Interacao com segundos vizinhos. Por padrao a interacao é nula.\n"
          "  -W 9.999   Interacao com primeiros vizinhos. Por padrao a interacao é e1=-1.0.\n"
          "  -h         Imprime esta ajuda.\n\n");
      exit(1);
      break;
  }
    
  if (sim.IMP > sim.PMC || ( (sim.IMP > sim.PGRAVA || sim.PGRAVA > sim.PMC) && sim.PGRAVA ) )
  {
    printf("Ha algum problema com os tempos de gravacao\n"
        "tempo de simulacao: %ld\n"
        "tempo de gravacao das configuracoes: %ld\n"
        "tempo de impressao %ld\n", sim.PMC, sim.PGRAVA, sim.IMP);
    exit(1);
  }

  sim.PMC /= sim.IMP;  
  if (sim.PGRAVA) sim.PGRAVA /= sim.IMP;
    
  srand (sim.seed);
  
  sis.r = (long int *) calloc( sis.L + 1, sizeof( long int ) );
  sis.p = (long int *) calloc( sis.N + 2, sizeof( long int ) );  
  sis.d = (long int *) calloc( sis.N + 2, sizeof( long int ) );
}

// OK
void configuracao_aleatoria()
{
  long int i, x;
  for ( i=1; i<=sis.N; i++ ) 
  {
    do 
    {
      x = floor( ( double ) sis.L*RAND );
    } while ( sis.r[x] );
    
    sis.r[x] = i;
    sis.p[i] = x;
  }
}

// OK
void imprime_configuracao(long int *s, long int l, FILE *arq)
{
  long int i;
  for ( i=0; i<l; i++ ) fprintf(arq,"%ld ",s[i]);
  fprintf(arq,"\n");
}

// implementar
int le_configuracao_inicial( long int *s, char *arq_inicial)
{
  return 0;
}

double energia_particula(long int x) 
// mudar o parametro de entrada... pedir o número da partícula, ao invés de pedir a posição dela
{
    double ene = 0.0;
    long int xv, i; 
    
	for ( i=-1; i<2; i++ )
	{
		if (!i) continue;
		xv =  (x + i + sis.L )%sis.L;
		if ( sis.r[xv] ) 
			ene +=  sis.e1;
		else if ( sis.e2 != 0.000000000 ) // se o sitio estiver vazio passa para o segundo vizinho, se houver interação
		{
			xv =  (x + 2*i + sis.L )%sis.L;
			if (sis.r[xv]) ene += sis.e2;
		} 
	}
	
    return ene;
}

// Energia na vizinhança da partícula no sítio x, considerando a região onde a partícula pode interagir após o salto...
double energia_local( long int x) // Checar se não seria possível reduzir a vizinhança em que é necessário calcular as energias... para otimizar o código (?).
{
	double ene = 0.0;
	long int i, xv;
	for ( i=-ALCANCE_SALTO-ALCANCE_INT; i<ALCANCE_SALTO+ALCANCE_INT+1; i++ )
	{
		xv =  (x + i + sis.L )%sis.L;    
		if ( sis.r[xv] ) ene += energia_particula(xv)/2.0;
	}
	return ene;
}

void calcula_dados_sistema()
{
    int part;    
    sis.energia = 0.0;
    
    for (part=1;part<sis.N+1;part++) sis.energia += energia_particula( sis.p[part] )/2.0; 
}
