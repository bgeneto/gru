/* GRU (Gás de rede unidimensional)
 * Autor: MARCO AURÉLIO A. BARBOSA
 * e-mail: aureliobarbosa@gmail.com
 *
  * FUNCIONALIDADE
 * O programa realiza simulações de Monte Carlo usando o critério de Metrópolis em uma rede unidimensional
 * com partículas com caroço duro e interações de primeiros e segundos vizinhos.
 *
 * DETALHES DO CÓDIGO
 * - Note que as posições na rede estão indexadas de 0 a L-1, mas que as partículas estão indexadas de 1 a N
 *
 * * Checar o arquivo detalhes-codigo.txt para mais informações
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "mytools.c"

#include "gru.h"

void prepara_simulacao(int argc, char *argv[]);
void gera_relatorio(int argc, char *argv[]);

void gera_configuracao_aleatoria();
int le_configuracao_inicial( long int *s, char *nome_arq_inicial);  // fazer esta funcao
void imprime_configuracao(long int *s, long int l, FILE *arq);

void  calcula_dados_sistema(); //

double energia_particula( long int x );
double energia_local( long int x );

int main (int argc, char *argv[])
{
  t_sistema *sis;
  t_simulacao *sim;
  
  sis = ( t_sistema * )   calloc( 1, sizeof ( t_sistema ) ) ;
  sim = ( t_simulacao * ) calloc( 1, sizeof ( t_simulacao ) ) ;
  
  double enei, enef, dene, x2d;
  long int i, part, x, xv, eta, etav, salto, x2;

  parametros( argc, argv );

  // inicio da simulacao
  // calculo dos dados iniciais
  calcula_dados_sistema();

  for (sim.pmc=0;sim.pmc<sim.PMC;sim.pmc++)
  {
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

      xv =  (x + salto + sis.L )%sis.L;

      if ( !sis.r[xv] ) // NÃO COLIDIU
      {
	  // criterio de Metropolis
	  enei = energia_local(x);
	  sis.r[xv] = sis.r[x];
	  sis.r[x]  = 0;
	  enef = energia_local(x);
	  dene = enef - enei;

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
      }  


    } // fim do laco de impressao

    x2 = 0; // Deslocamento quadrático médio
    for ( i=1;i<sis.N+1;i++)
    {
	x2 += sis.d[i]*sis.d[i];
    }
    x2d = (double) x2/sis.N;
    fprintf(sim.arq_simulacao,"%ld %.6f %.6f \n", sim.pmc*sim.IMP, sis.energia, x2d);

    if ( sim.PGRAVA && !( sim.pmc % (sim.PGRAVA+1) ) ) imprime_configuracao( sis.r, sis.L, sim.arq_simulacao); // ver impressao de arquivo

  } // fim do laco da corrida
}

void parametros(int argc, char *argv[])
{
  int i, habilita_relatorio=0;
  
  char *nome_arq_simulacao = NULL;  
  sim.arq_simulacao = stdout;
  
  sis.e1 = -1.00000000000;
  sis.e2 = 0.00000000000;
  sis.t  = 1.00000000000;
  sis.L = 100;
  sis.N = 50;  
  
  sim.semente = time(NULL);
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
      case 'S': // random semente
        sim.semente = atoi( argv[i+1] );
        break;
      case 'r':
	habilita_relatorio = 1;
	i--;
	break;	
//       case 'i':  
// 	indicou_configuracao_inicial=1; // arquivo de inicialização   - ESTA FUNÇÃO NÃO FOI IMPLEMENTADA NO PROGRAMA!!!!!!!!!!!!!!!!!!
//         nome_arq_inicial  = ( char * ) calloc(strlen(argv[i+1]),sizeof(char) );
//         strcpy( nome_arq_inicial, argv[i+1] ); 
// 	
//         break;
	  
      case 'a':
        nome_arq_simulacao  = ( char * ) calloc(strlen(argv[i+1]),sizeof(char) );
        strcpy( nome_arq_simulacao, argv[i+1] );
        sim.arq_simulacao  = fopen(nome_arq_simulacao, "w");
        break;
      case 'P': // interacao com o segundo vizinho
        sis.e2 = atof(argv[i+1]);
        break;
      case 'W': // interacao com o primeiro vizinho
        sis.e1 = atof(argv[i+1]);
        break;
    case 'h':
      printf("\nMonte Carlo - Canonico - gas de rede unidimensional (Última atualização: 16/10/2015)\n"
          "Autor: Marco A. A. Barbosa - aureliobarbosa@gmail.com\n\n"
	 
          "Opções do sistema:\n"
          "  -l 99999   Tamanho da rede (Por padrao é 100).\n"
          "  -N 99999   Numero de particulas na rede .\n"
          "  -P 9.999   Interacao com segundos vizinhos (Pode ser atrativa ou repulsiva. Por padrao é nula, 0.0).\n"
          "  -W 9.999   Interacao com primeiros vizinhos (Pode ser atrativa ou repulsiva. Por padrao é -1.0).\n"
          "  -T 9.999   Temperatura (Por padrão é 1.0).\n\n"
          ""
          "Opções da simulação:\n"   
          "  -S 999     Semente Aleatoria (Por padrao utiliza o tempo do sistema).\n \n"	  
          "  -t 99999   Numero total de passos.\n"
          "  -s 999     Passos entre impressoes (Por padrão imprime sempre). \n"
          "  -g 999     Passos entre a gravacao das configuracao (Por padrao está desabilitado).\n\n"
          "" 
	  "Opções de exibição:"
          "  -a x.dat   Nome do arquivo de saida.\n"
	  "  -r         Gera um relatório no início da simulação.\n"
//          "  -i arq     Arquivo com a conformacao inicial ( por ser implementada ). Conformacao aleatoria por padrao.\n"
          "  -h         Imprime esta ajuda.\n\n");
      exit(1);
      break;
      default:
	printf("Erro: O parâmetro de entrada %s %s não existe! Encerrando o programa.\n",argv[i],argv[i+1]);
	exit(1);
  }

  if (sim.IMP > sim.PMC || ( (sim.IMP > sim.PGRAVA || sim.PGRAVA > sim.PMC) && sim.PGRAVA ) )
  {
    printf("Há algum problema com os tempos de gravacao\n"
        "tempo de simulacao: %ld\n"
        "tempo de gravacao das configuracoes: %ld\n"
        "tempo de impressao %ld\n", sim.PMC, sim.PGRAVA, sim.IMP);
    exit(1);
  }

  sim.PMC /= sim.IMP;
  if (sim.PGRAVA) sim.PGRAVA /= sim.IMP;

  srand (sim.semente);

  sis.r = (long int *) calloc( sis.L + 1, sizeof( long int ) );
  sis.p = (long int *) calloc( sis.N + 2, sizeof( long int ) );
  sis.d = (long int *) calloc( sis.N + 2, sizeof( long int ) );
  
  gera_configuracao_aleatoria( );
  if ( habilita_relatorio ) gera_relatorio( argc, argv );  
}

void gera_relatorio(int argc, char *argv[])
{
  int i;  fprintf(sim.arq_simulacao,"\n\n# Relatório de Simulação.\n");  
  fprintf(sim.arq_simulacao,"# GRU - Gás de Rede Unidimensional\n");
  fprintf(sim.arq_simulacao,"# Parâmetros da simulação:\n# ");  
  print_program_opts(sim.arq_simulacao,argc,argv);
  fprintf(sim.arq_simulacao,"\n\n");
}


void gera_configuracao_aleatoria()
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

void imprime_configuracao(long int *s, long int l, FILE *arq)
{
  long int i;
  for ( i=0; i<l; i++ ) fprintf(arq,"%ld ",s[i]);
  fprintf(arq,"\n");
}

// implementar
int le_configuracao_inicial( long int *s, char *nome_arq_inicial)
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

//   TRECHO DE CODIGO PARA GERAR IMPORTAR CONFIGURAÇÃO INICIAL

//   if ( !indicou_configuracao_inicial )
//    gera_configuracao_aleatoria( );
//   else
//   {
//     if ( !le_configuracao_inicial(sis.r, nome_arq_inicial) )
//     {
//       printf("O arquivo com a configuracao inicial: %s, esta com algum problema!\n", nome_arq_inicial);
//       exit(0);
//     }  
//   }
