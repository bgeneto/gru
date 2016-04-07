#ifndef _MC_H
#include "mc.h"
#endif

// TODO Será que vale a pena fazer checagem de erro na alocação de memória

//Manipulação de estruturas básicas
t_simulacao *aloca_memoria_simulacao();
void inicializa_simulacao_padrao(t_simulacao *sim);

t_dados *aloca_memoria_dados(long int n);

t_sistema *aloca_memoria_sistema();
void inicializa_sistema_padrao ( t_sistema *sis );

void aloca_vetores_sistema ( t_sistema *sis );
void libera_vetores_sistema ( t_sistema *sis );

// Entrada e saída...
void gera_configuracao_aleatoria(t_sistema *sim);
void imprime_configuracao(t_sistema *sis, FILE *arq);
void imprime_dados_simulacao(t_simulacao *sim);

t_simulacao *aloca_memoria_simulacao()
{
	t_simulacao *sim;
	sim = ( t_simulacao * ) calloc( 1, sizeof ( t_simulacao ) );
	return (sim);
}

void inicializa_simulacao_padrao(t_simulacao *sim)
{	
	// arquivo de saída e relatório
	sim -> arq_simulacao = stdout;
	sim -> relatorio = 0;
	
	// Configurações padrão da simulação
	sim -> semente = time(NULL);
	srand(sim -> semente);
	
	sim -> pmc = 10000;
    sim -> imp = 10;
    sim -> imprime_configuracao = 0; 
}

t_sistema *aloca_memoria_sistema()
{
	t_sistema *sis;
	// alocação de memória
	sis = ( t_sistema * )   calloc( 1, sizeof ( t_sistema ) );
	
	return (sis);
}

void inicializa_sistema_padrao ( t_sistema *sis )
{
	// Configurações padrão do sistema - Gas de rede apenas com interações atrativas entre primeiros vizinhos
	sis->e1 = -1.0;  // TODO Checar notação científica -1E0 ?
	sis->e2 = 0.0;
	sis->t  = 1.0;
	sis->L = 100;
	sis->N = 50;
}

t_dados *aloca_memoria_dados(long int n)
{
    t_dados *dados;
    
    dados          = ( t_dados * ) calloc( 1, sizeof ( t_dados ) );    
    
    dados->np      = n;
    dados->tempo   = (long int *) calloc( n, sizeof (long int) );
    dados->energia = (double *) calloc( n, sizeof (double) );
    dados->dqm     = (double *) calloc( n, sizeof (double) );
    
    return (dados);    
}

void aloca_vetores_sistema ( t_sistema *sis )
{
	sis->r = (long int *) calloc( sis->L + 1, sizeof( long int ) ); // Checar se dá para trocar o +1 por 0 e o +2 por +1 (+1 é necessário porque o indice 0 refere-se aos buracos), neste código. 
	sis->p = (long int *) calloc( sis->N + 2, sizeof( long int ) );
	sis->d = (long int *) calloc( sis->N + 2, sizeof( long int ) ); // Este é o único vetor que realmente precisa ser ( long int *). Checar se é possível converter os outros para (int).
}

void libera_vetores_sistema ( t_sistema *sis )
{
	free(sis->r);
	free(sis->p);
	free(sis->d);
}

// Note que esta função é extremamente não otimizada... para altas densidades ela roda for - ever
void gera_configuracao_aleatoria(t_sistema *sis)
{
	long int i, x;

	for ( i=1; i<=sis->N; i++ )
	{
		do
		{
			x = floor( ( double ) sis->L*RAND ); // transformar em MACRO ou em função inline!!!
		} while ( sis->r[x] );
		
		sis->r[x] = i;
		sis->p[i] = x;
	}
}

// Impressão de dados
void imprime_dados_simulacao(t_simulacao *sim)
{
    int i;
    
    for ( i=0; i < sim->dados->np; i++ )
    {
        fprintf(sim->arq_simulacao,"%ld %.6f %.6f \n", sim->dados->tempo[i], sim->dados->energia[i], sim->dados->dqm[i]); // TODO trocar acesso matricial por acesso via incremento de ponteiros            
        //  if ( sim->imprime_configuracao ) print_lint_vec( sim->arq_simulacao, sis->r, sis->L ); // TODO Reimplementar impressão de configuração
    }
}



//   TRECHO DE CODIGO PARA IMPORTAR CONFIGURAÇÃO INICIAL

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
