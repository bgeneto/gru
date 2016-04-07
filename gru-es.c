// GRU ENTRADA E SAÍDA.
#ifndef _MYTOOLS_C
#include "mytools.c"
#endif

#define GRU_VERSAO "2015.09.30"

// específico GRU
t_simulacao *prepara_simulacao(int argc, char *argv[]);
void le_parametros_gru(int argc, char *argv[], t_simulacao *sim);
int detecta_erro_parametros_simulacao(t_simulacao *sim);
void ajuda_gru();
void versao_gru();
void gera_relatorio(int argc, char *argv[], t_simulacao *sim);


t_simulacao *prepara_simulacao(int argc, char *argv[])
{
	t_simulacao *sim;
	int i;
	
	sim = aloca_memoria_simulacao ();
	sim->sis = aloca_memoria_sistema();
	
	inicializa_simulacao_padrao ( sim ); // Parâmetros padronizados
	inicializa_sistema_padrao ( sim->sis );
    
	le_parametros_gru ( argc, argv, sim );
    
	aloca_vetores_sistema ( sim->sis );
    
	gera_configuracao_aleatoria ( sim->sis );
	
	if ( sim-> relatorio ) gera_relatorio( argc, argv, sim );  
    
    /* ATENÇÃO: cada passo corresponde a uma tentativa de movimento por partícula, por isso o número de passos de Monte Carlo e o número de passos entre impressões são  multiplicados por N  */
    sim->dados = aloca_memoria_dados( sim->pmc / sim->imp );
    sim->pmc = sim->pmc*sim->sis->N;         
    sim->imp = sim->imp*sim->sis->N;    
	
	return (sim);
}

// TODO Utilizar a função getopts_long, da biblioteca <getopts.h> do Linux.
void le_parametros_gru(int argc, char *argv[], t_simulacao *sim)
{   
	int i;
    t_sistema *sis;
    sis = sim->sis;
	
	for ( i = 1; i < argc ; i += 2)
		switch(argv[i][1])
		{
            // Opções do sistema
			case 'L':  // tamanho da rede
				sis->L = atoi( argv[i+1] );
				break;
			case 'N': // numero de particulas na rede
				sis->N = atoi( argv[i+1] );
				break;
            case 'P': // interacao com o segundo vizinho
                sis->e2 = atof(argv[i+1]);
                break;
            case 'W': // interacao com o primeiro vizinho
                sis->e1 = atof(argv[i+1]);
                break;
            case 'T': // temperatura
                sis->t = atof( argv[i+1] );
                break;

            // Opções de simulação
            case 'S': // random semente
                sim->semente = atoi( argv[i+1] );
                srand (sim->semente);
                break;
			case 't': // tempo de simulacao
				sim->pmc = atoi( argv[i+1] );
				break;
			case 's': // tempo entre gravações
				sim->imp = atoi( argv[i+1] );
				break;
			case 'g': // habilita a impressão das configuracões
				sim->imprime_configuracao = 1; //
                i--;
				break;
            case 'a':
                sim->arq_simulacao = fopen(argv[i+1], "w");
                break;                
               
            // Opções de exibição:
			case 'r':
				sim->relatorio = 1;
				i--;
				break;	
			case 'h':
				ajuda_gru();
				exit(EXIT_SUCCESS);
				break;
			default:
                versao_gru(stderr);
				fprintf(stderr, "Erro: O parâmetro de entrada %s não existe! Encerrando o programa.\n",argv[i]);
				exit(EXIT_FAILURE);
		} 
		
		detecta_erro_parametros_simulacao ( sim ); 
}

// TODO melhorar a nomenclatura dos parâmetros
void ajuda_gru()
{
    versao_gru(stdout);
	fprintf( stdout,
	"Opções do sistema:\n"
	"  -L 99999   Tamanho da rede (Por padrão é 100).\n"
	"  -N 99999   Numero de particulas na rede (Por padrão é 50).\n"
	"  -P 9.999   Interação com segundos vizinhos (Pode ser atrativa ou repulsiva. Por padrão é nula, 0.0).\n"
	"  -W 9.999   Interação com primeiros vizinhos (Pode ser atrativa ou repulsiva. Por padrão é -1.0).\n"
	"  -T 9.999   Temperatura (Por padrão é 1.0).\n\n"
	""
	"Opções da simulação:\n"   
	"  -S 999     Semente Aleatória (Por padrão utiliza o tempo do sistema).\n"	  
	"  -t 99999   Número total de passos.\n"
	"  -s 999     Passos entre impressões (Por padrão imprime sempre). \n"
	"  -g         Habilita a impressão das configurações da rede (Por padrao é desabilitada).\n"
    "  -a x.dat   Nome do arquivo de saida.\n\n"  
    ""  
	"Opções de exibição:\n"
	"  -r         Gera um relatório no início da simulação.\n"
	"  -h         Imprime esta ajuda.\n\n");
    
    //        TODO ""  -i arq     Arquivo com a conformacao inicial ( por ser implementada ). Conformacao aleatoria por padrao.\n"
}

void versao_gru(FILE *arq)
{
    fprintf(arq,
            "# GRU - Simulação de Monte Carlo para o Gas de Rede Unidimensional. LCC/PPG-CIMA/FUP/UnB.\n"
            "# Versão: %s\n\n", GRU_VERSAO);
//            "Autor: Marco A. A. Barbosa - aureliobarbosa@gmail.com\n\n"  );   // Em breve teremos mais um autor!
}

int detecta_erro_parametros_simulacao(t_simulacao *sim)
{
	if ( sim->imp > sim->pmc ) 
	{
        versao_gru(stderr);
		fprintf(stderr,"Erro: Tempo de gravação maior do que o tempo de impressão!\n");
		exit(EXIT_FAILURE);
	}
	
	if (sim->sis->N >= sim->sis->L+1) 
	{
        versao_gru(stderr);
		fprintf(stderr, "Erro: O número de partículas é maior ou igual ao tamanho da rede!\n");
		exit(EXIT_FAILURE);
	}
	
	return 1;
}

// TODO Incluir data, separa relatório inicial e final... indicar tempo de execução!
void gera_relatorio(int argc, char *argv[], t_simulacao *sim)
{
	int i;  
	
    versao_gru(sim->arq_simulacao);
	fprintf(sim->arq_simulacao,	"# Parâmetros da simulação:\n# ");
	
	print_program_opts(sim->arq_simulacao,argc,argv); // REIMPRIME OS PARAMETROS DE SIMULAÇÃO
	
	fprintf(sim->arq_simulacao,"\n\n");
}