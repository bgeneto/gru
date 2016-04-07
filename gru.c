/* GRU (Gás de rede unidimensional)
 * Autor do código principal: MARCO AURÉLIO A. BARBOSA (aureliobarbosa@gmail.com)
 * Paralelização: BERNHARD ENDERS
 * 
 *
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

#include "mc.h"
#include "mc-ies.c"
#include "mc-calc.c"

#include "gru-es.c" // funções de entrada e saída - inicialização do programa.

int main (int argc, char *argv[])
{
    t_simulacao *sim;
    t_sistema *sis;
    
    double enei, enef, dene;
    long int i, l, part, x, xv, salto;
    
    sim = prepara_simulacao( argc, argv );
    sis = sim->sis;
    
    // inicio da simulacao
    calcula_energia_inicial(sis);
    
    for (i=0; i < sim->pmc; i++) // corrida 
    {
        // sorteio da particula
        part = (long int) sorteia_particula(sis->N);
        x = sis->p[part]; // descobre a posição de part
        
        // sorteio do salto a ser dado
        salto = (long int) sorteia_salto(ALCANCE_SALTO);
        
        // posição do salto, com condições periódicas de contorno
        xv =  (x + salto + sis->L )%sis->L; 
        
        if ( !sis->r[xv] ) // NÃO COLIDIU
        {
            // criterio de Metropolis
            enei = energia_local( sis, x );
            sis->r[xv] = sis->r[x];
            sis->r[x]  = 0;
            enef = energia_local(sis, x);
            dene = enef - enei;
            
            if ( RAND > exp(-dene/sis->t) )  	// Não passou por Metrópolis
            {
                sis->r[x] = sis->r[xv];
                sis->r[xv] = 0;
            } 
            else 				                // Passou por Metrópolis  -- > Atualizacao dos valores
            {
                sis->p[part] = xv;
                sis->d[part] += salto;
                sis->energia += dene;
            }
        }
        
        if ( !( i % sim->imp ) ) // Condição para armazenador informações do sistema em dados.
        {
            l = i / sim->imp;

            sim->dados->tempo[l] = i/sis->N;
            sim->dados->dqm[l] = calcula_deslocamento_quadratico (sis); // TODO trocar acesso matricial por acesso via incremento de ponteiros            
            sim->dados->energia[l] = sis->energia;             
        }
        
    } // fim do laco da corrida

    imprime_dados_simulacao(sim);
    
    exit(EXIT_SUCCESS);
    
}
