#ifndef _MC_H
#include "mc.h"
#endif

//Cálculos no sistema
double energia_particula( t_sistema *sis, long int x );
double energia_local( t_sistema *sis, long int x );
void calcula_energia_inicial(t_sistema *sis);
double calcula_deslocamento_quadratico(t_sistema *sis);

// double energia_particula( t_sistema *sis, long int part ) // CHECAR O QUE DEU ERRADO AQUI!!!
// // TODO Generalizar a interação com um vetor de interações... 
// // TODO Checar se o algoritmo nessa versão não atrapalha a paralelização em CUDA!
// {
// 	double ene = 0.0;
// 	long int x, xv, i, sinal;
// 	
//     x = sis->p[part];
//     
//     for (sinal=-1; sinal<2;sinal+=2) // Caminhar para a esquerda (-1) ou para a direita (+1)
//     {
//         i = 0;
//         do 
//         {
//             i++;
//             xv =  (x + sinal*i + sis->L )%sis->L;            
//         } while ( !sis->r[xv] ); // Possivelmente o while é problemático em CUDA...
//         
//         if ( i == 1 ) ene +=  sis->e1; // Esse trecho pode ser trocado por uma linha se for implementado um vetor de interações...
//         else 
//         {
//             if ( i == 2 ) ene += sis->e2;
//         }        
//     }    
// 	return ene;
// }

//VERSÃO ANTIGA, QUE DÁ CERTO!!!!
double energia_particula( t_sistema *sis, long int x )
// mudar o parametro de entrada... pedir o número da partícula, ao invés de pedir a posição dela
// Generalizar a interação com um vetor de interações... é só retirar o for, colocando um while... 
{
    double ene = 0.0;
    long int xv, i;
    
//    x = sis->p[part];
    
    for ( i=-1; i<2; i++ )
    {
        if (!i) continue;
        xv =  (x + i + sis->L )%sis->L;
        if ( sis->r[xv] )
            ene +=  sis->e1;
        else if ( sis->e2 != 0.000000000 ) // se o sitio estiver vazio passa para o segundo vizinho, se houver interação
        {
            xv =  (x + 2*i + sis->L )%sis->L;
            if (sis->r[xv]) ene += sis->e2;
        }
    }
    
    return ene;
}

// Energia na vizinhança da partícula no sítio x, considerando a região onde a partícula pode interagir após o salto...
// TODO CUDA: Esta função pode ser modificada para fazer uma leitura linear, calculando a energia em torno de uma vizinhança
double energia_local( t_sistema *sis, long int x ) // Checar se não seria possível reduzir a vizinhança em que é necessário calcular as energias... para otimizar o código (?).
{
	double ene = 0.0;
	long int i, xv;
	for ( i=-ALCANCE_SALTO-ALCANCE_INT; i<ALCANCE_SALTO+ALCANCE_INT+1; i++ )
	{
		xv =  (x + i + sis->L )%sis->L;
		if ( sis->r[xv] ) ene += energia_particula(sis,xv)/2.0;
	}
	return ene;
}

void calcula_energia_inicial( t_sistema *sis )
{
	int part;
	sis->energia = 0.0;
	
	for ( part = 1; part < (sis->N+1); part++ ) sis->energia += energia_particula( sis, sis->p[part] )/2.0;
}


// TODO A princípio está tudo ok, mas vale a pena checar as condições para o cálculo do deslocamento quadrático médio explodir!!!!
double calcula_deslocamento_quadratico(t_sistema *sis)
{
    int i, x2;
    
    x2 = 0; // Deslocamento quadrático médio
    for ( i=1; i < sis->N+1 ; i++ )
    {
        x2 += sis->d[i]*sis->d[i]; // será que nao estamos tendo erros numericos explosões aqui??
    }  
    return (double) x2 / sis->N;
  
}
