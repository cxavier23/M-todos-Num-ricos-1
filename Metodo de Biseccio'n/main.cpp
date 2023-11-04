#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//author: Christopher Xavier Sanchez Duran
//UAMC

//Define una funcion donde encuentra la ecuacion a resolver con el
//Metodo de Biseccion
double ecuacion(double x)
{
	double f;
	f = -8*exp((1-x))+(7/x);
	return f;

  }
//Define los intervalos dados dependiento de la exactitud que quieras definirlo
int main() {
	int i = 1;
	double intervaloInicial, intervaloFinal, pM, funcion;

	intervaloInicial =  1.201538;
	intervaloFinal = 2.0000;
  printf("El intervalo Inicial es: %lf\n\n", intervaloInicial);
  printf("El intervalo final es: %lf\n\n", intervaloFinal);
	pM = (intervaloInicial + intervaloFinal)/ 2;
	funcion = ecuacion(pM);

	// tolerancia
	double tolerancia =  0.00002;
	double nuevoIntervalo, to, fn, BeseccionM;
	nuevoIntervalo = intervaloInicial;
	to = intervaloFinal;
	fn = funcion;
	BeseccionM = pM;
	//Imprime nustros datos datos al intervalos a partir de la solucion dada
	while (abs(to-nuevoIntervalo) > tolerancia)
	{
		printf("n: %i\n|nuevoIntervalo: %.5f|\n|to: %.5f, BeseccionM: %.5f, funcion(BeseccionM): %.15f|\n", i , nuevoIntervalo, to, BeseccionM, fn);
		
		if(ecuacion(nuevoIntervalo)*ecuacion(BeseccionM) < 0){
			to = BeseccionM;
		}
		else{
			nuevoIntervalo = BeseccionM;
		}
		BeseccionM = (nuevoIntervalo + to)/2;
		fn = ecuacion(BeseccionM);
		i++;
	}
	return 0;
}
// author: Christopher Xavier Sanchez Duran
// Github: https://github.com/cxavier23
// Portafolio: https://cxavier23.github.io/

