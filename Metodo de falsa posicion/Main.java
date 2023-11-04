import java.util.Scanner;
import java.lang.Math;
import java.text.DecimalFormat;
 //author: Christopher Xavier Sanchez Duran
//UAMC

public class Main
{
	//instrucciones auxiliares
	private static final double pi = Math.PI, L = 10, r = 1, V = 12.4; //constantes
	private static final int MAX_IT = 30; //maximo de iteraciones
	private static final String [] funciones = {"10(arccos(x) - x*sqrt(1 - x^2)) - 12.4", "2x*cos(2x) - (x+1)^2"};
	private static final String [] metodos = {"punto fijo", "falsa posicion"};
	private static final Scanner scan = new Scanner(System.in); 
	private static double sqrt(double x) {return Math.sqrt(x);}	
	private static String round(double x, double tol)
	{
		String numDecimales = "0.";
		int maxDecimales = (int) Math.ceil(-Math.log10(tol)+2);
		for(int i = 0; i < maxDecimales; i++)  {numDecimales += "0";}
		DecimalFormat df = new DecimalFormat(numDecimales);
	
		return df.format(x);
	}
		
		
	//instrucciones metodos numericos
	private static double f(double x, int preg) //funciones originales
	{
		if(preg == 1)
		{
			return L*(r*r * Math.acos(x/r) - x*sqrt(r*r - x*x)) - V;
		}
		else //if(preg == 2)
		{
			return 2*x*Math.cos(2*x) - (x+1)*(x+1);
		}
	}
	
	private static double g(double x0, double a, double b, int preg, int met) //funciones auxiliares
	{
		if(met == 1) //punto fijo
		{
			if(preg == 1)
			{
				return sqrt(x0*(r * Math.cos((V/L + x0*sqrt(r*r - x0*x0))/(r*r))));
			}
			else //if(preg == 2)
			{
				return 0.5 * (Math.acos((x0+1)*(x0+1) / (2*x0)) - 2*pi);
			}
		}
		else //if(met == 2)
		{
			double fb = f(b,preg), fa = f(a,preg);
			return b - fb * ((b-a) / (fb-fa)); //falsa posicion
		}
	}

	private static void primera_Iteracion(double a, double b, double tol, int preg, int met)
	{
		System.out.println("     # Itera. \t\t    Aprox.  \t\t   Error\n");
		int i = 0;
		double b0 = b, a0 = a, x0 = b;
		
		if(met == 2) 
		{
			x0 = g(x0, a0, b0, preg, met); //primera iteracion falsa posicion
			
			if(f(x0,preg)*f(b,preg) < 0)
			{
				a0 = x0;
			}
			else //if(f(a,preg)*f(x,preg) < 0)
			{
				b0 = x0;
			}
		}

		System.out.println("\t " + i + " \t\t  " + round(x0,tol));
		siguiente_Iteracion(x0, a0, b0, tol, preg, met, i+1, tol+1); //recursion 
	}
	
		private static void siguiente_Iteracion(double x0, double a, double b, double tol, int preg, int met, int i, double err)
	{
		if(i < MAX_IT && err >= tol)
		{
			double x = g(x0, a, b, preg, met); //nueva aproximacion
			err = Math.abs((x - x0)/ x); //error relativo 
			System.out.println("\t " + i + " \t\t  " + round(x,tol) + "\t\t  " + round(err,tol));
			
			if(met == 2 && f(x,preg)*f(b,preg) < 0)
			{
				siguiente_Iteracion(x, x, b, tol, preg, met, i+1, err);
			}
			else //if(met == 1 || f(a,preg)*f(x,preg) < 0)
			{
				siguiente_Iteracion(x, a, x, tol, preg, met, i+1, err);
			}
		}
		else if(i >= MAX_IT) //casos base
		{
			System.out.println("\n No se alcanzo la tolerancia en " + MAX_IT + " iteraciones.");
		}
		else //if(err < tol)
		{
			System.out.println("\n\n La aproximacion encontrada a una raiz de f(x)" +
				" con una tolerancia de " + tol +  " es " + round(x0,tol) + ".");
		}
	}	
	
	
	//instrucciones menu interactivo
	public static void select_Preg() 
	{
		int preg = -1;
	
		while(preg!=1 && preg!=2)
		{	
			System.out.println("\n\n Seleccione una funcion:\n");
			System.out.println("\t (1) \t f(x) = " + funciones[0] + "\n\t (2) \t f(x) = " + funciones[1]);
			System.out.println("\t (0) \t Salir del programa\n");
			System.out.println("\n Inserte una opcion:");
			preg = scan.nextInt();
			
			if(preg == 0)
			{
				return;
			}
			else if(preg!=1 && preg!=2)
			{
				System.out.println("\n\n ERROR: Opcion invalida. Intente de nuevo.\n");
			}
			else //if(preg == 1 || preg == 2)
			{
				System.out.println("\n\n Ha seleccionado la funcion:\t f(x) = " + funciones[preg-1] + "\n\n");
			}
		}
		
		int met = select_Metodo();
		
		if(met != 0)
		{
			menu(preg, met);
		}
		else //if(met == 0)
		{
			select_Preg();
		}
	}
	
	private static int select_Metodo()
	{
		int met = -1;
	
		while(met!=0 && met!=1 && met!=2)
		{
			System.out.println("\n Seleccione el metodo a emplear:\n");
			System.out.println("\t (1) \t Metodo de " + metodos[0] + "."); //punto fijo
			System.out.println("\t (2) \t Metodo de " + metodos[1] + "."); //falsa posicion
			System.out.println("\t (0) \t Volver al menu anterior.\n");
			System.out.println("\n Inserte una opcion:");
			met = scan.nextInt();
			
			if(met!=0 && met!=1 && met!=2)
			{
				System.out.println("\n\n ERROR: Opcion invalida. Intente de nuevo.\n");
			}
			else if(met != 0)
			{
				System.out.println("\n\n Ha seleccionado el metodo de " + metodos[met-1] + ".\n");
			}
		}
		
		return met;
	}
	
	private static void menu(int preg, int met)
	{		
		double a = -1, b = -1, fa = -1, fb = -1;
		
		if(met == 1)
		{	
			System.out.println("\n\n Inserte una aproximacion inicial:");
			b = scan.nextDouble(); //preg 1: (0, 0.2),   preg 2: [-2.5,-2]
		}
		else //if(met == 2)
		{
			while(fa*fb > 0)
			{
				System.out.println("\n\n Inserte el extremo izquierdo del intervalo:");
				a = scan.nextDouble(); //preg 1: [-1, -0.1],    preg 2:  [-4, -2.2]
				fa = f(a,preg);
				
				if(fa == 0) //si x=a es la solucion
				{
					b = a;
					fb = fa;
				}
				else //if(fa != 0)
				{
					System.out.println("\n Inserte el extremo derecho del intervalo:");
					b = scan.nextDouble(); //preg 1:[0.2, 1],    preg 2: [-2, -1]
					fb = f(b,preg);
				}
				
				if(fa*fb > 0) //si no hay cambio de signo
				{
					System.out.println("\n\n ERROR: El intervalo no satisface la hipotesis" +
							" del teorema de Bolzano. Intente de nuevo.\n");
				}
			}
		}
		
		if(fb == 0) //si x=a o x=b es la solucion
		{
			System.out.println("\n\n\n La solucion a la ecuacion " + funciones[preg-1] + " = 0  es:\t" + b);
			return;
		}
		
		System.out.println("\n\n Inserte una tolerancia:");
		double tol = scan.nextDouble(); //Preg 1: tol <= 0.01,   Preg 4: tol <= 0.00001
		tol = Math.abs(tol); //ignora el signo
		
		System.out.println("\n\n\n\n\n Funcion:\t\tf(x) = " + funciones[preg-1] + "\n Metodo:\t\t" + metodos[met-1]);
		if(met == 1)
		{
			System.out.println(" Aprox. Inicial:\t" + b);
		}
		else //if(met == 2)
		{
			System.out.println(" Intervalo:\t\t" + "(" + a + ", " + b + ")");
		}
		System.out.println(" Tolerancia:\t\t" + tol + "\n\n");
		
		primera_Iteracion(a, b, tol, preg, met); 
		select_Preg();
	}
	
	private static void bienvenida()
	{
		System.out.println("\n\n\n Este programa implementa los metodos numericos de " + metodos[0] + 
				" y de " + metodos[1] +  "\n para aproximar raices de funciones. \n" +  
				"\n\n\n Este programa fue creado por: \n\n \t Christopher Xavier Sanchez Duran " + 
				"Presione ENTER para continuar...");
		scan.nextLine();
	}
	
	public static void main(String [] args)
	{	
		bienvenida();
		select_Preg();
	}
}
// author: Christopher Xavier Sanchez Duran
// Github: https://github.com/cxavier23
// Portafolio: https://cxavier23.github.io/
