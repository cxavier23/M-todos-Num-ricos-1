import java.lang.Math;
//author: Christopher Xavier Sanchez Duran
//UAMC

public class ProyectoMN1
{
	public static double volumenAbrevadero(double L, double r, double h)
	{
		return L*(Math.pow(r,2) * (Math.acos(h/r) - h*Math.sqrt(Math.pow(r,2)-Math.pow(h,2))));
	}
	
	public static double gAuxAbrevadero(double L, double r, double h, double V)
	{
		return r * Math.cos((V/L + h*Math.sqrt(Math.pow(r,2)-Math.pow(h,2))) / Math.pow(r,2));
	}
	
	public static void puntoFijoAbrevadero(double L, double r, double h0, double V, double tolerancia)
	{
		System.out.println("\n\n\t h[0] = \t" + h0);
		PF_abrevaderoAUX(L, r, h0, V, tolerancia, 0, 100*tolerancia);
	}
	
	public static void PF_abrevaderoAUX(double L, double r, double h0, double V, double tol, int i, double err)
	{
		if(i<1000 && err>tol)
		{
			i++;
			double h1 = gAuxAbrevadero(L, r, h0, V);
			err = Math.abs((h1-h0)/h1);
			System.out.println("\n\t h" + "[" + i + "] = \t" + h1);
			PF_abrevaderoAUX(L, r, h1, V, tol, i, err);
		}
	}
	
	
	public static void main(String [] args)
	{
		double L = 10, r = 1, V = 12.4, tol = 0.001, h0 = 0;
		puntoFijoAbrevadero(L, r, h0, V, tol);
	}
}
//author: Christopher Xavier Sanchez Duran
// Github: https://github.com/cxavier23
// Portafolio: https://cxavier23.github.io/
