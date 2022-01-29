#include <iostream>
#include <cmath>

using namespace std;

double PolyFunction(double root) {
	return powl(root, 3) - 2 * powl(root, 2) - 4 * root + 7;
}
double TranscendentalFunction(double root) {
	double e = 2.718281828;
	return powl(e, powl(root, 2)) - 1/(root - 1);
}

double PolyDerivative(double root) {
	return 3 * powl(root, 2) - 4 * root - 4;
}

double PolySecDerivative(double root) {
	return 6 * root - 4;
}

double TranscendentalDerivative(double root) {
	double e = 2.718281828;
	return 2 * root * powl(e, powl(root, 2)) + 1 / (powl(root - 1, 2));
}

double TranscendentalSecDerivative(double root) {
	double e = 2.718281828;
	return 2 * powl(e, powl(root, 2)) + 4 * powl(root, 2) * powl(e, powl(root, 2)) - 2 / powl((root - 1), 3);
}


double BisectionMethod(double(*Function)(double root), double error, double a, double b) {

	if (!(Function(a) * Function(b) < 0)) {
		return INFINITY;
	}

	double Fa = Function(a);
	double c = (a + b) / 2;
	double Fc = Function(c);
	
	while (abs(Fc) > error) {
        
		if (Fa * Fc > 0) {
			a = c;
		}
		else if (Fa * Fc < 0){
			b = c;
		}

		c = (a + b) / 2;
		cout << c << endl;
		Fc = Function(c);
	}

	return c;
}



double NewtonMethodWithSec(double(*Function)(double root), double(*Derivative)(double root), double(*SecDerivative)(double root), double error, double a, double b) {
	double approximateRoot;
	if (SecDerivative(a) < 0) {
		approximateRoot = a;
	}
	else if (SecDerivative(a) > 0){
		approximateRoot = b;
	}
	double FPrevC = Function(approximateRoot);
	double FDPrevC = Derivative(approximateRoot);
	double prevC = approximateRoot;
	double c = prevC - FPrevC / FDPrevC;
	double Fc = Function(c);

	while (abs(Fc) > error) { //abs(c+1  -  c) < error

		FPrevC = Function(c);
		FDPrevC = Derivative(c);
		prevC = c;
		c = prevC - FPrevC / FDPrevC;
		cout << c << endl;
		Fc = Function(c);
	}
	return c;
}

int main(){

	double error = 0.001;
	
	cout.precision(5);
	
	//negative roots
	//cout << NewtonMethodWithSec(PolyFunction, PolyDerivative, PolySecDerivative, error, -4, -0.6) << endl;
	
	//positive roots
    //cout << NewtonMethodWithSec(PolyFunction, PolyDerivative, PolySecDerivative, error, 0.6, 2) << endl;
	//cout << NewtonMethodWithSec(PolyFunction, PolyDerivative, PolySecDerivative, error, 2, 5) << endl;
	
	//-------------------
	//cout << NewtonMethodWithSec(TranscendentalFunction, TranscendentalDerivative, TranscendentalSecDerivative, error, 1.1, 1.29) << endl;
    
	cout << BisectionMethod(PolyFunction, error, -4, -0.6) << endl;
	cout << "------------" << endl;
	cout << BisectionMethod(PolyFunction, error, 0.6, 2) << endl;
	cout << "------------" << endl;
	cout << BisectionMethod(PolyFunction, error, 2, 5) << endl;
	cout << "------------" << endl;
	cout << BisectionMethod(TranscendentalFunction, error, 1.1, 1.29) << endl;

	return 0;
}

