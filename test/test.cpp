#include <iostream>
#include <complex>

using namespace std;
int main(){
	std::complex<double>* x = (std::complex<double> *)malloc(sizeof(std::complex<double>)*3);
	x[0] = std::complex<double>(1,2);
	x[1] = std::complex<double>(2,2);
	x[2] = std::complex<double>(3,2);
	cout << 2.0*x[0] << endl;
	cout << x[1] << endl;
	cout << x[2] << endl;

	free(x);
}


