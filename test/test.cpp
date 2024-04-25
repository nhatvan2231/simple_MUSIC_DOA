#include <iostream>

using namespace std;
int test(int* &x){
	x[0] = 1;
	x[1] = 2;
	return 0;
}
int main(){
	int* x = (int *)malloc(sizeof(int)*2);
	cout << x[0] << endl;
	cout << x[1] << endl;

	test(x);

	cout << x[0] << endl;
	cout << x[1] << endl;
	free(x);
}


