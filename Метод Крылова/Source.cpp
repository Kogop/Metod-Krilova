#include <iostream>



using namespace std;
int const n = 4;
double A[4][4] = { 2.2,1,0.5,2,1,1.3,2,1,0.5,2,0.5,1.6,2,1,1.6,2 };
double y_0[4] = { 1,0,0,0 };
double Y[4] = { 0,0,0,0 };
double YItoe[n][n];

int k = -1;
double q[n][n];
double x[n][n];
double q1[n][n+1];
int N = n + 1;
double p1[n];
const double eps = 0.0001, delta = 0.001, Min = -100, Max = 100, H = 0.001;
double K_Pol[n], S_korn[n];


void Gauss(int k, double Matrix[n][n + 1]) {
	if (Matrix[k][k] != 1) {
		double T = Matrix[k][k];
		for (int j = k; j < n+1; j++) {
			Matrix[k][j] = Matrix[k][j] / T;
		}
	}
	for (int i = 0; i < n; i++) {
		if ((Matrix[i][k] != 0) && (i != k)) {
			double T = Matrix[i][k];
			Matrix[i][k] = 0;
			for (int j = k + 1; j < n+1; j++) {
				Matrix[i][j] -= Matrix[k][j] * T;
			}
		}
	}
	if (k < n - 1) {
		Gauss(k + 1, Matrix);
	}
}

double Zamena(int k) {

	for (int i = 0; i < n; i++)
	{
		y_0[i] = Y[i];
	
		YItoe[i][k] = Y[i];
		Y[i] = 0;
	}
	return 1;
}

double Yitie() {

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Y[i] += A[i][j] * y_0[j];
		}
		//cout << Y[i] << "\n";
		
	}
	k++;
	Zamena(k);
	return 1;
}

double ZagatovkaDlyaGaussa() {
	for (int i = 0; i < n; i++)
	{
		q1[i][0] = YItoe[i][2];
		q1[i][1] = YItoe[i][1];
		q1[i][2] = YItoe[i][0];
		q1[i][4] = YItoe[i][3];
	}
	q1[0][3] = { 1 };
	q1[1][3] = { 0 };
	q1[2][3] = { 0 };
	q1[3][3] = { 0 };

	return 1;
}

//значение полинома в точке xf
double Pol(double* A, double xf) {
	double Sum = pow(xf, N);
	for (size_t i = 1; i <= N; i++) {
		Sum -= K_Pol[i - 1] * pow(xf, N - i);
	}
	return pow(-1, N) * Sum;
}
//производная полинома в точке xf
double dx(double xf) {
	return (Pol(K_Pol, xf + delta) - Pol(K_Pol, xf - delta)) / 2 / delta;
}

//помещаем начальное значение, получаем уточненное значение
double Nyuton(double x_0) {
	double Mod = 1, xi = 0;
	while (eps < Mod) {
		xi = x_0 - Pol(K_Pol, x_0) / dx(x_0);
		Mod = fabs(xi - x_0);
		x_0 = xi;
	}
	return xi;
}

//отделение корней, помещает примерное значение корней
void Search(double* A) {
	int i = 0;
	for (double t = Min; t < Max; t += H) {
		if (Pol(K_Pol, t) * Pol(K_Pol, t + H) < 0) {
			A[i] = t + H / 2;
			i++;
		}
	}
}


int main() {

	for (int i = 0; i < n; i++)
	{
	Yitie();
	}
	ZagatovkaDlyaGaussa();
	Gauss(0, q1);
	for (int i = 0; i < n; i++)
	{
		p1[i] = q1[i][4];
	}

	// C помощью метода гауса нашли все p, они равны соответственно	
	double p[n] = { 6,0.2,-12.735,2.7616 };
	//имея эти значения можем найти Лямбды
	double l[n] = { -1.4201,0.2226,1.5454,5.652 };

	for (int i = 0; i < n; i++)
	{
		q[0][i] = 1;
		for (int j = 1; j < n; j++)
		{	
			q[j][i] = l[i] * q[j-1][i] - p[j-1];
		}
		
	}

	cout << "\n";

		y_0[0] = { 1 };
		y_0[1] = { 0 };
		y_0[2] = { 0 };
		y_0[3] = { 0 };

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				x[i][j] = YItoe[i][2] + q[1][j] * YItoe[i][1] + q[2][j] * YItoe[i][0] + q[3][j] * y_0[i];
			}
				
		}
	
	

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) {
			cout << q[i][j] << " ";
		}
		cout << "\n";
	}
	cout << endl;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) {
			cout << YItoe[i][j] << " ";
		}
		cout << "\n";
	}

	cout << endl;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) {
			cout << x[i][j] << " ";
		}
		cout << "\n";
	}

	cout << endl;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < N; j++) {
			cout << q1[i][j] << " ";
		}
		cout << "\n";
	}

	for (size_t i = 0; i < n; i++) {
		K_Pol[i] = p[i];
	}

	Search(S_korn);
	for (size_t i = 0; i < n; i++) {
		cout << endl << Nyuton(S_korn[i]) << endl;
		cout << p1[i] << endl;
	}

	return 1;
}