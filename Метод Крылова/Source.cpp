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




int main() {

	for (int i = 0; i < n; i++)
	{
	Yitie();
	}
	
	// C помощью метода гауса нашли все p, они равны соответственно	
	double p[n] = { 6,0.2,-12.735,2.7616 };
	//име€ эти значени€ можем найти Ћ€мбды
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
	
	return 1;
	
}