// ConsoleApplication16.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
using namespace std;


double f(double x, int param) {
	switch (param) {
	case 1:
		return x * x;
	case 2:
		return 1 / (1 + x * x);
	case 3:
		return 1 / (atan(1 + 10 * x * x));
	}
	return x * x;
}
vector<vector<double>> grid(double a, double b, int n, int param) {
	vector<vector<double>>grid(n, vector<double>(2));
	double x;
	double step = (b - a) / (n - 1);
	for (int i = 0; i < n; i++) {
		x = a + i * step;
		grid[i][0] = x;
		grid[i][1] = f(x, param);
	}
	return grid;
}
vector<vector<double>> gridch(double a, double b, int n, int param) {
	vector<vector<double>>grid(n, vector<double>(2));
	double x;
	double plus = (b + a) * 0.5;
	double minus = (b - a) * 0.5;
	for (int i = 0; i < n; i++) {
		x = plus + minus * cos((2 * (n - 1 - i) + 1) * 3.14159265358 * 0.5 / n);
		grid[i][0] = x;
		grid[i][1] = f(x, param);
	}
	return grid;
}

double ck(int n, vector<vector<double>> grid, int k, double x) { //n degree of polinom
	double p = 1;
	double xk = grid[k][0];
	for (int i = 0; i < n; i++) {
		if (i != k) {
			p *= (x - grid[i][0]) / (xk - grid[i][0]);
		}
	}
	return p;
}
double lagrange(int n, vector<vector<double>>grid, double x) {
	double ans = 0;
	for (int i = 0; i < n; i++) {
		ans += ck(n, grid, i, x) * grid[i][1];
	}
	return ans;
}
vector<double> progonka(int n, vector<vector<double>>A, vector<double> B) {
	vector<double>a(n);
	vector<double>b(n);
	vector<double>res(n);
	double y = A[0][0];
	a[0] = -A[0][1] / y;
	b[0] = B[0] / y;
	for (int i = 1; i < n - 1; i++) {
		y = A[i][i] + A[i][i - 1] * a[i - 1];
		a[i] = -A[i][i + 1] / y;
		b[i] = (B[i] - A[i][i - 1] * b[i - 1]) / y;
	}
	res[n - 1] = (B[n - 1] - A[n - 1][n - 2] * b[n - 2]) / (A[n - 1][n - 1] + A[n - 1][n - 2] * a[n - 2]);
	for (int i = n - 2; i >= 0; i--) {
		res[i] = a[i] * res[i + 1] + b[i];
	}
	return res;
}
int poisk(int n, vector<vector<double>>grid, double x) {
	int l = 0;
	for (int i = 0; i < n - 1; i++) {
		if (x >= grid[i][0] and x <= grid[i + 1][0]) {
			return i;
		}
	}

}
double splain(int n, vector<vector<double>>grid, double x) {
	//vector<vector<double>>grid(n, vector<double>(n));
	vector<vector<double>>t(n, vector<double>(n));
	vector<double>y(n);
	vector<double>a(n);
	vector<double>b(n);
	vector<double>c(n);
	vector<double>d(n);
	vector<double>h(n);
	for (int i = 0; i < n - 1; i++) {
		h[i] = grid[i + 1][0] - grid[i][0];
	}

	y[0] = 0;
	y[n - 1] = 0;
	for (int i = 1; i < n - 1; i++) {
		y[i] = 6 * ((grid[i + 1][1] - grid[i][1]) / h[i] - (grid[i][1] - grid[i - 1][1]) / h[i - 1]);
	}
	t[0][0] = 1;
	t[n - 1][n - 1] = 1;
	for (int i = 1; i < n - 1; i++) {
		t[i][i] = 2 * (h[i - 1] + h[i]);//
		t[i][i - 1] = h[i - 1];//h[0]
		t[i][i + 1] = h[i];//h[1]
	}
	vector<double> m = progonka(n, t, y);

	for (int i = 0; i < n - 1; i++) {
		a[i] = grid[i][1];
		b[i] = (grid[i + 1][1] - grid[i][1]) / h[i] - h[i] * (m[i + 1] - m[i]) / 6 - h[i] * m[i] / 2;
		c[i] = 0.5 * m[i];
		d[i] = (m[i + 1] - m[i]) / (6 * h[i]);
		//cout << "i=" << i <<" abcd:"<< endl;
		//cout << a[i] << " " << b[i] << " " << c[i] << " " << d[i] << endl;
	}
	int ind = poisk(n, grid, x);
	//cout << "ind= " << ind << endl;
	//cout << "gridind= " << grid[ind][0] << " "<< grid[ind+1][0] <<endl;
	double ans = a[ind] + (x - grid[ind][0]) * b[ind] + pow((x - grid[ind][0]), 2) * c[ind] + pow((x - grid[ind][0]), 3) * d[ind];
	return ans;
}
void write(int n, vector<vector<double>>grid) {
	ofstream f2;
	f2.open("C:\\Users\\vfrcb\\Downloads\\2.txt");
	f2 << n << endl;
	for (int i = 0; i < n; i++) {
		f2 << grid[i][0] << " " << grid[i][1] << endl;
	}
	double step = 0.01;
	for (double i = grid[0][0]; i <= grid[n - 1][0]; i += step) {
		f2 << lagrange(n, grid, i) << " ";
	}
	f2.close();
}
int main() {
	int n = 5;
	vector<vector<double>>mygrid = gridch(-1, 1, n, 2);
	//cout<<splain(10, mygrid, -1)<<endl;
	write(n, mygrid);
	/*for (int i = 0; i < 10; i++) {
		cout << mygrid[i][0] << " " << mygrid[i][1] << endl;
	}*/
	//cout<<lagrange(10, mygrid, 1);
	//vector<vector<double>>t = { {,1,0,0},{1,1,1},{0,1,0} };
	//vector<double>b = { 2,2,2 };
	/*b=progonka(3, t, b);
	for (int i = 0; i < 3; i++) {
		cout << b[i]<<endl;
	}*/


}

