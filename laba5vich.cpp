// ConsoleApplication17.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
#include <iostream>
using namespace std;
#include <vector>
#include<cmath>
#include <fstream>
#include<string>
double eps = 0.00000001;


double norm0(int n, vector<double> x) {
	double ans = 0;
	for (int i = 0; i < n; i++) {
		ans += (x[i] * x[i]);
	}
	return sqrt(ans);
}

template<typename t>
t norm01(int n, vector<t> x, vector<t> x1) {
	t ans = 0;
	for (int i = 0; i < n; i++) {
		ans += ((x[i] - x1[i]) * (x[i] - x1[i]));
	}
	return sqrt(ans);
}
//template<typename t>
//vector<t> step(int n, vector<vector<t>>a, vector<t>b, vector<t>& x, int ceo) {
//	vector<t>zero(n, 0);
//	vector<t>temp = x;
//	switch (ceo) {
//	case 1:
//		for (int i = 0; i < n; i++) {
//			x[i] = 1 / tay * b[i] / (1 / tay * a[i][i] + 1);
//		}
//		for (int i = 0; i < n; i++) {
//			for (int j = 0; j < n; j++) {
//				if (i != j) {
//					x[i] += (1 / tay * a[i][j] * temp[j]) / (-1 / tay * a[i][i] + 1);
//				}
//			}
//		}
//		return x;
//	case 2:
//		for (int i = 0; i < n; i++) {
//			x[i] = b[i] / a[i][i];
//		}
//		for (int i = 0; i < n; i++) {
//			for (int j = 0; j < n; j++) {
//				if (i != j) {
//					x[i] = x[i] - temp[j] * (a[i][j] / a[i][i]);
//				}
//			}
//		}
//		return x;
//	case 3:
//		for (int i = 0; i < n; i++) {
//			x[i] = b[i] / a[i][i];
//		}
//		for (int i = 0; i < n; i++) {
//			for (int j = 0; j < n; j++) {
//				if (j < i) {
//					x[i] = x[i] - x[j] * a[i][j] / a[i][i];
//				}
//				else if (j != i) {
//					x[i] = x[i] - temp[j] * a[i][j] / a[i][i];
//				}
//			}
//		}
//		return x;
//	case 4:
//		for (int i = 0; i < n; i++) {
//			x[i] = w * b[i] / a[i][i];
//		}
//
//		for (int i = 0; i < n; i++) {
//			for (int j = 0; j < n; j++) {
//				if (j < i) {
//					x[i] = x[i] - w * x[j] * a[i][j] / a[i][i];
//				}
//				else if (j != i) {
//					x[i] = x[i] - w * temp[j] * a[i][j] / a[i][i];
//				}
//				else {
//					x[i] = x[i] + (1 - w) * temp[j];
//				}
//			}
//		}
//		return x;
//	}
//}
//template<typename t>
//t universalnorm(int n, vector<vector<t>>a, vector<t>b, int ceo, double tay, double w) {
//	vector<vector<t>> e;
//	vector<t> x(n, 0.0);
//	vector<t>z(n, 0);
//	//vector<t> b(n, 0.0);
//	t sum = 0.0;
//	for (int i = 0; i < n; i++) {
//		x[i] = 1.0;
//		e.push_back(x);
//		x[i] = 0.0;
//	}
//	vector<t>c(n, 0);
//	double max = 0;
//	vector<t>nor(n, 0);
//	switch (ceo) {
//	case 1://простой итерации
//		for (int i = 0; i < n; i++) {
//			x = e[i];
//			c = step(n, a, z, x, 1);
//			for (int j = 0; j < n; j++) {
//				nor[j] += abs(c[j]);
//			}
//			//if (nor > max) { max = nor; }
//		}
//		break;
//	case 2://якоби
//		for (int i = 0; i < n; i++) {
//			nor[i] = 0;
//		}
//		for (int i = 0; i < n; i++) {
//			x = e[i];
//			x = step(n, a, z, x, 2);
//			for (int j = 0; j < n; j++) {
//				nor[j] += abs(x[j]);
//			}
//		}
//		break;
//	case 3://зейделя
//		for (int i = 0; i < n; i++) {
//			x = e[i];
//			c = step(n, a, z, x, 3);
//			for (int j = 0; j < n; j++) {
//				nor[j] += abs(c[j]);
//			}
//		}
//		break;
//	case 4://релаксация
//		if (w == 1.0) { return universalnorm(n, a, b, 3, tay, w); }
//		for (int i = 0; i < n; i++) {
//			x = e[i];
//			c = step(n, a, z, x, 4);
//			for (int j = 0; j < n; j++) {
//				nor[j] += abs(c[j]);
//			}
//		}
//		break;
//	}
//	for (int i = 0; i < n; i++) {
//		if (nor[i] > max) { max = nor[i]; }
//	}
//	return max;
//}
//template<typename t>
//vector<t>seid(int n, vector<vector<t>>a, vector<t>b) {
//	vector<t>x(n);
//	int k = 0;
//	vector<t>temp(n);
//	t normc = universalnorm(n, a, b, 3, tay, w);
//	do {
//		k++;
//		temp = x;
//		x = step(n, a, b, x, 3);
//	} while ((norm0(n, x, temp)) > eps);
//	/*for (int i = 0; i < n; i++) {
//		cout << x[i] << " ";
//	}*/
//	return x;
//	//cout << endl << "iteration = " << k << endl;
//	//cout << "Norm C = " << normc << endl << endl;
//}

double f(double x, int param) {
	switch (param) {
	case 1:
		return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
	case 2:
		return sqrt(x + 1) - 1;
	case 3:
		return 35 * x * x * x - 67 * x * x - 3 * x + 3;
	}
	return x * x;
}
double df(double x, int param) {
	switch (param) {
	case 1:
		return (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75) + (x - 0.1) * (x - 0.55) * (x - 0.7) * (x - 0.75) + (x - 0.1) * (x - 0.22) * (x - 0.7) * (x - 0.75) + (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.75) + (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7);
	case 2:
		return 1 / (2 * sqrt(x + 1));
	case 3:
		return 3 * 35 * x * x - 67 * 2 * x - 3;
	}
	return x * x;
}
double dfch(double x, int param) {
	double h = 0.001;
	return (f((x + h), param) - f(x, param)) / h;
}
vector<double> f2(double x1, double x2, int param) {
	vector<double>x(2);
	switch (param) {
	case 4:
		x[0] = x1 * x1 - x2 * x2 - 15;
		x[1] = x1 * x2 + 4;
		return x;
	case 5:
		x[0] = x1 * x1 + x2 * x2 + x1 + x2 - 8;
		x[1] = x1 * x1 + x2 * x2 + x1 * x2 - 7;
		return x;
	}
	return x;
}
vector<vector<double>> df2(double x1, double x2, int param) {
	vector<vector<double>>x(2, vector<double>(2));
	switch (param) {
	case 4:
		x[0][0] = 2 * x1;
		x[0][1] = -2 * x2;
		x[1][0] = x2;
		x[1][1] = x1;
		return x;
	case 5:
		x[0][0] = 2 * x1 + 1;
		x[0][1] = 2 * x2 + 1;
		x[1][0] = 2 * x1 + x2;
		x[1][1] = 2 * x2 + x1;
		return x;
	}
	return x;
}
vector<vector<double>> df2ch(double x1, double x2, int param) {
	vector<vector<double>>x(2, vector<double>(2));
	double h = 0.0001;
	x[0][0] = (f2(x1 + h, x2, param)[0] - f2(x1 - h, x2, param)[0]) / (2 * h);
	x[0][1] = (f2(x1, x2 + h, param)[0] - f2(x1, x2 - h, param)[0]) / (2 * h);
	x[1][0] = (f2(x1 + h, x2, param)[1] - f2(x1 - h, x2, param)[1]) / (2 * h);
	x[1][1] = (f2(x1, x2 + h, param)[1] - f2(x1, x2 - h, param)[1]) / (2 * h);
	return x;
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

double bis(double a, double b, int param) {
	int k = 0;
	double x;
	while (abs(b - a) > eps) {
		k++;
		x = (a + b) / 2;
		if (f(x, param) * f(a, param) < 0) {
			b = x;
		}
		else {
			a = x;
		}
	}if (abs(f(x, param)) < 0.001) {
		if (abs(x - 0.75) < 0.03) { x = 0.75; }
		cout << "Root: " << x << endl;
		cout << k << endl;
	}
	return (a + b) / 2;
}
double newton(double a, double b, int param) {
	//cout << a << " " << b << endl;
	int k = 0;
	double xk = 0.0;//(a + b) / 2
	//double xk = (f(a, param) * b - f(b, param) * a) / (f(a, param) - f(b, param));

	double x = xk - f(xk, param) / df(xk, param);
	while (abs(xk - x) > eps) {
		k++;
		xk = x;
		x = xk - f(xk, param) / df(xk, param);
		/*if (x > b) { x = (f(b, param) * xk - f(xk, param) * b) / (f(b, param) - f(xk, param)); }
		if (x < a) { x = (f(a, param) * xk - f(xk, param) * a) / (f(a, param) - f(xk, param)); }*/
	}
	if (abs(f(x, param)) < eps) {
		cout << "Root: " << x << endl;
		cout << k << endl;
	}
	return x;
}
vector<double> newton2(vector<double> a, vector<double> b, int param) {
	int k = 0;
	vector<double> xk(2);
	vector<double> x(2);
	x[0] = (a[0] + b[0]) / 2 + eps;
	x[1] = (a[1] + b[1]) / 2 + eps;
	//cout << "intial " << x[0] << " " << x[1] << endl;
	//double xk[0] = (f2(a[0],a[1], param) * b - f2(b[0],b[1], param) * a) / (f2(a[0],a[1], param)[0] - f2(b[0],b[1], param)[0]);
	//x[0] = xk[0] - f(xk, param) / df(xk, param);
	do {
		k++;
		xk = x;
		vector<double>f = f2(xk[0], xk[1], param);
		vector<vector<double>>jac = df2(xk[0], xk[1], param);
		double d = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
		//if (d < eps) { cout << "zero division" << endl; }
		x[0] = xk[0] - (jac[1][1] * f[0] - jac[0][1] * f[1]) / d;
		x[1] = xk[1] - (-jac[1][0] * f[0] + jac[0][0] * f[1]) / d;
		if (x[0] > b[0]) { x[0] = b[0]; }
		if (x[1] > b[1]) { x[1] = b[1]; }
		if (x[0] < a[0]) { x[0] = a[0]; }
		if (x[1] < a[1]) { x[1] = a[1]; }
		if (k > 40) {
			break;
		}
	} while (norm01(2, xk, x) > eps * norm0(2, xk));
	if (norm0(2, f2(x[0], x[1], param)) < 1) {
		cout << " Root  " << x[0] << " ; " << x[1] << endl;
		cout << k << endl;
	}
	return x;
}
void split(int n, vector<vector<double>>grid, int param, int m) {
	double x, fx;
	//int n = 26;
	double xp = grid[0][0];
	double fp = grid[0][1];
	for (int i = 1; i < n; i++) {
		fx = f(grid[i][0], param);
		if (fp * fx <= 0) {
			if (m == 1) {
				newton(xp, grid[i][0], param);
			}
			else { bis(xp, grid[i][0], param); }
			fp = fx;
			xp = grid[i][0];
		}
	}
}
void grid2(double xR, double yR, int N, int M, int param) {
	double hx = 2 * xR / (double)N, hy = hx;//2 * yR / (double)M
	vector<double> xNext(2, 0), xPrev(2, 0);
	for (size_t i = 0; i < N + 1; i++) {
		for (size_t j = 0; j < M + 1; j++) {
			if ((i + 1) <= N and (j + 1) <= M) {
				xPrev[0] = -yR + hy * j;
				xPrev[1] = -xR + hx * i;
				xNext[0] = -yR + hy * (j + 1);
				xNext[1] = -xR + hx * (i + 1);
				//cout << "a= " << xPrev[0] << " " << xPrev[1] << endl;
				//cout << "b= " << xNext[0] << " " << xNext[1] << endl;
				newton2(xPrev, xNext, param);
			}
		}
	}
}
//void write(double a, double b, int param, int nt) {
//	ofstream f2;
//	f2.open("C:\\Users\\vfrcb\\Downloads\\2.txt");
//	if (nt) {
//		double xk = (f(a, param) * b - f(b, param) * a) / (f(a, param) - f(b, param));
//		//ans.push_back(xk);
//		f2 << xk << " ";
//		double x = xk - f(xk, param) / df(xk, param);
//		//ans.push_back(x);
//		f2 << x << " ";
//		int k = 0;
//		while (abs(xk - x) > eps) {
//			xk = x;
//			x = xk - f(xk, param) / df(xk, param);
//			if (x > b) { x = b; }
//			if (x < a) { x = a; }
//			//ans.push_back(x);
//			f2 << x << " ";
//		}
//
//	}
//	else {
//		double x;
//		while (abs(b - a) > eps) {
//			x = (a + b) / 2;
//			//ans.push_back(x);
//			f2 << x << " ";
//			if (f(x, param) * f(a, param) < 0) { b = x; }
//			else { a = x; }
//		}
//	}
//	f2.close();
//}
//void write2(vector<double> a, vector<double> b, int param) {
//	ofstream file;
//	file.open("C:\\Users\\vfrcb\\Downloads\\2.txt");
//	int k = 0;
//	vector<double> xk(2); vector<double> x(2);
//	x[0] = (a[0] + b[0]) / 2 + eps; x[1] = (a[1] + b[1]) / 2 + eps;
//	file << x[0] << " ";
//	file << x[1] << " ";
//	do {
//		k++; xk = x;
//		vector<double>f = f2(xk[0], xk[1], param);
//		vector<vector<double>>jac = df2(xk[0], xk[1], param);
//		double d = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
//		x[0] = xk[0] - (jac[1][1] * f[0] - jac[0][1] * f[1]) / d;
//		x[1] = xk[1] - (-jac[1][0] * f[0] + jac[0][0] * f[1]) / d;
//		if (x[0] > b[0]) { x[0] = b[0]; }
//		if (x[1] > b[1]) { x[1] = b[1]; }
//		if (x[0] < a[0]) { x[0] = a[0]; }
//		if (x[1] < a[1]) { x[1] = a[1]; }
//		if (k > 40) { break; }
//		file << x[0] << " ";
//		file << x[1] << " ";
//	} while (norm01(2, xk, x) > eps * norm0(2, xk));
//	if (norm0(2, f2(x[0], x[1], param)) < 1) {//cout << " Root  " << x[0] << " ; " << x[1] << endl;
//	}
//	file.close();
//}
int main() {
	int n = 20;
	int param = 1;
	vector<vector<double>>mygrid;

	//split(n, mygrid, param,1);
	int i = 1;
	cout << endl << " Test " << i << endl;
	cout << "Newton:" << endl;
	mygrid = grid(0, 1, n, i);
	split(n, mygrid, i, 1);
	cout << "Biss:" << endl;
	split(n, mygrid, i, 2);

	i = 2;
	cout << endl << " Test " << i << endl;
	cout << "Newton:" << endl;
	mygrid = grid(-1, 10, n, i);
	split(n, mygrid, i, 1);
	cout << "Biss:" << endl;
	split(n, mygrid, i, 2);

	i = 3;
	cout << endl << " Test " << i << endl;
	cout << "Newton:" << endl;
	mygrid = grid(0, 1, n, i);
	split(n, mygrid, i, 1);
	cout << "Biss:" << endl;
	split(n, mygrid, i, 2);

	cout << " Test 4" << endl << endl;
	grid2(10, 10, 2, 2, 4);
	cout << " Test 5" << endl << endl;
	grid2(10, 10, 7, 7, 5);
	cout << endl << endl;

	//newton(0, 1, 1);
	/*vector<double>a = { -5,-5 }; vector<double>b = { 0,0 };
	newton2(a, b, 4);
	a = { 0,0 }; b = { 5,5 };
	newton2(a, b, 4);
	a = { -5,0 }; b = { 0,5 };
	newton2(a, b, 4);
	a = { 0,-5 }; b = { 5,0 };
	newton2(a, b, 4);*/

	//write(0, 1, 1, 0);
	//write2({ -10,-10 }, { 10,10 }, 5);
}
