#include <iostream>
using namespace std;
#include <vector>
#include <math.h>
#include <fstream>
double eps = 0.00000001;
double w = 0.5;
double tay = 100;
template<typename t>
vector<t> step(int n, vector<vector<t>>a, vector<t>b, vector<t>& x, int ceo) {
	vector<t>zero(n, 0);
	vector<t>temp = x;
	switch (ceo) {
	case 1:
		for (int i = 0; i < n; i++) {
			x[i] = tay * b[i] / (tay * a[i][i] + 1);
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					x[i] += (tay * a[i][j] * temp[j]) / (-tay* a[i][i] + 1);
				}
			}
		}
		return x;
	case 2:

		for (int i = 0; i < n; i++) {
			x[i] = b[i] / a[i][i];
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j) {
					x[i] = x[i] - temp[j] * (a[i][j] / a[i][i]);
				}
			}

		}
		return x;
	case 3:
		for (int i = 0; i < n; i++) {
			x[i] = b[i] / a[i][i];
		}
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (j < i) {
					x[i] = x[i] - x[j] * a[i][j] / a[i][i];
				}
				else if (j != i) {
					x[i] = x[i] - temp[j] * a[i][j] / a[i][i];
				}
			}
		}
		return x;
	case 4:
		for (int i = 0; i < n; i++) {
			x[i] = w * b[i] / a[i][i];
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (j < i) {
					x[i] = x[i] - w * x[j] * a[i][j] / a[i][i];
				}
				else if (j != i) {
					x[i] = x[i] - w * temp[j] * a[i][j] / a[i][i];
				}
				else {
					x[i] = x[i] + (1 - w) * temp[j];
				}
			}
		}
		return x;
	case 5:
		x[0] = -temp[1] / a[0][0] + 6 / a[0][0];
		for (int i = 1; i < n - 1; i++) {
			x[i] = -(x[i - 1] / a[0][i] + temp[i + 1] / a[0][i]) + (10 - 2 * ((i + 1) % 2)) / a[0][i];
		}
		x[n - 1] = -x[n - 2] / a[0][n - 1] + (9 - 3 * (n % 2)) / a[0][n - 1];
		return x;
	case 6:
		t d = a[0][0];
		x[0] = -w * temp[1] / d + (1 - w) * x[0] + w * 6 / d;
		for (int i = 1; i < n - 1; i++) {
			x[i] = -(w * x[i - 1] / d + w * temp[i + 1] / d) + (1 - w) * temp[i] + w * (10 - 2 * ((i + 1) % 2)) / d;
		}
		x[n - 1] = -w * x[n - 2] / d + (1 - w) * temp[n - 1] + w * (9 - 3 * (n % 2)) / d;
		return x;
	}
	//return x;
}

template<typename t>
t norm0(int n, vector<t> x, vector<t>x2) {
	t ans = 0;
	for (int i = 0; i < n; i++) {
		ans += (x[i] - x2[i]) * (x[i] - x2[i]);
	}
	return sqrt(ans);
}
template<typename t>
t norm00(int n, vector<t> x) {
	t ans = 0;
	for (int i = 0; i < n; i++) {
		ans += (x[i]) * (x[i]);
	}
	return sqrt(ans);
}
template<typename t>
void printLDU(int n, vector <vector<t>> a, t w) {
	cout << "U = " << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i < j)
			{
				cout << -w * a[i][j] / a[i][i] << " ";
			}
			cout << " ";

		}
		cout << endl;
	}
	cout << "D = " << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				cout << 1 - w << " ";
			}
			cout << " ";

		}
		cout << endl;
	}
	cout << "L = " << endl;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i > j)
			{
				cout << -w * a[i][j] / a[i][i] << " ";
			}
			cout << " ";

		}
		cout << endl;
	}
}
template<typename t>
t normrelaxb(int n, vector <vector<t>> a, t w) {
	t max = abs(-w * a[2][0] / a[0][0]) + abs(1 - w);
	if (abs(-w * a[1][n - 2] / a[0][n - 1]) + abs(1 - w) > max) {
		max = abs(-w * a[1][n - 2] / a[0][n - 1]) + abs(1 - w);
	}


	for (int j = 1; j < n - 1; ++j) {
		t sum = abs(-w * a[1][j] / a[0][j]) + abs(-w * a[2][j] / a[0][j]) + abs(1 - w);
		if (sum > max)
			max = sum;
	}
	return max;
}
template<typename t>
t normseidb(int n, vector <vector<t>> a) {
	t max = a[2][0] / a[0][0];
	if (a[1][n - 2] / a[0][n - 1] > max) {
		max = a[1][n - 2] / a[n - 1][n - 1];
	}
	for (int j = 1; j < n - 1; j++) {
		t sum = 0;

		sum = abs(a[1][j] / a[0][j]) + abs(a[2][j] / a[0][j]);

		if (sum > max)
			max = sum;
	}
	return max;
}

template<typename t>
t universalnorm(int n, vector<vector<t>>a, vector<t>b, int ceo, double tay, double w) {
	vector<vector<t>> e;
	vector<t> x(n, 0.0);
	vector<t>z(n, 0);
	//vector<t> b(n, 0.0);
	t sum = 0.0;
	for (int i = 0; i < n; i++) {
		x[i] = 1.0;
		e.push_back(x);
		x[i] = 0.0;
	}
	vector<t>c(n, 0);
	double max = 0;
	vector<t>nor(n, 0);
	switch (ceo) {
	case 1://простой итерации
		for (int i = 0; i < n; i++) {
			x = e[i];
			c = step(n, a, z, x, 1);
			for (int j = 0; j < n; j++) {
				nor[j] += abs(c[j]);
			}
			//if (nor > max) { max = nor; }
		}
		break;


	case 2://якоби
		for (int i = 0; i < n; i++) {
			nor[i] = 0;
		}
		for (int i = 0; i < n; i++) {
			x = e[i];
			x = step(n, a, z, x, 2);
			for (int j = 0; j < n; j++) {
				nor[j] += abs(x[j]);
			}
		}
		break;

	case 3://зейделя
		for (int i = 0; i < n; i++) {
			x = e[i];
			c = step(n, a, z, x, 3);
			for (int j = 0; j < n; j++) {
				nor[j] += abs(c[j]);
			}
		}
		break;

	case 4://релаксация
		if (w == 1.0) { return universalnorm(n, a, b, 3, tay, w); }
		for (int i = 0; i < n; i++) {
			x = e[i];
			c = step(n, a, z, x, 4);
			for (int j = 0; j < n; j++) {
				nor[j] += abs(c[j]);
			}
		}

		break;
	case 5:
		return normseidb(n, a);
		//for (int i = 0; i < n; i++) {
		//	x = e[i];
		//	c = step(n, a, z, x, 5);
		//	//double nor = 0;
		//	for (int j = 0; j < n; j++) {
		//		nor[j] += abs(c[j]);
		//		cout << nor[j] << " ";
		//	}
		//	//if (nor > max) { max = nor; }
		//}
		//break;
	case 6:
		return normrelaxb(n, a, w);
		//for (int i = 0; i < n; i++) {
		//	x = e[i];
		//	c = step(n, a, z, x, 6);
		//	//double nor = 0;
		//	for (int j = 0; j < n; j++) {
		//		nor[j] += abs(c[j]);
		//	}
		//	//if (nor > max) { max = nor; }
		//}
		//break;
	}for (int i = 0; i < n; i++) {
		if (nor[i] > max) { max = nor[i]; }
	}
	return max;
}

template<typename t>
void prost(int n, vector<vector<t>>a, vector<t>b, t tay) {
	vector<t>x(n);
	int k = 0;
	vector<t>temp(n);
	t normc = universalnorm(n, a, b, 1, tay, w);
	do {
		k++;
		temp = x;
		x = step(n, a, b, x, 1);
	} while (norm0(n, x, temp) > eps * (1 - normc) / normc);
	for (int i = 0; i < n; i++) {
		cout << x[i] << " ";
	}
	cout << endl << "norm C = " << normc;
	cout << endl << "iteration = " << k << endl << endl;
}
template<typename t>
void jac(int n, vector<vector<t>>a, vector<t>b) {
	vector<t>x(n);
	int k = 0;
	for (int i = 0; i < n; i++) {
		x[i] = b[i] / a[i][i];
	}
	vector<t>temp(n);
	t normc = universalnorm(n, a, b, 2, tay, w);
	do {
		k++;
		temp = x;
		x = step(n, a, b, x, 2);
	} while (norm0(n, x, temp) > eps * (1 - normc) / normc);
	for (int i = 0; i < n; i++) {
		cout << x[i] << " ";
	}
	cout << endl;
	cout << "norm C = " << normc << endl;
	cout << "iteration = " << k << endl << endl;

}
template<typename t>
void seid(int n, vector<vector<t>>a, vector<t>b) {
	vector<t>x(n);
	int k = 0;

	if (n < 100) {
		vector<t>temp(n);
		t normc = universalnorm(n, a, b, 3, tay, w);
		do {
			k++;
			temp = x;
			x = step(n, a, b, x, 3);
		} while ((norm0(n, x, temp)) > eps);
		for (int i = 0; i < n; i++) {
			cout << x[i] << " ";
		}
		cout << endl << "iteration = " << k << endl;
		cout << "Norm C = " << normc << endl << endl;
	}
	else {
		t normc = universalnorm(n, a, b, 5, tay, w);
		t d = a[0][0];
		x[0] = 6 / a[0][0];
		for (int i = 1; i < n - 1; i++) {
			x[i] = (10 - 2 * ((i + 1) % 2)) / a[0][i];
		}
		x[n - 1] = (9 - 3 * (n % 2)) / a[0][n - 1];
		vector<t>temp(n);
		do {
			k++;
			temp = x;
			x = step(n, a, b, x, 5);
		} while (norm0(n, x, temp) > (1 - normc) / normc * eps);
		for (int i = 0; i < n; i++) {
			cout << x[i] << " ";
		}
		cout << endl << "iteration = " << k << endl;
		cout << "Norm C = " << normc << endl;
	}
}
template<typename t>
void relax(int n, vector<vector<t>>a, vector<t>b, t w) {
	if (w == 1.0) { return seid(n, a, b); }
	vector<t>x(n);
	int k = 0;
	if (n < 100) {
		for (int i = 0; i < n; i++) {
			x[i] = b[i] / a[i][i];
		}
		vector<t>temp(n);
		t normc = universalnorm(n, a, b, 4, tay, w);
		do {
			k++;
			temp = x;
			x = step(n, a, b, x, 4);
		} while (norm0(n, x, temp) > eps * (1 - normc) / normc);
		for (int i = 0; i < n; i++) {
			cout << x[i] << " ";
		}
		cout << endl << "iteration = " << k << endl;
		cout << "Norm C = " << normc << endl << endl;
	}
	else {
		t d = a[0][0];
		t normc = universalnorm(n, a, b, 6, tay, w);
		x[0] = 6 * w / d;
		for (int i = 1; i < n - 1; i++) {
			x[i] = w * (10 - 2 * ((i + 1) % 2)) / a[0][i];
		}
		x[n - 1] = w * (9 - 3 * (n % 2)) / a[0][n - 1];
		vector<t>temp(n);
		do {
			k++;
			temp = x;
			x = step(n, a, b, x, 6);
		} while ((norm0(n, x, temp)) > (1 - normc) / normc * eps);
		for (int i = 0; i < n; i++) {
			cout << x[i] << " ";
		}
		cout << endl << "iteration = " << k << endl;
		cout << "Norm C = " << normc << endl;
	}
	//printLDU(n, a, w);
}

void read_double(int& n, vector<vector<double>>& a, vector<double>& b, bool flag) {
	ifstream f1;
	f1.open("1.txt");
	if (flag == 1)
	{
		f1.seekg(102, ios_base::beg);
	}

	f1 >> n;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			f1 >> a[i][j];
		}
	}
	for (int j = 0; j < n; j++) {
		f1 >> b[j];
	}


	f1.close();
}
void write(int n, vector<double>x) {
	ofstream f2;
	f2.open("2.txt");
	for (int i = 0; i < n; i++) {
		f2 << x[i];
	}
	f2.close();

}
int main() {
	setlocale(LC_ALL, "russian");

	int N = 200;
	int n;
	double tay = 0.001;
	double w = 1.0;
	vector<vector<double>>test1 = { {15,2,-3,7}, {-5,11,2,-3},{0,-1,7,4},{12,0,-6,20} };
	vector<double>b1 = { 53,-90,107,68 };
	vector<vector<double>>test2 = { {86.00,-8.93,-9.59, -3.91}, {4.05,-100.00,-9.10,-8.14},{0.26,3.61,-71.80,-4.28},{-4.03,-6.88,6.57,-198.60} };
	vector<double>b2 = { 818.58,898.74,-912.22,-687.06 };
	//vector<vector<double>>test3 = { {4} };																																																							
																																																																																																					//tay=tay*100000;
	vector<double>left(N + 7, 1.0);
	vector<double>right(N + 7, 1.0);
	vector<double>diag(N + 8, 4.0);
	vector<vector<double>>test3 = { diag,left,right };
	vector<double>x(4);
	bool flag = 0;
	int test;
	vector<double>e = { 1,1,1,1 };
	vector<vector<double>>a = { e,e,e,e };
	vector<double>b(4);

	/*cout << "Введите номер теста (1 или 2) = ";
	cin >> test;
	cout << endl;
	if (test == 2)
	{
		flag = 1;
	}
	read_double(n, a, b, flag);*/
	//prost(n, a, b, tay);
	//prost(4, test2, b2, tay);
	//prost(4, test2, b2, tay);
	/*prost(4, test1, b1, tay);
	jac(4, test1, b1);
	seid(4, test1, b1);
	relax(4, test1, b1, w);*/

	prost(4, test2, b2, tay);
	/*jac(4, test2, b2);
	seid(4, test2, b2);
	relax(4, test2, b2, w);*/
	cout << endl << universalnorm(4, test1, b1, 1, tay, w);
	/*cout<<endl<<universalnorm(4, test1,b1, 2, tay,w);
	cout << endl << universalnorm(4, test1, b1, 3, tay, w);
	cout << endl << universalnorm(4, test1, b1, 4, tay, w);
	cout << endl << universalnorm(4, test2, b2, 1, tay, w);
	cout << endl << universalnorm(4, test2, b2, 2, tay, w);
	cout << endl << universalnorm(4, test2, b2, 3, tay, w);
	cout << endl << universalnorm(4, test2, b2, 4, tay, w) << endl;*/
	//seid(208, test3, b2);
	//relax(208, test3, b2, w);
	//cout << endl << universalnorm(208, test3, b2, 5, tay, w);
	//cout << endl << universalnorm(208, test3, b2, 6, tay, w) << endl;

}