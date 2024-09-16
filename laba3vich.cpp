
#include <iostream>
using namespace std;
#include <vector>
#include <math.h>
#include <fstream>
double eps = 0.00001;

template<typename t>
vector<vector<t>> one(int n, vector<vector<t>>& c) {
	//vector<vector<t>>c(n, vector<t>(n));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) { c[i][j] = 1; }
			else { c[i][j] = 0; }
		}
	}
	return c;
}
double null_round(double num) {
	if (abs(num) < eps) {
		return 0;
	}
	return num;
}
template<typename t>
void print(int n, vector <vector<t>>& a) {
	vector<t>x(n);
	for (int i = 0; i < n; i++) {
		cout << endl;
		for (int j = 0; j < n; j++) {
			if (abs(a[i][j]) < 0.0001) { a[i][j] = 0; }
			cout << a[i][j] << " ";
		}
	}
	cout << endl;
}
template<typename t>
void trans(int n, vector <vector<t>>& a) {
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			t temp = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = temp;
		}
	}
}

template<typename t>
vector<t> inverse_step(int n, vector <vector<t>> a, vector<t>b) {
	vector<t>x(n);
	for (int k = n - 1; k >= 0; k--) {
		x[k] = b[k] / a[k][k];
		for (int j = 0; j < k; j++) {
			b[j] = b[j] - a[j][k] * x[k];
		}
	}
	return x;
}

template<typename t>
t norm0(int n, vector<t> x) {
	t ans = 0;
	for (int i = 0; i < n; i++) {
		ans += x[i] * x[i];
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
template<typename t>

//Октаэдрическая норма
t norm1(int n, vector <vector<t>> a) {
	t max = 0.0;
	for (int j = 0; j < n; ++j) {
		t sum = 0;
		for (size_t i = 0; i < n; ++i)
			sum += abs(a[i][j]);
		if (sum > max)
			max = sum;
	}
	return max;
}

template<typename t>

t norm2(int n, vector <vector<t>> a) {
	t max = 0.0;
	for (int j = 0; j < n; ++j) {
		t sum = 0;
		for (int i = 0; i < n; ++i)
			sum += abs(a[j][i]);
		if (sum > max)
			max = sum;
	}
	return max;
}
template<typename t>
vector<t> qr2(int n, vector <vector<t>> a, vector <vector<t>>& r, vector <vector<t>>& q, bool
	param) {
	vector<vector<t>>c(n, vector<t>(n));
	vector<vector<t>>s(n, vector<double>(n));
	vector<vector<t>>e(n, vector<double>(n));

	if (param == 1) { print(n, a); }
	one(n, e);
	vector<vector<t>>a1 = a;
	vector<vector<t>>orig(n, vector<double>(n));
	vector<vector<t>>orig2(n, vector<double>(n));
	//vector<t>orig3(n);
	for (int k = 0; k < n - 1; k++) {
		for (int l = k + 1; l < n; l++) {
			if (abs(a[l][k]) < eps) { continue; }
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					orig[i][j] = a[i][j]; orig2[i][j] = e[i][j]; //orig3[j] = b[j];
				}
			}

			t co; t so;
			if (max(abs(a[k][k]), abs(a[l][k])) < eps) { co = 1;  so = 0; }
			else {
				t temp = sqrt(a[k][k] * a[k][k] + a[l][k] * a[l][k]);
				co = a[k][k] / temp;
				so = a[l][k] / temp;
			}
			t ckk = co;
			t ckl = so;
			t clk = -so;
			t cll = co;

			for (int j = 0; j < n; j++) {
				e[k][j] = 0; a[k][j] = 0; e[l][j] = 0; a[l][j] = 0;
				a[k][j] += ckk * orig[k][j]; e[k][j] += ckk * orig2[k][j];
				a[l][j] += clk * orig[k][j]; e[l][j] += clk * orig2[k][j];
				a[k][j] += ckl * orig[l][j]; e[k][j] += ckl * orig2[l][j];
				a[l][j] += cll * orig[l][j]; e[l][j] += cll * orig2[l][j];
			}
			/*for (int j = 0; j < n; j++) {b[j] = 0;for (int z = 0; z < n; z++) { b[j] += c[z][j]
			* orig3[j]; }}*/
		}
	}
	r = a;
	vector<vector<t>>q = e;
	trans(n, q);
}

void sign_vec(vector<double>& vec, int n) {
	if (vec[0] < 0) {
		for (int i = 0; i < n; i++) {
			vec[i] *= -1;
		}
	}
}

template<typename t>
vector<t> gaus(int n, int m, vector <vector<t>> a, vector<t>b, bool param) {
	if (param) { print(n, a); }
	vector<t>b1 = b;
	vector<vector<t>>a1 = a;
	int i = 0;
	while (i < m) {
		//ищем максимальный элемент в i столбце
		t max = abs(a[i][i]);
		int ind = i;
		for (int j = i; j < m; j++) {
			if (abs(a[j][i]) > max) {
				max = abs(a[j][i]);
				ind = j;
			}
		}
		if ((max < eps) and (param == 1)) { cout << "zero colom, not one solution" << endl; }

		//при необходимости меняем строки местами
		if (ind != i) {
			vector<t>temp(m);
			for (int k = 0; k < m; k++) {
				temp[k] = a[i][k];
				a[i][k] = a[ind][k];
			}
			for (int k = 0; k < m; k++) {
				a[ind][k] = temp[k];
			}
			t temp2 = b[i];
			b[i] = b[ind];
			b[ind] = temp2;
		}
		for (int k = i; k < m; k++) {
			t temp = a[k][i];
			if (abs(temp) < eps) continue;
			for (int j = 0; j < n; j++) {
				a[k][j] = a[k][j] / temp;
			}
			b[k] = b[k] / temp;
			if (k == i) { continue; }
			for (int j = 0; j < n; j++) {
				a[k][j] = a[k][j] - a[i][j];
			}
			b[k] = b[k] - b[i];
		}
		i++;
	}
	//обратная подстановка
	vector<t>x = inverse_step(n, a, b);
	if (param == 0) { return x; }
	cout << "Result: " << endl;
	for (int i = 0; i < n; i++) {
		cout << x[i] << endl;
	}
	//ищем вектор невязки
	//if (param) { disjoint(n, a1, b1, x); }

	return x;
}

void read_double(int& n, vector<vector<double>>& a, vector<double>& b, int k) {
	ifstream f1;
	f1.open("C:\\Users\\vfrcb\\Downloads\\1.txt");
	for (int u = 0; u <= k; u++) {
		f1 >> n;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				f1 >> a[i][j];
			}
		}
	}
	f1.close();
}
void write(int n, vector<double>x) {
	ofstream f2;
	f2.open("C:\\Users\\vfrcb\\Downloads\\2.txt");
	for (int i = 0; i < n; i++) {
		f2 << x[i] << endl;
	}
	f2.close();

}

vector<double> mulvm(int n, vector<vector<double>> a, vector<double>x) {
	vector<double>temp(n, 0);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			temp[i] += a[i][j] * x[j];
		}
	}
	return temp;
}
template<typename t>
void hessenberg(vector<vector<t>>& matrix, int N) {
	double c = 0, s = 0, sq = 0;
	for (int k = 1; k < N - 1; k++) {
		for (int l = k + 1; l < N; l++) {
			sq = null_round(sqrt(matrix[k][k - 1] * matrix[k][k - 1] + matrix[l][k - 1] *
				matrix[l][k - 1]));
			if (sq == 0) {
				c = 1, s = 0;
			}
			else {
				c = matrix[k][k - 1] / sq;
				s = matrix[l][k - 1] / sq;
			}
			double temp;
			for (int j = 0; j < N; j++) {
				temp = matrix[k][j];
				matrix[k][j] = (matrix[k][j] * c + matrix[l][j] * s);
				matrix[l][j] = (matrix[l][j] * c - temp * s);
			}
			for (int j = k - 1; j < N; j++) {
				temp = matrix[j][k];
				matrix[j][k] = (matrix[j][k] * c + matrix[j][l] * s);
				matrix[j][l] = (matrix[j][l] * c - temp * s);
			}
		}
	}
}

template<typename t>
void shift(int n, vector <vector<t>>& a, double sigm) {
	for (int i = 0; i < n; i++) {
		a[i][i] -= sigm;
	}
}
template<typename t>
void qr(int n, vector <vector<t>>& a, vector <vector<t>>& r, vector <vector<t>>& q) {
	vector<vector<t>>e(n, vector<double>(n));
	one(n, e);
	vector<vector<t>>orig(n, vector<double>(n));
	vector<vector<t>>orig2(n, vector<double>(n));
	for (int k = 0; k < n - 1; k++) {
		if (abs(a[k + 1][k]) < 0.0000001) { continue; }
		int l = k + 1;
		orig = a; orig2 = e;
		t co; t so;
		if (max(abs(a[k][k]), abs(a[l][k])) < eps) { co = 1;  so = 0; }
		else {
			t temp = null_round(sqrt(a[k][k] * a[k][k] + a[l][k] * a[l][k]));
			co = a[k][k] / temp;
			so = a[k + 1][k] / temp;
		}
		for (int j = 0; j < n; j++) {
			e[k][j] = 0; a[k][j] = 0; e[k + 1][j] = 0; a[k + 1][j] = 0;
			a[k][j] += co * orig[k][j]; e[k][j] += co * orig2[k][j];
			a[l][j] += (-so) * orig[k][j]; e[l][j] += (-so) * orig2[k][j];
			a[k][j] += so * orig[l][j]; e[k][j] += so * orig2[l][j];
			a[l][j] += co * orig[l][j]; e[l][j] += co * orig2[l][j];
		}
	}
	r = a;
	trans(n, e);
	q = e;
}

template<typename t>
void mult(int n, vector <vector<t>>& r, vector <vector<t>>& q, vector <vector<t>>& pr) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			pr[i][j] = 0;
			for (int z = 0; z < n; z++) {
				pr[i][j] += r[i][z] * q[z][j];
			}
		}
	}
}
template<typename t>
bool qr_stop_case(vector<vector<t>> a, int N) {
	bool flag = 0;
	for (int i = 0; i < N - 1; i++) {
		for (int j = 0; j < N - 1; j++) {
			if ((i > j) && a[i][j] > eps) { flag = 1; }
		}
	}
	return flag;
}
template<typename t>
void find(int n, vector <vector<t>>& test) {
	int k = 0;
	vector<vector<double>>r = test;
	vector<vector<double>>q;
	double s;
	while (qr_stop_case(test, n)) {
		k++;
		s = test[3][3];
		shift(4, test, s);
		qr(4, test, r, q);
		mult(4, r, q, test);
		shift(4, test, -s);
	}
	cout << k << endl;
}
template<typename t>
void inviter(int n, vector <vector<t>> a, double lambda) {
	vector<t>x(n, 0);
	vector<t>xa(n, 0);
	vector<t>xprev(n, 0);
	vector<t>y(n, 0);
	int k = 0;
	double e = 0.0000001;
	double norm = 0;
	for (int i = 0; i < n; i++) {
		x[i] = 0;
		y[i] = 0;
		xprev[i] = 0;
	}
	x[0] = 1.0;
	for (int i = 0; i < n; i++) {//A-lambda_E
		a[i][i] -= lambda;
	}
	while (norm01(n, x, xprev) > e * norm0(n, xprev) and k < 10) {
		xprev = x;
		y = gaus(n, n, a, x, 0);
		x = y;
		norm = norm0(n, x);
		for (int j = 0; j < n; j++) {
			x[j] /= norm;
		}
		xa = x;
		sign_vec(x, n);
		k++;
	}
	cout << endl << k << endl;
	for (int i = 0; i < n; i++) {
		cout << xa[i] << " ";
	}
}
template<typename t>
void rel(vector <vector<t>> a, int n) {
	vector<t>xprev(n, 0);
	vector<t>y(n, 0);
	vector<t>x(n, 0);
	vector<t>xmi(n, 0);
	vector<t>xma(n, 0);
	double lambda = 0;
	int k = 0;
	double minim = 100000;
	double maxim = -10000;
	double e = 0.0001;
	double norm = 0;
	for (int i = 0; i < n; i++) {
		y[i] = 0;
		x[i] = 1;
		xprev[i] = 0;
	}
	x[0] = 1.0;

	/*x[0] = 0.50;
	x[1] = 0.00;
	x[2] = -0.43;
	x[3] = -0.75;*/
	norm = norm0(n, x);
	for (int j = 0; j < n; j++) {
		x[j] /= norm;
	}
	while (norm01(n, x, xprev) > eps or k < 100) {
		xprev = x;
		vector<t> Ax = mulvm(n, a, x);
		lambda = 0;
		for (int i = 0; i < n; i++) {
			lambda += Ax[i] * x[i];
		}
		if (lambda < minim) { minim = lambda; xmi = gaus(n, n, a, x, 0); }
		if (lambda > maxim) { maxim = lambda; xma = gaus(n, n, a, x, 0); }

		for (int i = 0; i < n; i++) {
			a[i][i] -= lambda;
		}
		y = gaus(n, n, a, x, 0);
		/*if (norm0(n,y)<eps) {
		break;
		}*/
		x = y;
		norm = norm0(n, x);
		for (int j = 0; j < n; j++) {
			x[j] /= norm;
		}
		for (int i = 0; i < n; i++) {
			a[i][i] += lambda;
		}
		k++;
	}
	cout << "lambda: " << lambda << endl;
	for (int i = 0; i < n; i++) {
		cout << -x[i] << " ";
	}
	/*cout << "minimum: " << minim << endl;
	cout << endl;
	for (int i = 0; i < n; i++) {
	cout << xmi[i] << " ";
	}
	cout <<endl<< "maximum " << maxim << endl;
	for (int i = 0; i < n; i++) {
	cout << xma[i] << " ";
	}*/
}
int main() {
	vector<vector<double>>test = { {1.50,0.0,-0.43,-0.75},{0.0,3.0,0.87,-0.50},{-
	0.43,0.87,2.90,-0.22},{-0.75,-0.50,-0.22,2.60} };
	print(4, test);
	vector<vector<double>>a = test;
	hessenberg(test, 4);
	print(4, test);
	find(4, test);
	print(4, test);

	for (int i = 0; i < 4; i++) {
		inviter(4, a, test[i][i]);
	}


	cout << endl << "Reley:" << endl;
	rel(a, 4);

}// 