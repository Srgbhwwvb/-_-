// metodvich.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
#include <iostream>
using namespace std;
#include <vector>
#include <math.h>
#include <fstream>
double eps = 0.0000001;

template<typename t>
void print(int n, vector <vector<t>> a) {
	vector<t>x(n);
	for (int i = 0; i < n; i++) {
		cout << endl;
		for (int j = 0; j < n; j++) {
			cout << a[i][j] << " ";
		}
	}
	cout << endl;
}
template<typename t>
void trans(int n, vector <vector<t>> a) {
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			t temp = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = temp;
		}
	}
}
template<typename t>
void mult(int n, vector <vector<t>> a, vector <vector<t>> a1, vector <vector<t>>& pr) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			pr[i][j] = 0;
			for (int z = 0; z < n; z++) {
				pr[i][j] += a1[z][i] * a[z][j];
			}
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
template<typename t>
t norm0(int n, vector<t> x) {
	t ans = 0;
	for (int i = 0; i < n; i++) {
		ans += x[i] * x[i];
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
void disjoint(int n, vector <vector<t>> a1, vector<t>b1, vector<t>x) {
	vector<t>b2(n, 0.0);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			b2[i] += a1[i][j] * x[j];
		}
	}
	vector<t>dif(n);
	t acq = 0.0;
	for (int k = 0; k < n; k++) {
		dif[k] = abs(b1[k] - b2[k]);
		acq += (b1[k] - b2[k]) * (b1[k] - b2[k]);
	}
	cout << "the norm of the disjoint vector is " << sqrt(acq) << endl;

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
	if (param) { disjoint(n, a1, b1, x); }

	return x;
}

template<typename t>
vector<t> qr(int n, int m, vector <vector<t>> a, vector<t>b, bool param) {
	vector<vector<t>>c(n, vector<t>(n));
	vector<vector<t>>s(n, vector<double>(n));
	vector<vector<t>>e(n, vector<double>(n));
	vector<t>b1 = b;
	//
	if (param == 1) { print(n, a); }
	one(n, e);
	vector<vector<t>>a1 = a;

	vector<vector<t>>orig(n, vector<double>(n));
	vector<vector<t>>orig2(n, vector<double>(n));
	for (int k = 0; k < n - 1; k++) {
		for (int l = k + 1; l < n; l++) {
			if (abs(for (int z = 0; z < n; z++) {
					a[k][j] += c[k][z] * orig[z][j];
					e[k][j] += c[k][z] * orig2[z][j];
					a[l][j] += c[l][z] * orig[z][j];
					e[l][j] += c[l][z] * orig2[z][j];
				}a[l][k]) < eps) { continue; }
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					orig[i][j] = a[i][j]; orig2[i][j] = e[i][j];
				}
			}
			one(n, c);
			t co; t so;
			if (max(abs(a[k][k]), abs(a[l][k])) < eps) { co = 1;  so = 0; }
			else {
				co = a[k][k] / (sqrt(a[k][k] * a[k][k] + a[l][k] * a[l][k]));
				so = a[l][k] / (sqrt(a[k][k] * a[k][k] + a[l][k] * a[l][k]));
			}
			c[k][k] = co;
			c[k][l] = so;
			c[l][k] = -so;
			c[l][l] = co;

			//меняем(умножаем) только k и l строки матрицы 
			for (int j = 0; j < n; j++) {
				e[k][j] = 0; a[k][j] = 0; e[l][j] = 0; a[l][j] = 0;
				
			}
			
			
			/*for (int j = 0; j < n; j++) {b[j] = 0;for (int z = 0; z < n; z++) { b[j] += c[z][j] * orig3[j]; }}*/
		}
	}
	t sum = 0.0;
	//R = TA, R - верхне треугольная

	//проверка на невырожденность(единственность решения)
	for (int i = 0; i < n; i++) {
		sum = 0.0;
		for (int j = 0; j < n; j++) {
			sum += a[i][j];
		}
		if ((abs(sum) < eps) and (param == 1)) { cout << endl << endl << "Not one solution!" << endl; }
	}

	//Qb*=b  ->  b* = Tb
	for (int j = 0; j < n; j++) {
		b[j] = 0;
		for (int z = 0; z < n; z++) { b[j] += e[j][z] * b1[z]; }
	}

	vector<t>x = inverse_step(n, a, b);//Rx=b*
	if (param == 0) { return x; }

	cout << endl << "R:";
	print(n, a);

	//Q=transp(T), Q - ортогональная
	vector<vector<t>>q = e;
	trans(n, q);
	cout << "Q:";
	print(n, q);

	cout << endl << "Result: " << endl;
	for (int i = 0; i < n; i++) {
		cout << x[i] << endl;
	}

	//ищем вектор невязки
	if (param) { disjoint(n, a1, b1, x); }
	return x;
}
template<typename t>
void cond(int n, vector <vector<t>> a) {
	vector <vector<t>>e(n, vector<double>(n));
	one(n, e);
	vector <vector<t>>a1(n, vector<double>(n));
	for (int i = 0; i < n; i++) {
		a1[i] = gaus(4, 4, a, e[i], 0);
	}
	cout << "Inverse matrix:" << endl;
	print(n, a1);
	t ans = norm1(4, a) * norm1(4, a1);
	cout << endl << "condA " << ans << endl;
	ans = norm2(4, a) * norm2(4, a1);
	cout << "condA " << ans << endl << endl;

	cout << "A*inverse(A):" << endl;
	vector<vector<t>>pr(n, vector<double>(n));
	mult(n, a1, a, pr);
	print(n, pr);
}
template<typename t>
void cond2(int n, vector <vector<t>> a, vector<t>b) {
	t eps2 = 0.01;
	vector<t>xg = gaus(n, n, a, b, 0);
	vector<t>xq = qr(n, n, a, b, 0);
	vector<t>lb(n, eps);
	vector<t>lxg = gaus(n, n, a, lb, 0);
	vector<t>lxq = qr(n, n, a, lb, 0);
	t dxg = norm0(n, lxg) / norm0(n, xg);
	t dxq = norm0(n, lxq) / norm0(n, xq);
	t db = norm0(n, lb) / norm0(n, b);
	t maxg = dxg / db;
	t maxq = dxq / db;
	for (int i = 0; i < 100; i++) {
		eps2 += 0.001;
		for (int j = 0; j < n; j++) { lb[j] = eps2; }
		lxg = gaus(n, n, a, lb, 0);
		lxq = qr(n, n, a, lb, 0);
		dxg = norm0(n, lxg) / norm0(n, xg);
		dxq = norm0(n, lxq) / norm0(n, xq);
		db = norm0(n, lb) / norm0(n, b);
		if (dxg / db > maxg) { maxg = dxg / db; }
		if (dxq / db > maxq) { maxq = dxq / db; }
	}
	cout << "CondA gause >= " << maxg << endl;
	cout << "CondA qr >= " << maxq << endl;

}
template<typename t>
void influence1(int n, vector <vector<t>> a, vector<t>b)
{
	t eps2 = 0.01;
	vector<t>xg = gaus(n, n, a, b, 0);
	vector<t>xq = qr(n, n, a, b, 0);
	vector<t>lb(n, 0);
	vector<t>max(n, 0);
	//vector<t>lxg = gaus(n, n, a, lb, 0);
	//vector<t>lxq = qr(n, n, a, lb, 0);
	for (int k = 0;k < n;k++)
	{
		
		vector<t>lxg = gaus(n, n, a, lb, 0);
		t dxg = norm0(n, lxg) / norm0(n, xg);
		t db = norm0(n, lb) / norm0(n, b);
		t maxg = 0;
		for (int i = 0; i < 100; i++) {
			eps2 += 0.01;
			lb[k] = eps2;
			for (int j = 0; j < n; j++) {
				lxg = gaus(n, n, a, lb, 0);
				dxg = norm0(n, lxg) / norm0(n, xg);
				db = norm0(n, lb) / norm0(n, b);
				if (dxg / db > maxg) 
				{	 
					 maxg = dxg / db;
					 max[k] = maxg;
				}
				
			}
		}
		lb[k] = 0;
	}
	t maxim = 0.0;
	t index;
	for (int i = 0;i < n;i++)
	{
		if (max[i] > maxim)
		{
			maxim = max[i];
			index = i;
		}
	}
	cout << index << endl;
	/*cout << "CondA gause >= " << maxg << endl;
	cout << "CondA qr >= " << maxq << endl;*/
}
void read_double(int& n, vector<vector<double>>& a, vector<double>& b, int k) {
	ifstream f1;
	f1.open("C:\\Users\\vfrcb\\Downloads\\1.txt");
	/*if (k > 0) {
	for (int i = 0; i < k; i++) {
	f1 >> n;
	f1.seekg((n*(n+1)),ios::cur);
	}
	}*/
	for (int u = 0; u <= k; u++) {
		f1 >> n;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				f1 >> a[i][j];
			}
		}
		for (int j = 0; j < n; j++) {
			f1 >> b[j];
		}
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
	vector<vector<double>>test1(4, vector<double>(4));
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i <= j) { test1[i][j] = 1; }
			else { test1[i][j] = 0; }
		}
	}
	vector<double>b1{ 4.0, 3.0, 2.0, 1.0 };

	vector<vector<double>>test2(4, vector<double>(4));
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i >= j) { test2[i][j] = 1; }
			else { test2[i][j] = 0; }
		}
	}

	vector<double>b2{ 1.0, 2.0, 3.0, 4.0 };
	vector<vector<double>>test21{ {0,0, 0, 1.0},{0,0, 1.0, 1.0}, {0,1.0,1.0,1.0},{1.0,1.0,1.0,1.0} };
	vector<vector<double>>test3{ {1.0,1.0, 1.0, 1.0},{2.0,3.0, 3.0, 3.0}, {2.0,4.0,4.0,4.0},{4.0,5.0,6.0,7.0} };
	vector<double>b3{ 4.0, 11.0, 15.0, 22.0 };
	vector<vector<double>>test4{ {10.0,6.0, 2.0, 0.0},{5.0,1.0, -2.0, 4.0}, {3.0,5.0,1.0,-1.0},{0.0,6.0,-2.0,2.0} };
	vector<double>b4{ 25.0, 14.0, 10.0, 8.0 };
	vector<vector<double>>test5{ {28.859,-0.008, 2.406, 19.24},{14.436,-0.001, 1.203, 9.624}, {120.204,-
	0.032,10.024,80.144},{-57.714,0.016,-4.812,-38.478} };
	vector<double>b5{ 30.459, 18.248, 128.156, -60.908 };



	vector<vector<double>>test111(4, vector<double>(4));
	vector<double>b111(4);
	int m = 4;
	//read_double(m, test111, b111, 3);
	gaus(4, 4, test4, b4, 1);


	//cout << "QR" << endl;

	qr(4, 4, test1, b1,1);
	qr(4, 4, test21, b2, 1);
	qr(4, 4, test3, b3, 1);
	qr(4, 4, test4, b4, 1);
	qr(4, 4, test5, b5, 1);
	//cond(4, test4);
	//cond2(4, test4, b4);
	//influence1(4, test4, b4);
}

