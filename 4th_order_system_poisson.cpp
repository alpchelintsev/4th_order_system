#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <mpreal.h>
using namespace mpfr;

#include <sys/time.h>

typedef boost::numeric::ublas::matrix<mpreal> matrix;
typedef boost::numeric::ublas::vector<mpreal> vector_col;
typedef boost::numeric::ublas::matrix<double> matrix_double;

// Точность по степенному ряду
#define eps_pw  "1e-20"

// Количество бит под мантиссу вещественного числа
#define b_m     128

// Шаг по времени для отслеживания устойчивости по Пуассону
#define Delta_tau 1e-7

int m;
matrix A;
vector<matrix> Q;
double mu, sum_nA_2mu, sum_nA_mu;

// Функция возвращает длину отрезка сходимости степенного ряда
mpreal get_delta_t(const vector_col &Upsilon0)
{
	mpreal h1 = boost::numeric::ublas::norm_1(Upsilon0);
	mpreal h2 = h1 > 1 ? mu * h1 * h1 + sum_nA_2mu * h1 : sum_nA_mu;
	return 1 / (h2 + "1e-10");
}

// Функция вычисления значений фазовых координат в конечный момент времени
// T - длина отрезка интегрирования;
// way - направление поиска решений: 1 - вперед по времени, -1 - назад по времени
void calc(vector_col &X, const mpreal &T, int way = 1)
{
	mpreal t = 0, delta_t, L, pr;
	bool fl_rp;
	vector_col Psi(m);
	do
	{
		delta_t = get_delta_t(X);
		t += delta_t;
		if(t < T)
			fl_rp = true;
		else if(t > T)
		{
			delta_t -= t-T;
			fl_rp = false;
		}
		else
			fl_rp = false;

		vector<vector_col> Upsilon;
		Upsilon.push_back(X);

		L = boost::numeric::ublas::norm_1(X);
		int i = 0;
		pr = way * delta_t;
		while(L > eps_pw)
		{
			// Вычисляем новые коэффициенты степенного ряда
			mpreal sum_prod;
			for(int p = 0; p < m; p++)
			{
				sum_prod = 0;
				for(int j = 0; j <= i; j++)
				{
					vector_col prod_matrix_vector = prod(Q[p], Upsilon[j]);
					sum_prod += inner_prod(prod_matrix_vector, Upsilon[i - j]);
				}
				Psi(p) = sum_prod;
			}
			Upsilon.push_back((prod(A, Upsilon[i]) + Psi) / (i + 1));

			i++;

			X += Upsilon[i] * pr;
			L = fabs(pr) * boost::numeric::ublas::norm_1(Upsilon[i]);
			pr *= way * delta_t;
		}
	}
	while(fl_rp);
}

inline mpreal dist(vector_col &X, vector_col &X0)
{
	return boost::numeric::ublas::norm_2(X - X0);
}

unsigned long GetTickCount()
{
	timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

int main()
{
	mpreal::set_default_prec(b_m);
	cout << "Машинный эпсилон = " << machine_epsilon() << endl;

	mpreal T;
	cout << "\nВведите длину отрезка времени > ";
	cin >> T;

	ifstream f_ini_val("initial_vals.txt");
	if(!f_ini_val)
		throw "Input File initial_vals.txt is not opened\n";
	f_ini_val >> m;
	vector_col X(m);
	string s;
	for(int p = 0; p < m; p++)
	{
		f_ini_val >> s;
		X(p) = s;
	}
	f_ini_val.close();

	ifstream matrix_A("A.txt");
	if(!matrix_A)
		throw "Input File A.txt is not opened\n";
	matrix_double A_double;
	matrix_A >> A_double;
	matrix_A.close();
	A = A_double;

	ifstream matrices_Q("Q.txt");
	if(!matrices_Q)
		throw "Input File Q.txt is not opened\n";
	double maxnormQ, normQ;
	for(int p = 0; p < m; p++)
	{
		matrix_double Qp;
		matrices_Q >> Qp;
		Q.push_back(Qp);
		if(!p)
			maxnormQ = boost::numeric::ublas::norm_1(Qp);
		else
		{
			normQ = norm_1(Qp);
			if(normQ > maxnormQ)
				maxnormQ = normQ;
		}
	}

	double nA  = boost::numeric::ublas::norm_1(A_double);
	mu         = m * maxnormQ;
	sum_nA_mu  = nA + mu;
	sum_nA_2mu = sum_nA_mu + mu;

	mpreal t = 3 * Delta_tau;
	vector_col X0 = X;
	calc(X, Delta_tau);
	vector_col X1 = X;
	calc(X, Delta_tau);
	vector_col X2 = X;
	calc(X, Delta_tau);
	vector_col X3 = X;

	unsigned long t1 = GetTickCount();
	mpreal eps1 = dist(X1, X0), eps2 = dist(X2, X0), eps3 = dist(X3, X0);
	ofstream points("points.txt");
	while(t <= T)
	{
		if(eps1 > eps2 && eps2 < eps3 && eps2 < 1)
		{
			points << "t = " << t - Delta_tau << endl;
			points << "eps = " << eps2 << endl;
			points << "X(t) = " << X2 << endl << endl;
			cout << "t = " << t - Delta_tau << endl;
			cout << "eps = " << eps2 << endl;
			cout << "X(t) = " << X2 << endl << endl;
		}

		X1 = X2;
		eps1 = eps2;
		X2 = X3;
		eps2 = eps3;
		calc(X3, Delta_tau);
		eps3 = dist(X3, X0);

		t += Delta_tau;
	}
	cout << "\nTime = " << (GetTickCount() - t1) / 60000.0 << " min.\n\n";
	points << "Time = " << (GetTickCount() - t1) / 60000.0 << " min.\n\n";
	points.close();

	return 0;
}
