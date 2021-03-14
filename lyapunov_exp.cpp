#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
using namespace std;

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <mpreal.h>
using namespace mpfr;

#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <sys/time.h>

typedef boost::numeric::ublas::matrix<mpreal> matrix;
typedef boost::numeric::ublas::vector<mpreal> vector_col;
typedef boost::numeric::ublas::matrix<double> matrix_double;

// Точность по степенному ряду
#define eps_pw   "1e-20"

// Количество бит под мантиссу вещественного числа
#define b_m      128

// Количество шагов алгоритма для определения показателей Ляпунова
#define M        20000

// Для генерации псевдослучайных чисел
#define NUM_GEN  35

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
	vector_col Psi(2 * m);
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
			for(int p = 0; p < 2 * m; p++)
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

void orthogonal_norm(vector<vector_col> &arr_X, vector_col &norm_arr_X)
{
	mpreal v;
	vector<vector_col> orth_norm_arr_X;
	for(int i = 0; i < m; i++)
	{
		vector_col sum_vect = arr_X[i];
		for(int j = 0; j < i; j++)
		{
			v = inner_prod(arr_X[i], orth_norm_arr_X[j]);
			sum_vect -= v * orth_norm_arr_X[j];
		}
		v = boost::numeric::ublas::norm_2(sum_vect);
		norm_arr_X(i) = v;
		sum_vect /= v;
		orth_norm_arr_X.push_back(sum_vect);
	}
	for(int i = 0; i < m; i++)
		arr_X[i] = orth_norm_arr_X[i];
}

void calc_new_perturbation(vector_col &X, const vector_col &X_base,
                           vector_col &X_perturb, const mpreal &tau, int way = 1)
{
	for(int i = 0; i < m; i++)
	{
		X(i)     = X_base(i);
		X(i + m) = X_perturb(i);
	}

	calc(X, tau, 1);

	for(int i = 0; i < m; i++)
		X_perturb(i) = X(i + m);
}

string toString(int t)
{
	ostringstream s;
	s << t;
	return s.str();
}

void lyapunov_exp(vector_col &X, const mpreal &T)
{
	mpreal tau = T / M;

	vector_col lambda(m), X_base(m);
	vector<vector_col> X_perturb(m);
	ofstream *fL = new ofstream[m];
	string fn = "L";
	for(int p = 0; p < m; p++)
	{
		lambda(p) = 0;
		X_base(p) = X(p);
		X_perturb[p].resize(m);
		for(int l = 0; l < m; l++)
			X_perturb[p](l) = (rand() % NUM_GEN) / (double)(NUM_GEN + 1);

		fL[p].open((fn + toString(p + 1) + ".txt").c_str());
 		fL[p] << "set term eps\nset output \"L" << p + 1 << ".eps\"\n";
		fL[p] << "set xlabel \"t\"\nset ylabel \"lambda" << p + 1 << "\"" << endl;
		fL[p] << "plot '-' with line\n";
	}

	vector_col norm_arr_X(m);
	orthogonal_norm(X_perturb, norm_arr_X);
	ofstream f_out("lambda.txt");
	for(int k = 1; k <= M; k++)
	{
		for(int p = 0; p < m; p++)
			calc_new_perturbation(X, X_base, X_perturb[p], tau);
		for(int p = 0; p < m; p++)
			X_base(p) = X(p);

		orthogonal_norm(X_perturb, norm_arr_X);
		for(int p = 0; p < m; p++)
		{
			lambda(p) += log(norm_arr_X(p));
			fL[p] << k * tau << " " << (lambda(p) / (k * tau)).toString(10) << endl;
		}
		cout << "t = " << k * tau << " : lambda = " << lambda / (k * tau) << endl;
		f_out << "t = " << k * tau << " : lambda = " << lambda / (k * tau) << endl;
	}
	f_out.close();
	for(int p = 0; p < m; p++)
		fL[p].close();
	delete []fL;
}

unsigned long GetTickCount()
{
	timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

int main()
{
	srand(time(NULL));
	mpreal::set_default_prec(b_m);
	cout << "Машинный эпсилон = " << machine_epsilon() << endl;

	mpreal T;
	cout << "\nВведите длину отрезка времени > ";
	cin >> T;

	ifstream f_ini_val("initial_vals.txt");
	if(!f_ini_val)
		throw "Input File initial_vals.txt is not opened\n";
	f_ini_val >> m;
	vector_col X(2 * m);
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
	A.resize(2 * m, 2 * m);
	for(int i = 0; i < 2 * m; i++)
		for(int j = 0; j < 2 * m; j++)
			if(i < m && j < m)
				A(i, j) = A_double(i, j);
			else if(i >= m && j >= m)
				A(i, j) = A_double(i - m, j - m);
			else
				A(i, j) = 0;

	ifstream matrices_Q("Q.txt");
	if(!matrices_Q)
		throw "Input File Q.txt is not opened\n";
	matrix_double Q_add(2 * m, 2 * m);
	double maxnormQ, normQ;
	for(int l = 0; l < 2 * m; l++)
	{
		if(l < m)
		{
			matrix_double Qp;
			matrices_Q >> Qp;
			for(int i = 0; i < 2 * m; i++)
				for(int j = 0; j < 2 * m; j++)
					if(i < m && j < m)
						Q_add(i, j) = Qp(i, j);
					else
						Q_add(i, j) = 0;
		}
		else
			for(int i = 0; i < 2 * m; i++)
				for(int j = 0; j < 2 * m; j++)
					if(i < m && j >= m)
						Q_add(i, j) = Q[l - m](i, j - m).toDouble() + Q[l - m](j - m, i).toDouble();
					else
						Q_add(i, j) = 0;
		Q.push_back(Q_add);
		if(!l)
			maxnormQ = boost::numeric::ublas::norm_1(Q_add);
		else
		{
			normQ = boost::numeric::ublas::norm_1(Q_add);
			if(normQ > maxnormQ)
				maxnormQ = normQ;
		}
	}

	double nA  = boost::numeric::ublas::norm_1(A_double);
	mu         = 2 * m * maxnormQ;
	sum_nA_mu  = nA + mu;
	sum_nA_2mu = sum_nA_mu + mu;

	unsigned long t1 = GetTickCount();
	lyapunov_exp(X, T);
	cout << "\nTime = " << (GetTickCount() - t1) / 60000.0 << " min.\n\n";

	return 0;
}
