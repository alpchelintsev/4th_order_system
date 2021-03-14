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
#define eps_pw  "1e-30"

// Количество бит под мантиссу вещественного числа
#define b_m     128

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
	mpreal t = 0, delta_t, L, pr, tmin, tmax, dmin, dmax, tdmin, tdmax;
	int l = 0, imax, imin, lmin, lmax = 0, lmindt = 0, lmaxdt = 0;
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

		l++;
		bool flagimin = false;
		if(l == 1)
		{
			imax = 0;
			flagimin = true;
			tdmax = tdmin = tmax = 0;
			dmax = dmin = delta_t;
		}
		else if(delta_t < dmin)
		{
			dmin = delta_t;
			tdmin = t;
			lmindt = l;
		}
		else if(delta_t > dmax)
		{
			dmax = delta_t;
			tdmax = t;
			lmaxdt = l;
		}

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
		if(i >= imax)
		{
			imax = i;
			lmax = l;
			tmax = t-delta_t;
		}
		if(flagimin)
		{
			imin = i;
			lmin = l;
			tmin = 0;
		}
		else if(i < imin)
		{
			imin = i;
			lmin = l;
			if(t >= T)
				tmin = T-delta_t;
			else
				tmin = t-delta_t;
		}
	}
	while(fl_rp);

	cout << "\nКоординаты в конечный момент времени:\n";
	for(int p = 0; p < m; p++)
		cout << "X(" << p + 1 << ") = " << X(p).toString() << endl;

	cout << "\nmin Степень полиномов = " << imin << endl;
	cout << "Соответствующий момент времени = " << tmin << endl;
	cout << "Разница = " << T-tmin << endl;
	cout << "Номер момента времени = " << lmin << endl;
	cout << "\nmax Степень полиномов = " << imax << endl;
	cout << "Соответствующий момент времени = " << tmax << endl;
	cout << "Разница = " << T-tmax << endl;
	cout << "Номер момента времени = " << lmax << endl;
	cout << "\nmin delta_t = " << dmin << endl;
	cout << "Соответствующий момент времени = " << tdmin << endl;
	cout << "Разница = " << T-tdmin << endl;
	cout << "Номер момента времени = " << lmindt << endl;
	cout << "\nmax delta_t = " << dmax << endl;
	cout << "Соответствующий момент времени = " << tdmax << endl;
	cout << "Разница = " << T-tdmax << endl;
	cout << "Номер момента времени = " << lmaxdt << endl;

	cout << "\nN = " << l << endl;
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

	unsigned long t1 = GetTickCount();
	calc(X, T);
	//calc(X, T, -1);
	cout << "\nTime = " << (GetTickCount() - t1) / 60000.0 << " min.\n\n";

	return 0;
}
