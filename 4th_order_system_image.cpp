#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

#include <cstdlib>

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

// Шаг по времени для построения дуги траектории
#define Delta_tau 1e-4

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

int main()
{
	mpreal::set_default_prec(b_m);
	cout << "Машинный эпсилон = " << machine_epsilon() << endl;

	mpreal T;
	cout << "\nВведите длину отрезка времени > ";
	cin >> T;

	cout << "\nВведите проекцию > ";
	string sp;
	cin >> sp;
	int x1 = sp[0] - 48;
	int x2 = sp[1] - 48;
	int x3 = sp[2] - 48;

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

	string ranges[3][2]; int ind_range = 0;
	ifstream f_ranges("ranges.txt");
	if(!f_ranges)
		throw "Input File ranges.txt is not opened\n";
	string current_st1, current_st2;
	for(int p = 1; p <= m; p++)
	{
		f_ranges >> current_st1 >> current_st2;
		if(p == x1 || p == x2 || p == x3)
		{
			ranges[ind_range][0] = current_st1;
			ranges[ind_range][1] = current_st2;
			ind_range++;
		}
	}
	f_ranges.close();

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

	ofstream f("image_gnuplot.txt");
	if(!f)
		throw "Output File image_gnuplot.txt is not opened\n";
	f << "set term eps\nset output \"traj_img.eps\"\n";
	f << "set xlabel \"x" << x1 << "\"\nset ylabel \"x" << x2 <<
	     "\"\nset zlabel \"x" << x3 << "\"\n";
	f << "set xrange [" << ranges[0][0] << ":" << ranges[0][1] <<
	     "]\nset yrange [" << ranges[1][0] << ":" << ranges[1][1] <<
	     "]\nset zrange [" << ranges[2][0] << ":" << ranges[2][1] << "]\n";
	f << "splot '-' with line\n";

	x1--; x2--; x3--;
	f << X(x1).toString(4) << " " << X(x2).toString(4) << " " <<
	     X(x3).toString(4) << endl;
	mpreal t = 0;
	while(t <= T)
	{
		calc(X, Delta_tau);
		f << X(x1).toString(4) << " " << X(x2).toString(4) << " " <<
		     X(x3).toString(4) << endl;
		t += Delta_tau;
	}
	f.close();

	system("gnuplot image_gnuplot.txt");

	return 0;
}
