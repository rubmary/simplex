#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include <set>
#include <fstream>
#include <limits>
#include <vector>
#include <ostream>
#include <iomanip>
#define EPS 1-10
using namespace Eigen;
using namespace std;

void print_set(const set <int> &S, ostream *os = &cout)
{
	int T = S.size();
	for (auto x : S)
		*os << x << ' ';
	*os << endl;
}

void print_steps(
			bool print,
			int k,
			const set <int> &B,
			const set <int> &N,
			const MatrixXd &AB,
			const MatrixXd &AN,
			const VectorXd &cB,
			const VectorXd &cN)
{
	if (!print)
		return;
	
	for (int i = 0; i < 50; i++) cout << '_';
		cout << endl;
	for (int i = 0; i < 18; i++) cout << '#';
	cout << " Iteracion  " << k << ' ';
	for (int i = 0; i < 18; i++) cout << '#';
	cout << endl;
	for (int i = 0; i < 50; i++) cout << '-';
		cout << endl << endl;

	cout << "B = ";
	print_set(B);
	cout << "N = ";
	print_set(N);
	cout << endl;
	cout << "AB =" << endl;
	cout << AB << endl << endl;
	cout << "AN =" << endl;
	cout << AN << endl << endl;
	cout << "cB =" << endl;
	cout << cB << endl << endl;
	cout << "cN =" << endl;
	cout << cN << endl << endl;
	cout << endl;
}

void print_steps(
			bool print,
			double z,
			const VectorXd &xB,
			const VectorXd &y,
			const VectorXd &cr)
{

	if (!print)
		return;
	cout << "Calcular xB resolviendo AB*xB = b" << endl;
	cout << "xB =" << endl;
	cout << xB << endl << endl;
	cout << "Calcular funcion objetivo" << endl;
	cout << "z = " << z << endl << endl;
	cout << "Calcular y resolviendo yt*AB = cBt" << endl;
	cout << "y =" << endl;
	cout << y << endl << endl;
	cout << "Calcular costos reducidos cr = cN - ANt*y" << endl;
	cout << "cr = " << endl;
	cout << cr << endl << endl;
}

void print_steps( 
			bool print,
			int e,
			const VectorXd &dB,
			const VectorXd &tB) 
{
	if (!print)
		return;
	cout << "Calcular la variable entrante a la base" << endl;
	cout << "e = " << e << endl << endl;

	cout << "Calcular dB resolviendo AB*dB = -cole(A)" << endl;
	cout << "dB =" << endl;
	cout << dB << endl << endl;

	cout << "Calcular la cota del parametro t, tB" << endl;
	cout << "tB =" << endl;
	cout << tB << endl << endl;
}

void print_steps(
			bool &print,
			int e,
			int s)
{
	if (!print)
		return;
	cout << "Calcular la variable saliente de la base" << endl;
	cout << "s = " << s << endl << endl;
	cout << "Se pivotea x" << e << " con x" << s << endl; 

	for (int i = 0; i < 40; i++) cout << '-';
	cout << endl << endl;

	cout << "Siguiente iteracion? ";
	string answer;
	cin >> answer;
	if (answer[0] == 'n' || answer[0] == 'N')
		print = false;
}

MatrixXd extract(	const MatrixXd &A,
					const set<int> &indices)
{
	int M = A.rows();
	int N = indices.size();
	MatrixXd B(M, N);
	int j = 0;
	for (auto k : indices) {
		for (int i = 0; i < M; i++)
			B(i, j) = A(i, k-1);
		j++;
	}
	return B;
}

VectorXd extract(	const VectorXd &v, 
					const set<int> &indices )
{
	int N = indices.size();
	VectorXd x(N);
	int k = 0;
	for (auto i : indices)
		x[k++] = v[i-1];
	return x;
}

int variable_entrante(	const VectorXd &cr,
						const set<int> &indices )
{

	int k = 0;
	for (auto i : indices) {
		if(cr[k] > 0)
			return i;
		k++;
	}
	return -1;
}

VectorXd cotas(	const VectorXd &xB,
				const VectorXd &dB)
{
	int M = xB.size();
	VectorXd t(M);
	for (int i = 0; i < M; i++) {
		if (dB[i] >= 0 )
			t[i] = std::numeric_limits<double>::infinity();
		else
			t[i] = -xB[i]/dB[i];
	}
	return t;
}

int variable_saliente(	const VectorXd &dB,
						const VectorXd &tB,
						const set <int> &B)
{
	int s = -1;
	int k = 0;
	double min = std::numeric_limits<double>::infinity();
	for (auto j : B)
	{
		if (dB[k] < 0 && tB[k] < min) {
			s = j;
			min = tB[k];
		}
		k++;
	}
	return s;
}

set <int> simplex(
				const VectorXd &c,
				const VectorXd &b,
				const MatrixXd &A,
				bool print=false,
				set<int>B=set<int>(),
				set<int>N=set<int>())
{
	int m = b.size();
	int n = (int)c.size() - m;
	VectorXd xB, y, cB, cN, cr, dB, tB;
	MatrixXd AB, AN;
	int e, s;
	double z;

	if (B.empty()) {
		for (int j = 1; j <= n; j++)
			N.insert(j);
		for (int i = n+1; i <= n+m; i++)
			B.insert(i);
	} else {
		for (int k = 1; k <= n+m; k++) {
		if (B.find(k) == B.end())
			N.insert(k);
		}
	}

	for (int k = 1; ; k++) {
		AB = extract(A, B);
		AN = extract(A, N);
		cB = extract(c, B);
		cN = extract(c, N);
		print_steps(print, k, B, N, AB, AN, cB, cN);

		xB = AB.fullPivLu().solve(b);
		z = cB.transpose()*xB;
		y = AB.transpose().fullPivLu().solve(cB);
		cr = cN - AN.transpose()*y;
		print_steps(print, z, xB, y, cr);
		
		e = variable_entrante(cr, N);
		if (e == -1)
			return B;

		dB = AB.fullPivLu().solve(-A.col(e-1));
		tB = cotas(xB, dB);
		print_steps(print, e, dB, tB);
		
		s = variable_saliente(dB, tB, B);
		if (s == -1)
			return set<int>();
		print_steps(print, e, s);

		B.erase(s);  B.insert(e);
		N.erase(e); N.insert(s);
	}
	return set<int>();
}


bool answer (string ans) {
	return !(ans[0] == 'n' || ans[0] == 'N'); 
}

void print_solution (	
				const set <int> &B,
				const VectorXd &c,
				const VectorXd &b,
				const MatrixXd &A,
				string name)
{
	VectorXd xB, cB, cN, cr, dB, tB;
	MatrixXd AB, AN;
	double z;
	int m = b.size();
	int n = c.size() - m;

	VectorXd x(n+m), y(m);
	set <int> N;

	for (int k = 1; k <= n+m; k++) {
		if (B.find(k) == B.end())
			N.insert(k);
		x[k-1] = 0;
	}

	AB = extract(A, B);
	AN = extract(A, N);
	cB = extract(c, B);
	cN = extract(c, N);

	xB = AB.fullPivLu().solve(b);
	int k = 0;
	for (auto j : B)
		x[j-1] = xB[k++];
	z = cB.transpose()*xB;
	y = AB.transpose().fullPivLu().solve(cB);
	cr = cN - AN.transpose()*y;

	bool unique_solution = true;
	for (int j = 0; j < n; j++)
		if (abs(cr[j]) < EPS)
			unique_solution = false;

	ostream *os;
	if (name == "")
		os = &cout;
	else
		os = new ofstream(name.c_str());

	for (int i = 0; i < 50; i++) *os << '_';
		*os << endl;
	for (int i = 0; i < 20; i++) *os << '#';
	*os << " SOLUCION ";
	for (int i = 0; i < 20; i++) *os << '#';
	*os << endl;
	for (int i = 0; i < 50; i++) *os << '-';
		*os << endl << endl;
	if (unique_solution)
		*os << "Solucion unica" << endl;
	else
		*os << "Infinitas soluciones" << endl;
	*os << "Las variables basicas son" << endl;
	*os << "B = ";
	print_set(B, os);
	*os << "Las variables no basicas son" << endl;
	*os << "N = ";
	print_set(N, os);
	*os << endl;
	*os << "Solucion primal" << endl;
	*os << "x = "  << endl << x  << endl << endl;
	*os << "Solucion dual" << endl;
	*os << "y = "  << endl << y  << endl << endl;
	*os << "Costos reducidos" << endl;
	*os << "cr = " << endl << cr << endl << endl;
	*os << "Funcion objetivo" << endl;
	*os << "z* = ";
	*os << fixed << setprecision(5) << z << endl;
	for (int i = 0; i < 50; i++) *os << '=';
		*os << endl << endl;
}

int main(int argc, char **argv) {
	ifstream file(argv[1]);
	int N, M;				// variables y restricciones
	file >> N >> M;			// leer parametros
	VectorXd c(N+M), b(M);	// crear vector de costos y recursos
	MatrixXd A(M, N+M);		// crear matriz extendida
	
	// Leer el vector de costos
	for (int j = 0; j < N; j++)
		file >> c[j];
	for (int j = N; j < N+M; j++)
		c[j] = 0;

	// Leer el vector de recursos
	for (int i = 0; i < M; i++)
		file >> b[i];

	// Leer la matriz
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			file >> A(i, j);

	// Crear matriz extendida
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			if (i == j)
				A(i, j+N) = 1;
			else
				A(i, j+N) = 0;
		} 
	}

	cout << "El numero de variables es:     " << N << endl;
	cout << "El numero de restricciones es: " << M << endl;
	cout << "El vector de costos es:" << endl;
	cout << c << endl;

	cout << "El vector de recursos b es:" << endl;
	cout << b << endl;

	cout << "La matriz A es: " << endl;
	cout << A << endl;
	cout << endl;

	bool print;
	string ans;
	set <int> B;
	cout << "Utilizar base inicial? ";
	cin >> ans;
	if (!answer(ans)) {
		cout << "Ingrese "  << M << " numeros: ";
		int x;
		for (int j = 0; j < M; j++) {
			cin >> x;
			B.insert(x);
		}
	}
	cout << endl << "Desea mostrar cada iteracion? ";
	cin >> ans;
	cout << endl;
	print = answer(ans);
	B = simplex(c, b, A, print, B);

	while (true) {
		cout << "Elegir una opcion" << endl;
		cout << "1. Mostrar solucion en consola"  << endl;
		cout << "2. Imprimir solucion en archivo" << endl;
		cout << "3. Analisis de sensitividad"     << endl;
		cout << "4. Salir"                        << endl;
		cout << endl;
		cout << "Opcion [1/2/3/4]: ";
		cin >> ans;
		cout << endl;
		if (ans == "1" || ans == "2") {
			string name = "";
			if (ans == "2") {
				cout << "Ingresa el nombre del archivo: ";
				cin >> name;
				cout << endl;
			}
			print_solution(B, c, b, A, name);
		}
		else if (ans == "3") {
			cout << "No esta listo :(" << endl;
			cout << endl;
		} else if (ans == "4") {
			cout << "Adios!" << endl;
			break;
		} else
			cout << "Opcion no valida" << endl;
	}
}