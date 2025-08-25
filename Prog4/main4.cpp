#include <fstream> 
#include <iostream>
#include <time.h>
#include <iomanip>//для вывода красивого в файлики
#include <cmath>
#define EPS 1e-6
#define MAX_ITER 100
using namespace std;

void pogreshnost1(double& e1_h05_abs, double& e1_h05_othn, double& e2_h05_abs, double& e2_h05_othn, double& einf_h05_abs, double& einf_h05_othn, int n, double* y, double* y05);

//для линейной задачи                                                                 
double f(double x, double y, double dy) {
	return dy * cos(x) + y * pow(x, 3) + tan(x);
}
void progonka(double* x, double* y, double h, double q_a, double q_b, int n, double* Mat1, double* Mat2, double* Mat3, double* Mat_prav, double* A, double* B) {
	{
		//граничные усл
		Mat1[0] = 1, Mat1[n - 1] = 1;//на главную диагональ
		Mat2[0] = 0;//первый эл-т верхней диагонали
		Mat3[n - 2] = 0;//последний эл-т нижней диагонали
		Mat_prav[0] = q_a;//гр условия y в а
		Mat_prav[n - 1] = q_b;//гр усл в b
		A[0] = 0;//прогоночные коэффициенты посчитаны на листочке
		B[0] = q_a;

		//формируем матрицу
		for (int i = 1; i < n - 1; i++) {//для внутренних узлов система
			Mat2[i] = 1 - cos(x[i]) * h / 2.0;//верхняя
			Mat3[i - 1] = 1 + cos(x[i]) * h / 2.0;//нижняя
			Mat1[i] = -(2 + pow(x[i], 3) * pow(h, 2));
			Mat_prav[i] = tan(x[i]) * pow(h, 2);
		}

		for (int i = 1; i < n - 1; i++) { // прямая прогонка(в целом n-1 коэф-ты не нужны мы же знаемзначение у(n-1)
			A[i] = -Mat2[i] / (Mat1[i] + Mat3[i - 1] * A[i - 1]);
			B[i] = (Mat_prav[i] - Mat3[i - 1] * B[i - 1]) / (Mat1[i] + Mat3[i - 1] * A[i - 1]);
		}
		// обратая прогонка поднимаемся вверх от n-2 
		y[n - 1] = q_b;//y(n-1) мы уже знаем
		for (int i = n - 2; i >= 0; i--){
			y[i] = B[i] + A[i] * y[i + 1];
		}
	}
}
///////////Рунге-Кутта 4го порядка для линейной задачи/////////////////
void Runge_4(double *x, double *y, double *dy, double h, int n) {
	double k1 = 0.0, k_1 = 0.0, k2 = 0.0, k_2 = 0.0, k3 = 0.0, k_3 = 0.0, k4 = 0.0, k_4 = 0.0;
	for (int i = 0; i < n - 1; i++) {
		k1 = h * dy[i];
		k_1 = h * f(x[i], y[i], dy[i]);

		k2 = h * (dy[i] + k_1 / 2.0);
		k_2 = h * f(x[i] + h / 2.0, y[i] + k1 / 2.0, dy[i] + k_1 / 2.0);

		k3 = h * (dy[i] + k_2 / 2.0);
		k_3 = h * f(x[i] + h / 2.0, y[i] + k2 / 2.0, dy[i] + k_2 / 2.0);

		k4 = h * (dy[i] + k_3);
		k_4 = h * f(x[i] + h, y[i] + k3, dy[i] + k_3);

		dy[i + 1] = dy[i] + 1.0 / 6.0 * (k_1 + 2.0 * k_2 + 2.0 * k_3 + k_4);
		y[i + 1] = y[i] + 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
	}
}
void strelba(double* x, double* y, double* dy, double* y_1, double* y_2, double h, int n, double q_a, double q_b) {//y_1 и y_2 -для первого и второго решений коши
	int i = 0;
	double t1 = 1.0;//параметр для 1го пристрела
	double t2 = 2.0;
	double l1 = 0.0;//здесь будет лежать значение полученного из задачи Коши у в точке b
	double l2 = 0.0;
	y_1[0] = q_a;//решение 1го прицела краевое знаечние на левом конце
	y_2[0] = q_a;
	dy[0] = t1;//формируем задачу Коши даём ну для y' в точке а
	Runge_4(x, y_1, dy, h, n);
	l1 = y_1[n - 1];//значение последнего у в 1 пристреле
	dy[0] = t2;
	Runge_4(x, y_2, dy, h, n);
	l2 = y_2[n - 1];//значение последнего у во 2 пристреле
	for (i = 0; i < n; i++){
		y[i] = ((q_b - l2) * y_1[i]) / (l1 - l2) + ((l1 - q_b) * y_2[i]) / (l1 - l2);
	}
}

//для нелинейной задачи                                                                
//double f2(double x, double y, double dy) {
	//return dy * cos(x) + pow(y, 2) * pow(x, 3) + tan(x);
//}
double f2(double x, double y, double dy) {
	return pow(dy,2) * cos(x) + pow(y, 2) * pow(x, 3) + tan(x);
}
///////////Рунге-Кутта 4го порядка для нелинейной задачи/////////////////
void Runge_42(double *x, double *y, double *dy, double h, int n) {
	double k1 = 0.0, k_1 = 0.0, k2 = 0.0, k_2 = 0.0, k3 = 0.0, k_3 = 0.0, k4 = 0.0, k_4 = 0.0;
	for (int i = 0; i < n - 1; i++) {
		k1 = h * dy[i];
		k_1 = h * f2(x[i], y[i], dy[i]);

		k2 = h * (dy[i] + k_1 / 2.0);
		k_2 = h * f2(x[i] + h / 2.0, y[i] + k1 / 2.0, dy[i] + k_1 / 2.0);

		k3 = h * (dy[i] + k_2 / 2.0);
		k_3 = h * f2(x[i] + h / 2.0, y[i] + k2 / 2.0, dy[i] + k_2 / 2.0);

		k4 = h * (dy[i] + k_3);
		k_4 = h * f2(x[i] + h, y[i] + k3, dy[i] + k_3);

		dy[i + 1] = dy[i] + 1.0 / 6.0 * (k_1 + 2.0 * k_2 + 2.0 * k_3 + k_4);
		y[i + 1] = y[i] + 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
	}
}

int Sekush(double* x, double* y, double* dy, double* Y, double* U , double h, int n, double q_a, double q_b) {
	int iter = 0;
	//double s = -3.1;
	//double s = 2.9;
	//double s = 1.8;
	double s = 0.1;
	double qb_iter = y[n - 1];//зачение на правом конце
	double delta = 1e-5;//в знаменателе производной
	double fi_strih;//производная

	while (fabs(qb_iter - q_b) > EPS && iter < MAX_ITER) {//пока значение на праом конце не станет боизко к нужному с достаточной точность.
		   iter++;
		   dy[0] = s;//недостающее начальное значение
		   y[0] = q_a;
		   Runge_42(x, y, dy, h, n);
		   qb_iter = y[n - 1];
		   U[0] = s + delta;//недостащее нач условие
		   Y[0] = q_a;
	       Runge_42(x, Y, U, h, n);
		   fi_strih = (Y[n - 1] - y[n - 1])/delta;
		   s = s - (y[n - 1] - q_b) / fi_strih;
	}
	return iter;
}
//double f_y_dy(double x, double y, double dy, double Y, double U) {
	//return 2 * y * pow(x, 3) * Y +  cos(x) * U;
//}

double f_y_dy(double x, double y, double dy, double Y, double U) {//для второй задачи Коши в методе Ньютона
return 2 * y * pow(x, 3) * Y + 2 * dy * cos(x) * U;
}

void Runge_43(double* x, double* y, double* dy, double* Y, double* U, double h, int n) {
	double k1 = 0.0, k_1 = 0.0, k2 = 0.0, k_2 = 0.0, k3 = 0.0, k_3 = 0.0, k4 = 0.0, k_4 = 0.0;
	for (int i = 0; i < n - 1; i++) {
		k1 = h * U[i];
		k_1 = h * f_y_dy(x[2*i], y[2 * i], dy[2 * i], Y[i], U[i]);

		k2 = h * (U[i] + k_1 / 2.0);
		k_2 = h * f_y_dy(x[2*i] + h / 2.0, y[i * 2 + 1], dy[i * 2 + 1], Y[i] + k1 / 2.0, U[i] + k_1 / 2.0);

		k3 = h * (U[i] + k_2 / 2.0);
		k_3 = h * f_y_dy(x[2*i] + h / 2.0, y[i * 2 + 1], dy[i * 2 + 1], Y[i] + k2 / 2.0, U[i] + k_2 / 2.0);

		k4 = h * (U[i] + k_3);
		k_4 = h * f_y_dy(x[2*i] + h, y[i * 2 + 2], dy[i * 2 + 2], Y[i] + k3, U[i] + k_3);

		U[i + 1] = U[i] + 1.0 / 6.0 * (k_1 + 2.0 * k_2 + 2.0 * k_3 + k_4);
		Y[i + 1] = Y[i] + 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
	}
}
int Newtoun(double* x, double* y, double* dy, double* Y, double* U, double h, int n, double q_a, double q_b) {
	int iter = 0;
	//double s = 2.9;
	//double s = 1.8;
	double s = 0.1;//выбирается значение производной на левом конце
	double qb_iter = y[n - 1];

	while (fabs(qb_iter - q_b) > EPS && iter < MAX_ITER) {//пока значение на правом конце рещения не будет меньше точности
		iter++;
		dy[0] = s;//кр.усл.
		y[0] = q_a;
		Y[0] = 0;
		U[0] = 1;
		Runge_42(x, y, dy, h, n);//рунге-кутта для нелинейной задачи
		qb_iter = y[n - 1];//значение у на правом конце
		double h2 = 2 * h;
		int n2 = n / 2.0;
		Runge_43(x, y, dy,Y,U, h2, n2);
		s = s - (y[n - 1] - q_b) / Y[n2-1];
	}
	return iter;
}

int main(){
	const int n = 20;
	const int n05 = (n - 1) * 2 + 1;
	double a = 0.0;
	double b = 1.0;
	double h = (b - a) / (n - 1);
	double h05 = (b - a) / (n05 - 1);
	double q_a = 0.5;//y(a)-краевае условия берем из отрезка[-1;1]
	//double q_a = 3;
	double q_b =0.1;//y(b)
	//double q_b = -2;

	///////////////////////////////////////конечно-разностный метод+метод прогонки/////////////////////////////////

	double* x = new double[n];//иксы пbо сетке h
	for (int i = 0; i < n; i++) {//равномерная сетка с шагом h
		x[i] = a + i * h;
	}
	double* y = new double[n];
	double* x05 = new double[n05];//для измельчённой сетки
	for (int i = 0; i < n05; i++) {//равномерная сетка с шагом h/2
		x05[i] = a + i * h05;
	}
	double* y05 = new double[n05];

	double* Mat1 = new double[n]; //средняя диагональ
	double* Mat2 = new double[n - 1]; //верхняя диагональ
	double* Mat3 = new double[n - 1]; //нижняя диагональ
	double* Mat_prav = new double[n]; //правые части
	double* A = new double[n];//прогоночные коэф-ты
	double* B = new double[n];

	double* Mat105 = new double[n05]; //средняя диагональ
	double* Mat205 = new double[n05 - 1]; //верхняя диагональ
	double* Mat305 = new double[n05 - 1]; //нижняя диагональ
	double* Mat_prav05 = new double[n05]; //правые части
	double* A05 = new double[n05];
	double* B05 = new double[n05];

	for (int i = 0; i < n; i++) {
		A[i] = 0.0;
		B[i] = 0.0;
	}
	for (int i = 0; i < n; i++) {
		A05[i] = 0.0;
		B05[i] = 0.0;
	}
	progonka(x, y, h, q_a, q_b, n, Mat1, Mat2, Mat3, Mat_prav, A, B);
	progonka(x05, y05, h05, q_a, q_b, n05, Mat105, Mat205, Mat305, Mat_prav05, A05, B05);
	cout << "phogonka_h" << endl;
	for (int i = 0; i < n; i++) {
		cout << y[i] << ' ';
	}
	cout << endl;
	cout << endl;
	cout << "phogonka_h/2" << endl;
	for (int i = 0; i < n05; i++) {
		cout << y05[i] << ' ';
	}
	cout << endl;
	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////метод стрельбы//////////////////////////////////
	double* y_1 = new double[n];
	double* y_2 = new double[n];
	double* y_str = new double[n];
	double* dy_str = new double[n];

	double* y05_1 = new double[n05];
	double* y05_2 = new double[n05];
	double* y05_str = new double[n05];
	double* dy05_str = new double[n05];

	strelba(x, y_str, dy_str, y_1, y_2, h, n, q_a, q_b);
	strelba(x05, y05_str, dy05_str, y05_1, y05_2, h05, n05, q_a, q_b);

	cout << endl;
	cout << "strelba_h" << endl;
	for (int i = 0; i < n; i++) {
		cout << y_str[i] << ' ';
	}
	cout << endl;
	cout << endl;
	cout << "strelba_h/2" << endl;
	for (int i = 0; i < n05; i++) {
		cout << y05_str[i] << ' ';
	}
	cout << endl;



	//////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////метод Ньютона///////////////////////////////////////
	int n2=n/2.0;
	//int n2 = (n - 1) / 2 + 1;
	double* y_N = new double[n];
	double* dy_N = new double[n];
	double* Y_N = new double[n2];//для вспомогательной задачи Коши
	double* U_N = new double[n2];
	double* y05_N = new double[n05];
	double* dy05_N = new double[n05];
	double* Y05_N = new double[n];
	double* U05_N = new double[n];//dY=U
	for (int i = 0; i < n; i++) {
		y_N[i] = 0.0;
		dy_N[i] = 0.0;
		Y05_N[i] = 0.0;
		U05_N[i] = 0.0;

	}
	for (int i = 0; i < n2; i++) {
		Y_N[i] = 0.0;
		U_N[i] = 0.0;
	}
	for (int i = 0; i < n05; i++) {
		y05_N[i] = 0.0;
		dy05_N[i] = 0.0;
	}

	int iter1 = 0;
	iter1 = Newtoun(x, y_N, dy_N, Y_N, U_N, h, n, q_a, q_b);

	cout << endl;
	cout << "Newtoun_h"<<endl;
	cout << "iter = " << iter1 << "\n";
	for (int i = 0; i < n; i++) {
		cout << y_N[i] << ' ';
	}
	cout << endl;

	int iter2 = 0;
	iter2 = Newtoun(x05, y05_N, dy05_N, Y05_N, U05_N, h05, n05, q_a, q_b);
	cout << endl;
	cout << "Newtoun_h/2"<<endl;
	cout << "iter = " << iter2 << "\n";
	for (int i = 0; i < n05; i++) {
		cout << y05_N[i] << ' ';
	}
	cout << endl;
	

	
	////////////////////////////////////Метод секущих//////////////////////////////////////
	double* y_sek = new double[n];
	double* dy_sek = new double[n];
	double* Y_sek = new double[n];
	double* U_sek = new double[n];
	double* y05_sek = new double[n05];
	double* dy05_sek = new double[n05];
	double* Y05_sek = new double[n05];
	double* U05_sek = new double[n05];

	for (int i = 0; i < n; i++){
		y_sek[i] = 0.0;
		dy_sek[i] = 0.0;
		Y_sek[i] = 0.0;
		U_sek[i] = 0.0;
	}
	for (int i = 0; i < n05; i++) {
		y05_sek[i] = 0.0;
		dy05_sek[i] = 0.0;
		Y05_sek[i] = 0.0;
		U05_sek[i] = 0.0;
	}

	int iter3 = 0;//счётчик итераций
	iter3 = Sekush(x, y_sek, dy_sek, Y_sek, U_sek, h, n, q_a, q_b);
	cout << endl;
	cout << "Sekush_h" << endl;
	cout << "iter = " << iter3 << "\n";
	for (int i = 0; i < n; i++) {
		cout << y_sek[i] << ' ';
	}
	cout << endl;

	int iter4 = 0;
	iter4 = Sekush(x05, y05_sek, dy05_sek, Y05_sek, U05_sek, h05, n05, q_a, q_b);
	cout << endl;
	cout << "Sekush_h/2" << endl;
	cout << "iter = " << iter4 << "\n";
	for (int i = 0; i < n05; i++) {
		cout << y05_sek[i] << ' ';
	}
	cout << endl;

///////////////////////////ПОГРЕШНОСТИ////////////////////////////////
	double e1_abs_1_05 = 0.0,  //  1 мтод прогонки e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2
		e1_otn_1_05 = 0.0,
		e2_abs_1_05 = 0.0,
		e2_otn_1_05 = 0.0,
		einf_abs_1_05 = 0.0,
		einf_otn_1_05 = 0.0,

		e1_abs_2_05 = 0.0,  //  2 метод стрельбы e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2
		e1_otn_2_05 = 0.0,
		e2_abs_2_05 = 0.0,
		e2_otn_2_05 = 0.0,
		einf_abs_2_05 = 0.0,
		einf_otn_2_05 = 0.0,

		e1_abs_3_05 = 0.0, //  3 метод Ньютона e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2
		e1_otn_3_05 = 0.0,
		e2_abs_3_05 = 0.0,
		e2_otn_3_05 = 0.0,
		einf_abs_3_05 = 0.0,
		einf_otn_3_05 = 0.0,

		e1_abs_4_05 = 0.0, //  4 метод секущих e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2
		e1_otn_4_05 = 0.0,
		e2_abs_4_05 = 0.0,
		e2_otn_4_05 = 0.0,
		einf_abs_4_05 = 0.0,
		einf_otn_4_05 = 0.0;
	//для y
	//между h и h/2	
	pogreshnost1(e1_abs_1_05, e1_otn_1_05, e2_abs_1_05, e2_otn_1_05, einf_abs_1_05, einf_otn_1_05, n, y, y05);
	pogreshnost1(e1_abs_2_05, e1_otn_2_05, e2_abs_2_05, e2_otn_2_05, einf_abs_2_05, einf_otn_2_05, n, y_str, y05_str);
	pogreshnost1(e1_abs_3_05, e1_otn_3_05, e2_abs_3_05, e2_otn_3_05, einf_abs_3_05, einf_otn_3_05, n, y_N, y05_N);
	pogreshnost1(e1_abs_4_05, e1_otn_4_05, e2_abs_4_05, e2_otn_4_05, einf_abs_4_05, einf_otn_4_05, n, y_sek, y05_sek);

	/////////////////////////////////////////////////////////////////////Файлы////////////////////////////////
	ofstream fout1;
	fout1.open("file1.txt");
	for (int i = 0; i < n; i++) {
		fout1 << fixed;
		fout1.precision(5);
		fout1.setf(ios::right);
		fout1.width(10);
		fout1 << x[i] << "\n";
	}
	fout1.close();
	ofstream fout2;
	fout2.open("file2.txt");
	for (int i = 0; i < n; i++) {
		fout2 << fixed;
		fout2.precision(5);
		fout2.setf(ios::right);
		fout2.width(10);
		fout2 << y[i] << "\n";
	}
	fout2.close();


	ofstream fout3;
	fout3.open("file3.txt");
	for (int i = 0; i < n; i++) {
		fout3 << fixed;
		fout3.precision(5);
		fout3.setf(ios::right);
		fout3.width(10);
		fout3 << y_str[i] << "\n";
	}
	fout3.close();

	ofstream fout4;
	fout4.open("file4.txt");
	for (int i = 0; i < n; i++) {
		fout4 << fixed;
		fout4.precision(5);
		fout4.setf(ios::right);
		fout4.width(10);
		fout4 << y_N[i] << "\n";
	}
	fout4.close();

	ofstream fout5;
	fout5.open("file5.txt");
	for (int i = 0; i < n; i++) {
		fout5 << fixed;
		fout5.precision(5);
		fout5.setf(ios::right);
		fout5.width(10);
		fout5 << y_sek[i] << "\n";
	}
	fout5.close();

	ofstream fout12;
	fout12.open("file12.txt");
	fout12 << "tabl for y, pogreshnost dlya h and h/2" << endl;
	fout12.setf(ios::scientific);

	fout12 << endl;
	fout12 << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(23) << left <<" "<< setw(30) << left << "|| . ||1" << setw(30) << left << "|| . || 2" << setw(30) << left << "|| . || inf " << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "metod" << setw(15) << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs " << setw(15) << left << "otn" << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "progonka" << setw(15) << left << e1_abs_1_05 << setw(15) << left << e1_otn_1_05 << setw(15) << left << e2_abs_1_05 << setw(15) << left << e2_otn_1_05 << setw(15) << left << einf_abs_1_05 << setw(15) << left << einf_otn_1_05 << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "strelba" << setw(15) << left << e1_abs_2_05 << setw(15) << left << e1_otn_2_05 << setw(15) << left << e2_abs_2_05 << setw(15) << left << e2_otn_2_05 << setw(15) << left << einf_abs_2_05 << setw(15) << left << einf_otn_2_05 << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "Newtoun" << setw(15) << left << e1_abs_3_05 << setw(15) << left << e1_otn_3_05 << setw(15) << left << e2_abs_3_05 << setw(15) << left << e2_otn_3_05 << setw(15) << left << einf_abs_3_05 << setw(15) << left << einf_otn_3_05 << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "sekushih" << setw(15) << left << e1_abs_4_05 << setw(15) << left << e1_otn_4_05 << setw(15) << left << e2_abs_4_05 << setw(15) << left << e2_otn_4_05 << setw(15) << left << einf_abs_4_05 << setw(15) << left << einf_otn_4_05 << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	
	fout12.close();
	




	delete[] Mat1;
	delete[] Mat2;
	delete[] Mat3;
	delete[] Mat_prav;
	delete[] Mat105;
	delete[] Mat205;
	delete[] Mat305;
	delete[] Mat_prav05;
	delete[] x;
	delete[] y;
	delete[] x05;
	delete[] y05;
	delete[] A;
	delete[] B;
	delete[] A05;
	delete[] B05;
	delete[] y_str;
	delete[] dy_str;
	delete[] y_1;
	delete[] y_2;
	delete[] y05_str;
	delete[] dy05_str;
	delete[] y05_1;
	delete[] y05_2;
	delete[] y_sek;
	delete[] dy_sek;
	delete[] Y_sek;
	delete[] U_sek;
	delete[] y05_sek;
	delete[] dy05_sek;
	delete[] Y05_sek;
	delete[] U05_sek;
	delete[] y_N;
	delete[] dy_N;
	delete[] Y_N;
	delete[] U_N;
	delete[] y05_N;
	delete[] dy05_N;
	delete[] Y05_N;
	delete[] U05_N;


	std::system("python PythonApplication1.py");
	return 0;
}
void pogreshnost1(double& e1_h05_abs, double& e1_h05_othn, double& e2_h05_abs, double& e2_h05_othn, double& einf_h05_abs, double& einf_h05_othn, int n, double* y, double* y05)
{
	for (int i = 0; i < n; i++)
	{
		e1_h05_abs += fabs(y[i] - y05[2 * i]);
		e1_h05_othn += fabs(y[i]);
		e2_h05_abs += pow(fabs(y[i] - y05[2 * i]), 2);
		e2_h05_othn += pow(fabs(y[i]), 2);
		einf_h05_abs = max(einf_h05_abs, fabs(y[i] - y05[2 * i]));
		einf_h05_othn = max(einf_h05_othn, fabs(y[i]));
	}

	e1_h05_othn = e1_h05_abs / e1_h05_othn;
	e2_h05_abs = sqrt(e2_h05_abs);
	e2_h05_othn = sqrt(e2_h05_othn);
	e2_h05_othn = e2_h05_abs / e2_h05_othn;
	einf_h05_othn = einf_h05_abs / einf_h05_othn;
}