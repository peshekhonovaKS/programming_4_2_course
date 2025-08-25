#include <fstream> 
#include <iostream>
#include <time.h>
#include <iomanip>//для вывода красивого в файлики
using namespace std;

double f(double x, double y, double dy){ //f(x,y,y')=y'cosx+yx^3+tgx
	return dy* cos(x) + y* pow(x, 3) + tan(x);
}

void pogreshnost1(double& e1_h05_abs, double& e1_h05_othn, double& e2_h05_abs, double& e2_h05_othn, double& einf_h05_abs, double& einf_h05_othn, int n, double* y, double* y05);
void pogreshnost2(double& e1_h2_abs, double& e1_h2_othn, double& e2_h2_abs, double& e2_h2_othn, double& einf_h2_abs, double& einf_h2_othn, int n, double* y, double* y2);

void Runge_2(double* x, double* y, double* dy, double h, int n) {
	double k1 = 0.0, k_1 = 0.0, k2 = 0.0, k_2 = 0.0;
	for (int i = 0; i < n - 1; i++) {
		k1 = h * dy[i];
		k_1 = h * f(x[i], y[i], dy[i]);

		k2 = h * (dy[i] + k_1);
		k_2 = h * f(x[i] + h, y[i] + k1, dy[i] + k_1);

		dy[i + 1] = dy[i] + (k_1 + k_2)/2.0;
		y[i + 1] = y[i] + (k1 + k2) /2.0;

	}
}
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
void Adams(double *x, double *y, double *dy, double h, int n) {
	for (int i = 2; i < n - 1; i++) {
		dy[i + 1] = dy[i] + h / 12.0 * (23.0 * f(x[i], y[i], dy[i]) - 16.0 * f(x[i - 1], y[i - 1], dy[i - 1]) + 5.0 * f(x[i - 2], y[i - 2], dy[i - 2]));
		y[i + 1] = y[i] + h / 12.0 * (23.0 * dy[i] - 16.0 * dy[i - 1] + 5.0 * dy[i - 2]);
	}
}
int main() {
	int n = 40;//для сетки h
	int n05 = (n - 1) * 2 + 1;//для сетки h/2
	int n2 = (n - 1) / 2 + 1;//для сетки 2h
	double a = 0.0;//отрезок начало
	double b = 1.0;//отрезок конец
	double h = (b - a) / (n - 1);
	double h05 = (b - a) / (n05 - 1);
	double h2 = (b - a) / (n2 - 1);
	//double q0 = -5.0;//y(0)
	//double q1 = 5.0;//dy(0)

	 double q0 = 0.5;//y(0) н.у.
	 double q1 = 0.5;//dy(0)

	//double q0 = 2.0;//y(0)
	//double q1 = 1.8;//dy(0)


	///////////////////////////////////////////////////МЕТОД ЭЙЛЕРА/////////////////////////////////////////////
	double* x = new double[n];//иксы пbо сетке h
	double* y = new double[n];
	double* dy = new double[n];
	for (int i = 0; i < n; i++) {//равномерная сетка с шагом h
		x[i] = a + i * h;
	}
	y[0] = q0;//н.у. y(0)
	dy[0] = q1;//н.у. y'(0)
	for (int i = 0; i < n - 1; i++) {
		dy[i + 1] = dy[i] + f(x[i], y[i], dy[i]) * h;
		y[i + 1] = y[i] + dy[i] * h;
	}

	double* x05 = new double[n05];//для измельчённой сетки
	double* y05 = new double[n05];
	double* dy05 = new double[n05];

	for (int i = 0; i < n05; i++) {//равномерная сетка с шагом h/2
		x05[i] = a + i * h05;
	}
	y05[0] = q0;
	dy05[0] = q1;
	for (int i = 0; i < n05 - 1; i++) {
		dy05[i + 1] = dy05[i] + f(x05[i], y05[i], dy05[i]) * h05;
		y05[i + 1] = y05[i] + dy05[i] * h05;
	}

	double* x2 = new double[n2];//для измельчённой сетки
	double* y2 = new double[n2];
	double* dy2 = new double[n2];
	for (int i = 0; i < n2; i++) {//равномерная сетка с шагом 2h
		x2[i] = a + i * h2;
	}
	y2[0] = q0;
	dy2[0] = q1;
	for (int i = 0; i < n2 - 1; i++) {
		dy2[i + 1] = dy2[i] + f(x2[i], y2[i], dy2[i]) * h2;
		y2[i + 1] = y2[i] + dy2[i] * h2;
	}
	////////////////////////////////////////////////////МЕТОД ЭЙЛЕРА С ПЕРЕСЧЁТОМ////////////////////////////////////////////////////
	double* y_sper = new double[n];
	double* dy_sper = new double[n];
	double u_i_zv = 0.0;//для пересчёта
	double y_i_zv = 0.0;//для пересчёта
	y_sper[0] = q0;
	dy_sper[0] = q1;
	for (int i = 0; i < n - 1; i++) {
		u_i_zv = dy_sper[i] + h * f(x[i], y_sper[i], dy_sper[i]);//по эйлеру ui*
		y_i_zv = y_sper[i] + h * dy_sper[i];//по эйлеру yi*
		y_sper[i + 1] = y_sper[i] + (h / 2.0) * (u_i_zv+ dy_sper[i]);
		dy_sper[i + 1] = dy_sper[i] + (h / 2.0) * (f(x[i], y_sper[i], dy_sper[i]) + f(x[i + 1], y_i_zv, u_i_zv));
	}
	 

	double* y05_sper = new double[n05];
	double* dy05_sper = new double[n05];
	y05_sper[0] = q0;
	dy05_sper[0] = q1;
	for (int i = 0; i < n05 - 1; i++) {
		u_i_zv = dy05_sper[i] + h05 * f(x05[i], y05_sper[i], dy05_sper[i]);//по эйлеру ui*
		y_i_zv = y05_sper[i] + h05 * dy05_sper[i];//по эйлеру yi*
		y05_sper[i + 1] = y05_sper[i] + (h05 / 2.0) * (u_i_zv + dy05_sper[i]);
		dy05_sper[i + 1] = dy05_sper[i] + (h05 / 2.0) * (f(x05[i], y05_sper[i], dy05_sper[i]) + f(x05[i + 1], y_i_zv, u_i_zv));
	}

	double* y2_sper = new double[n2];
	double* dy2_sper = new double[n2];
	y2_sper[0] = q0;
	dy2_sper[0] = q1;
	for (int i = 0; i < n2 - 1; i++) {
		u_i_zv = dy2_sper[i] + h2 * f(x2[i], y2_sper[i], dy2_sper[i]);//по эйлеру ui*
		y_i_zv = y2_sper[i] + h2 * dy2_sper[i];//по эйлеру yi*
		y2_sper[i + 1] = y2_sper[i] + (h2 / 2.0) * (u_i_zv + dy2_sper[i]);
		dy2_sper[i + 1] = dy2_sper[i] + (h2 / 2.0) * (f(x[i], y2_sper[i], dy2_sper[i]) + f(x[i + 1], y_i_zv, u_i_zv));
	}
	///////////////////////////////////////////////////////////////////Рунге-Кутта 2 порядка////////////////
	double* y_RK2 = new double[n];
	double* dy_RK2 = new double[n];
	y_RK2[0] = q0;
	dy_RK2[0] = q1;
	Runge_2(x, y_RK2, dy_RK2, h, n);
	double* y05_RK2 = new double[n05];
	double* dy05_RK2 = new double[n05];
	y05_RK2[0] = q0;
	dy05_RK2[0] = q1;
	Runge_2(x05, y05_RK2, dy05_RK2, h05, n05);
	double* y2_RK2 = new double[n2];
	double* dy2_RK2 = new double[n2];
	y2_RK2[0] = q0;
	dy2_RK2[0] = q1;
	Runge_2(x2, y2_RK2, dy2_RK2, h2, n2);
	////////////////////////////////////////////////////////////////////////////Рунге-Кутта 4 порядка////////////////
	double* y_RK4 = new double[n];
	double* dy_RK4 = new double[n];
	y_RK4[0] = q0;
	dy_RK4[0] = q1;
	Runge_4(x, y_RK4, dy_RK4, h, n);

	double* y05_RK4 = new double[n05];
	double* dy05_RK4 = new double[n05];
	y05_RK4[0] = q0;
	dy05_RK4[0] = q1;
	Runge_4(x05, y05_RK4, dy05_RK4, h05, n05);

	double* y2_RK4 = new double[n2];
	double* dy2_RK4 = new double[n2];
	y2_RK4[0] = q0;
	dy2_RK4[0] = q1;
	Runge_4(x2, y2_RK4, dy2_RK4, h2, n2);
	///////////////////////////////////////////////////////////////////АДАМС 3 порядка/////////////////////////////
	double* y_A = new double[n];
	double* dy_A = new double[n];
	for (int i = 0; i < n; i++) {//равномерная сетка с шагом 2h
		y_A[i] = 0;
		dy_A[i] = 0;
	}
	y_A[0] = q0;
	dy_A[0] = q1;
	int k = 3;
	Runge_4(x, y_A, dy_A, h, k);//1 и 2 находим по Р-К
	Adams(x, y_A, dy_A, h, n);//остальные по формулам 

	double* y05_A = new double[n05];
	double* dy05_A = new double[n05];
	for (int i = 0; i < n05; i++) {//равномерная сетка с шагом 2h
		y05_A[i] = 0;
		dy05_A[i] = 0;
	}
	y05_A[0] = q0;
	dy05_A[0] = q1;
	k = 3;
	Runge_4(x05, y05_A, dy05_A, h05, k);
	Adams(x05, y05_A, dy05_A, h05, n05);

	double* y2_A = new double[n2];
	double* dy2_A = new double[n2];
	for (int i = 0; i < n2; i++) {//равномерная сетка с шагом 2h
		y2_A[i] = 0;
		dy2_A[i] = 0;
	}
	y2_A[0] = q0;
	dy2_A[0] = q1;
	k = 3;
	Runge_4(x2, y2_A, dy2_A, h2, k);
	Adams(x2, y2_A, dy2_A, h2, n2);
	//////////////////////////////////////////ПОГРЕШНОСТИ////////////////////////////////
	double e1_abs_1_05 = 0, e1_abs_1_2 = 0, //  1 мтод Эйлера e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2, 2- шаг 2h
		e1_otn_1_05 = 0, e1_otn_1_2 = 0,
		e2_abs_1_05 = 0, e2_abs_1_2 = 0,
		e2_otn_1_05 = 0, e2_otn_1_2 = 0,
		einf_abs_1_05 = 0, einf_abs_1_2 = 0,
		einf_otn_1_05 = 0, einf_otn_1_2 = 0,

		e1_abs_2_05 = 0, e1_abs_2_2 = 0, //  2 мтод Эйлера e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2, 2- шаг 2h
		e1_otn_2_05 = 0, e1_otn_2_2 = 0,
		e2_abs_2_05 = 0, e2_abs_2_2 = 0,
		e2_otn_2_05 = 0, e2_otn_2_2 = 0,
		einf_abs_2_05 = 0, einf_abs_2_2 = 0,
		einf_otn_2_05 = 0, einf_otn_2_2 = 0,

		e1_abs_3_05 = 0, e1_abs_3_2 = 0, //  3 метод РК e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2, 2- шаг 2h
		e1_otn_3_05 = 0, e1_otn_3_2 = 0,
		e2_abs_3_05 = 0, e2_abs_3_2 = 0,
		e2_otn_3_05 = 0, e2_otn_3_2 = 0,
		einf_abs_3_05 = 0, einf_abs_3_2 = 0,
		einf_otn_3_05 = 0, einf_otn_3_2 = 0,

		e1_abs_4_05 = 0, e1_abs_4_2 = 0, //  4 метод РК e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2, 2- шаг 2h
		e1_otn_4_05 = 0, e1_otn_4_2 = 0,
		e2_abs_4_05 = 0, e2_abs_4_2 = 0,
		e2_otn_4_05 = 0, e2_otn_4_2 = 0,
		einf_abs_4_05 = 0, einf_abs_4_2 = 0,
		einf_otn_4_05 = 0, einf_otn_4_2 = 0,

		e1_abs_5_05 = 0, e1_abs_5_2 = 0, //  5 метод Адамса e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2, 2- шаг 2h
		e1_otn_5_05 = 0, e1_otn_5_2 = 0,
		e2_abs_5_05 = 0, e2_abs_5_2 = 0,
		e2_otn_5_05 = 0, e2_otn_5_2 = 0,
		einf_abs_5_05 = 0, einf_abs_5_2 = 0,
		einf_otn_5_05 = 0, einf_otn_5_2 = 0;

	//для dy
	double de1_abs_1_05 = 0, de1_abs_1_2 = 0, //  1 мтод Эйлера e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2, 2- шаг 2h
		de1_otn_1_05 = 0, de1_otn_1_2 = 0,
		de2_abs_1_05 = 0, de2_abs_1_2 = 0,
		de2_otn_1_05 = 0, de2_otn_1_2 = 0,
		deinf_abs_1_05 = 0, deinf_abs_1_2 = 0,
		deinf_otn_1_05 = 0, deinf_otn_1_2 = 0,

		de1_abs_2_05 = 0, de1_abs_2_2 = 0, //  2 мтод Эйлера e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2, 2- шаг 2h
		de1_otn_2_05 = 0, de1_otn_2_2 = 0,
		de2_abs_2_05 = 0, de2_abs_2_2 = 0,
		de2_otn_2_05 = 0, de2_otn_2_2 = 0,
		deinf_abs_2_05 = 0, deinf_abs_2_2 = 0,
		deinf_otn_2_05 = 0, deinf_otn_2_2 = 0,

		de1_abs_3_05 = 0, de1_abs_3_2 = 0, //  3 метод РК e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2, 2- шаг 2h
		de1_otn_3_05 = 0, de1_otn_3_2 = 0,
		de2_abs_3_05 = 0, de2_abs_3_2 = 0,
		de2_otn_3_05 = 0, de2_otn_3_2 = 0,
		deinf_abs_3_05 = 0, deinf_abs_3_2 = 0,
		deinf_otn_3_05 = 0, deinf_otn_3_2 = 0,

		de1_abs_4_05 = 0, de1_abs_4_2 = 0, //  4 метод РК e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2, 2- шаг 2h
		de1_otn_4_05 = 0, de1_otn_4_2 = 0,
		de2_abs_4_05 = 0, de2_abs_4_2 = 0,
		de2_otn_4_05 = 0, de2_otn_4_2 = 0,
		deinf_abs_4_05 = 0, deinf_abs_4_2 = 0,
		deinf_otn_4_05 = 0, deinf_otn_4_2 = 0,

		de1_abs_5_05 = 0, de1_abs_5_2 = 0, //  5 метод Адамса e1-||.||1  e2-||.||2 einf-||.||inf, 05-шаг h/2, 2- шаг 2h
		de1_otn_5_05 = 0, de1_otn_5_2 = 0,
		de2_abs_5_05 = 0, de2_abs_5_2 = 0,
		de2_otn_5_05 = 0, de2_otn_5_2 = 0,
		deinf_abs_5_05 = 0, deinf_abs_5_2 = 0,
		deinf_otn_5_05 = 0, deinf_otn_5_2 = 0;
	//для y
	//между h и h/2	
	pogreshnost1(e1_abs_1_05, e1_otn_1_05, e2_abs_1_05, e2_otn_1_05, einf_abs_1_05, einf_otn_1_05, n, y, y05);
	pogreshnost1(e1_abs_2_05, e1_otn_2_05, e2_abs_2_05, e2_otn_2_05, einf_abs_2_05, einf_otn_2_05, n, y_sper, y05_sper);
	pogreshnost1(e1_abs_3_05, e1_otn_3_05, e2_abs_3_05, e2_otn_3_05, einf_abs_3_05, einf_otn_3_05, n, y_RK2, y05_RK2);
	pogreshnost1(e1_abs_4_05, e1_otn_4_05, e2_abs_4_05, e2_otn_4_05, einf_abs_4_05, einf_otn_4_05, n, y_RK4, y05_RK4);
	pogreshnost1(e1_abs_5_05, e1_otn_5_05, e2_abs_5_05, e2_otn_5_05, einf_abs_5_05, einf_otn_5_05, n, y_A, y05_A);
	//между h и 2h	
	pogreshnost2(e1_abs_1_2, e1_otn_1_2, e2_abs_1_2, e2_otn_1_2, einf_abs_1_2, einf_otn_1_2, n2, y, y2);
	pogreshnost2(e1_abs_2_2, e1_otn_2_2, e2_abs_2_2, e2_otn_2_2, einf_abs_2_2, einf_otn_2_2, n2, y_sper, y2_sper);
	pogreshnost2(e1_abs_3_2, e1_otn_3_2, e2_abs_3_2, e2_otn_3_2, einf_abs_3_2, einf_otn_3_2, n2, y_RK2, y2_RK2);
	pogreshnost2(e1_abs_4_2, e1_otn_4_2, e2_abs_4_2, e2_otn_4_2, einf_abs_4_2, einf_otn_4_2, n2, y_RK4, y2_RK4);
	pogreshnost2(e1_abs_5_2, e1_otn_5_2, e2_abs_5_2, e2_otn_5_2, einf_abs_5_2, einf_otn_5_2, n2, y_A, y2_A);

	//для dy
	//между h и h/2	
	pogreshnost1(de1_abs_1_05, de1_otn_1_05, de2_abs_1_05, de2_otn_1_05, deinf_abs_1_05, deinf_otn_1_05, n, dy, dy05);
	pogreshnost1(de1_abs_2_05, de1_otn_2_05, de2_abs_2_05, de2_otn_2_05, deinf_abs_2_05, deinf_otn_2_05, n, dy_sper, dy05_sper);
	pogreshnost1(de1_abs_3_05, de1_otn_3_05, de2_abs_3_05, de2_otn_3_05, deinf_abs_3_05, deinf_otn_3_05, n, dy_RK2, dy05_RK2);
	pogreshnost1(de1_abs_4_05, de1_otn_4_05, de2_abs_4_05, de2_otn_4_05, deinf_abs_4_05, deinf_otn_4_05, n, dy_RK4, dy05_RK4);
	pogreshnost1(de1_abs_5_05, de1_otn_5_05, de2_abs_5_05, de2_otn_5_05, deinf_abs_5_05, deinf_otn_5_05, n, dy_A, dy05_A);
	//между h и 2h	
	pogreshnost2(de1_abs_1_2, de1_otn_1_2, de2_abs_1_2, de2_otn_1_2, deinf_abs_1_2, deinf_otn_1_2, n2, dy, dy2);
	pogreshnost2(de1_abs_2_2, de1_otn_2_2, de2_abs_2_2, de2_otn_2_2, deinf_abs_2_2, deinf_otn_2_2, n2, dy_sper, dy2_sper);
	pogreshnost2(de1_abs_3_2, de1_otn_3_2, de2_abs_3_2, de2_otn_3_2, deinf_abs_3_2, deinf_otn_3_2, n2, dy_RK2, dy2_RK2);
	pogreshnost2(de1_abs_4_2, de1_otn_4_2, de2_abs_4_2, de2_otn_4_2, deinf_abs_4_2, deinf_otn_4_2, n2, dy_RK4, dy2_RK4);
	pogreshnost2(de1_abs_5_2, de1_otn_5_2, de2_abs_5_2, de2_otn_5_2, deinf_abs_5_2, deinf_otn_5_2, n2, dy_A, dy2_A);
	

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
		fout3 << dy[i] << "\n";
	}
	fout3.close();

	ofstream fout4;
	fout4.open("file4.txt");
	for (int i = 0; i < n; i++) {
		fout4 << fixed;
		fout4.precision(5);
		fout4.setf(ios::right);
		fout4.width(10);
		fout4 << y_sper[i] << "\n";
	}
	fout4.close();

	ofstream fout5;
	fout5.open("file5.txt");
	for (int i = 0; i < n; i++) {
		fout5 << fixed;
		fout5.precision(5);
		fout5.setf(ios::right);
		fout5.width(10);
		fout5 << dy_sper[i] << "\n";
	}
	fout5.close();

	ofstream fout6;
	fout6.open("file6.txt");
	for (int i = 0; i < n; i++) {
		fout6 << fixed;
		fout6.precision(5);
		fout6.setf(ios::right);
		fout6.width(10);
		fout6 << y_RK2[i] << "\n";
	}
	fout6.close();

	ofstream fout7;
	fout7.open("file7.txt");
	for (int i = 0; i < n; i++) {
		fout7 << fixed;
		fout7.precision(5);
		fout7.setf(ios::right);
		fout7.width(10);
		fout7 << dy_RK2[i] << "\n";
	}
	fout7.close();

	ofstream fout8;
	fout8.open("file8.txt");
	for (int i = 0; i < n; i++) {
		fout8 << fixed;
		fout8.precision(5);
		fout8.setf(ios::right);
		fout8.width(10);
		fout8 << y_RK4[i] << "\n";
	}
	fout8.close();

	ofstream fout9;
	fout9.open("file9.txt");
	for (int i = 0; i < n; i++) {
		fout9 << fixed;
		fout9.precision(5);
		fout9.setf(ios::right);
		fout9.width(10);
		fout9 << dy_RK4[i] << "\n";
	}
	fout9.close();

	ofstream fout10;
	fout10.open("file10.txt");
	for (int i = 0; i < n; i++) {
		fout10 << fixed;
		fout10.precision(5);
		fout10.setf(ios::right);
		fout10.width(10);
		fout10 << y_A[i] << "\n";
	}
	fout10.close();

	ofstream fout11;
	fout11.open("file11.txt");
	for (int i = 0; i < n; i++) {
		fout11 << fixed;
		fout11.precision(5);
		fout11.setf(ios::right);
		fout11.width(10);
		fout11 << dy_A[i] << "\n";
	}
	fout11.close();

	ofstream fout12;
	fout12.open("file12.txt");
	fout12 << "tabl for y" << endl;
	fout12.setf(ios::scientific);

	fout12 << endl;
	fout12 << setw(40) << left << " " << setw(15) << left << " " << setw(15) << left << " h/2" << setw(15) << left << " " << setw(15) << left << " " << setw(15) << left << " " << setw(15) << left << " " << setw(15) << left << " " << setw(15) << left << "2h " << setw(15) << left << " " << setw(15) << left << " " << setw(15) << left << " " << endl;
	fout12 << endl;
	fout12 << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << " " << setw(30) << left << "         ||.||1" << setw(30) << left << "         ||.||2" << setw(30) << left << "         ||.||inf " << setw(30) << left << "         ||.||1 " << setw(30) << left << "         ||.||2" << setw(30) << left << "       ||.||inf" << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left<< "metod" << setw(15) << left << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs " << setw(15) << left << "otn" << setw(15) << left << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs" << setw(15) << left << "otn" << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "Euler" << setw(15) << left << e1_abs_1_05 << setw(15) << left << e1_otn_1_05 << setw(15) << left << e2_abs_1_05 << setw(15) << left << e2_otn_1_05 << setw(15) << left << einf_abs_1_05 << setw(15) << left << einf_otn_1_05 << setw(15) << left << e1_abs_1_2 << setw(15) << left << e1_otn_1_2 << setw(15) << left << e2_abs_1_2 << setw(15) << left << e2_otn_1_2 << setw(15) << left << einf_abs_1_2 << setw(15) << left << einf_otn_1_2 << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "Euler s peres" << setw(15) << left << e1_abs_2_05 << setw(15) << left << e1_otn_2_05 << setw(15) << left << e2_abs_2_05 << setw(15) << left << e2_otn_2_05 << setw(15) << left << einf_abs_2_05 << setw(15) << left << einf_otn_2_05 << setw(15) << left << e1_abs_2_2 << setw(15) << left << e1_otn_2_2 << setw(15) << left << e2_abs_2_2 << setw(15) << left << e2_otn_2_2 << setw(15) << left << einf_abs_2_2 << setw(15) << left << einf_otn_2_2 << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "Runge-Kutta 2" << setw(15) << left << e1_abs_3_05 << setw(15) << left << e1_otn_3_05 << setw(15) << left << e2_abs_3_05 << setw(15) << left << e2_otn_3_05 << setw(15) << left << einf_abs_3_05 << setw(15) << left << einf_otn_3_05 << setw(15) << left << e1_abs_3_2 << setw(15) << left << e1_otn_3_2 << setw(15) << left << e2_abs_3_2 << setw(15) << left << e2_otn_3_2 << setw(15) << left << einf_abs_3_2 << setw(15) << left << einf_otn_3_2 << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "Runge-Kutta 4" << setw(15) << left << e1_abs_4_05 << setw(15) << left << e1_otn_4_05 << setw(15) << left << e2_abs_4_05 << setw(15) << left << e2_otn_4_05 << setw(15) << left << einf_abs_4_05 << setw(15) << left << einf_otn_4_05 << setw(15) << left << e1_abs_4_2 << setw(15) << left << e1_otn_4_2 << setw(15) << left << e2_abs_4_2 << setw(15) << left << e2_otn_4_2 << setw(15) << left << einf_abs_4_2 << setw(15) << left << einf_otn_4_2 << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "Adams 3" << setw(15) << left << e1_abs_5_05 << setw(15) << left << e1_otn_5_05 << setw(15) << left << e2_abs_5_05 << setw(15) << left << e2_otn_5_05 << setw(15) << left << einf_abs_5_05 << setw(15) << left << einf_otn_5_05 << setw(15) << left << e1_abs_5_2 << setw(15) << left << e1_otn_5_2 << setw(15) << left << e2_abs_5_2 << setw(15) << left << e2_otn_5_2 << setw(15) << left << einf_abs_5_2 << setw(15) << left << einf_otn_5_2 << endl;
	
	fout12.close();

	ofstream fout13;
	fout13.open("file13.txt");
	fout13 << "tabl for dy" << endl;
	fout13.setf(ios::scientific);

	fout13 << endl;
	fout13 << setw(40) << left << " " << setw(15) << left << " " << setw(15) << left << " h/2" << setw(15) << left << " " << setw(15) << left << " " << setw(15) << left << " " << setw(15) << left << " " << setw(15) << left << " " << setw(15) << left << "2h " << setw(15) << left << " " << setw(15) << left << " " << setw(15) << left << " " << endl;
	fout13 << endl;
	fout13 << endl;
	fout13 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout13 << setw(15) << left << " " << setw(30) << left << "         ||.||1" << setw(30) << left << "         ||.||2" << setw(30) << left << "         ||.||inf " << setw(30) << left << "         ||.||1 " << setw(30) << left << "         ||.||2" << setw(30) << left << "       ||.||inf" << endl;
	fout13 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout13 << setw(15) << left << "metod" << setw(15) << left << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs " << setw(15) << left << "otn" << setw(15) << left << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs" << setw(15) << left << "otn" << endl;
	fout13 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout13 << setw(15) << left << "Euler" << setw(15) << left << de1_abs_1_05 << setw(15) << left << de1_otn_1_05 << setw(15) << left << de2_abs_1_05 << setw(15) << left << de2_otn_1_05 << setw(15) << left << deinf_abs_1_05 << setw(15) << left << deinf_otn_1_05 << setw(15) << left << de1_abs_1_2 << setw(15) << left << de1_otn_1_2 << setw(15) << left << de2_abs_1_2 << setw(15) << left << de2_otn_1_2 << setw(15) << left << deinf_abs_1_2 << setw(15) << left << deinf_otn_1_2 << endl;
	fout13 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout13 << setw(15) << left << "Euler s peres" << setw(15) << left << de1_abs_2_05 << setw(15) << left << de1_otn_2_05 << setw(15) << left << de2_abs_2_05 << setw(15) << left << de2_otn_2_05 << setw(15) << left << deinf_abs_2_05 << setw(15) << left << deinf_otn_2_05 << setw(15) << left << de1_abs_2_2 << setw(15) << left << de1_otn_2_2 << setw(15) << left << de2_abs_2_2 << setw(15) << left << de2_otn_2_2 << setw(15) << left << deinf_abs_2_2 << setw(15) << left << deinf_otn_2_2 << endl;
	fout13 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout13 << setw(15) << left << "Runge-Kutta 2" << setw(15) << left << de1_abs_3_05 << setw(15) << left << de1_otn_3_05 << setw(15) << left << de2_abs_3_05 << setw(15) << left << de2_otn_3_05 << setw(15) << left << deinf_abs_3_05 << setw(15) << left << deinf_otn_3_05 << setw(15) << left << de1_abs_3_2 << setw(15) << left << de1_otn_3_2 << setw(15) << left << de2_abs_3_2 << setw(15) << left << de2_otn_3_2 << setw(15) << left << deinf_abs_3_2 << setw(15) << left << deinf_otn_3_2 << endl;
	fout13 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout13 << setw(15) << left << "Runge-Kutta 4" << setw(15) << left << de1_abs_4_05 << setw(15) << left << de1_otn_4_05 << setw(15) << left << de2_abs_4_05 << setw(15) << left << de2_otn_4_05 << setw(15) << left << deinf_abs_4_05 << setw(15) << left << deinf_otn_4_05 << setw(15) << left << de1_abs_4_2 << setw(15) << left << de1_otn_4_2 << setw(15) << left << de2_abs_4_2 << setw(15) << left << de2_otn_4_2 << setw(15) << left << deinf_abs_4_2 << setw(15) << left << deinf_otn_4_2 << endl;
	fout13 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout13 << setw(15) << left << "Adams 3" << setw(15) << left << de1_abs_5_05 << setw(15) << left << de1_otn_5_05 << setw(15) << left << de2_abs_5_05 << setw(15) << left << de2_otn_5_05 << setw(15) << left << deinf_abs_5_05 << setw(15) << left << deinf_otn_5_05 << setw(15) << left << de1_abs_5_2 << setw(15) << left << de1_otn_5_2 << setw(15) << left << de2_abs_5_2 << setw(15) << left << de2_otn_5_2 << setw(15) << left << deinf_abs_5_2 << setw(15) << left << deinf_otn_5_2 << endl;

	fout13.close();







	delete[] x;
	delete[] y;
	delete[] dy;
	delete[] x05;
	delete[] y05;
	delete[] dy05;
	delete[] x2;
	delete[] y2;
	delete[] dy2;
	delete[] y_sper;
	delete[] dy_sper;
	delete[] y05_sper;
	delete[] dy05_sper;
	delete[] y2_sper;
	delete[] dy2_sper;
	delete[] y_RK2;
	delete[] dy_RK2;
	delete[] y05_RK2;
	delete[] dy05_RK2;
	delete[] y2_RK2;
	delete[] dy2_RK2;
	delete[] y_RK4;
	delete[] dy_RK4;
	delete[] y05_RK4;
	delete[] dy05_RK4;
	delete[] y2_RK4;
	delete[] dy2_RK4;
	delete[] y_A;
	delete[] dy_A;
	delete[] y05_A;
	delete[] dy05_A;
	delete[] y2_A;
	delete[] dy2_A;

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
void pogreshnost2(double& e1_h2_abs, double& e1_h2_othn, double& e2_h2_abs, double& e2_h2_othn, double& einf_h2_abs, double& einf_h2_othn, int n, double* y, double* y2)
{
	for (int i = 0; i < n; i++)
	{
		e1_h2_abs += fabs(y[i*2] - y2[i]);
		e1_h2_othn += fabs(y[i*2]);
		e2_h2_abs += pow(fabs(y[i*2] - y2[i]), 2);
		e2_h2_othn += pow(fabs(y[i*2]), 2);
		einf_h2_abs = max(einf_h2_abs, fabs(y[i*2] - y2[i]));
		einf_h2_othn = max(einf_h2_othn, fabs(y[i*2]));
	}

	e1_h2_othn = e1_h2_abs / e1_h2_othn;
	e2_h2_abs = sqrt(e2_h2_abs);
	e2_h2_othn = sqrt(e2_h2_othn);
	e2_h2_othn = e2_h2_abs / e2_h2_othn;
	einf_h2_othn = einf_h2_abs / einf_h2_othn;
}