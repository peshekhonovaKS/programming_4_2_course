#include <fstream> 
#include <iostream>
#include <time.h>
#include <iomanip>//для вывода красивого в файлики
using namespace std;

void Newtone(double& x, double& y, double e, int& j) {
	double j11, j12, j21, j22; //элементы Якобиана
	double f1, f2;             //правые части, которые будем брать с "-" для решения J(x(k))*delta_x(k)=-F(x(k)) delta_x(k)=x(k+1)-x(k)
	double delta_x, delta_y;
	f1 = -(x - pow(y, 3));
	f2 = -(log(x + 2) - y);
	j11 = 1;                   //считаем ручками якобиан
	j12 = -3 * pow(y, 2);
	j21 = 1 / (x + 2);
	j22 = -1;
	j = 1;                     //кол-во итераций
	delta_x = (f1 * j22 - f2 * j12) / (j11 * j22 - j12 * j21);   //выразила дельту из J(x(k))*delta_x(k)=-F(x(k))
	delta_y = (f2 * j11 - f1 * j21) / (j11 * j22 - j12 * j21);
	while((sqrt(pow(delta_x,2) + pow(delta_y,2)) > e) && (j<1000)){
		f1 = -(x - pow(y, 3));
		f2 = -(log((x + 2)) - y);
		j11 = 1;
		j12 = -3 * pow(y, 2);
		j21 = 1 / (x + 2);
		j22 = -1;
		delta_x = (f1 * j22 - f2 * j12) / (j11 * j22 - j12 * j21);
		delta_y = (f2 * j11 - f1 * j21) / (j11 * j22 - j12 * j21);
		x += delta_x;
		y += delta_y;
		j += 1;
	}
}
void Newtone2(double& x, double& y, double e, int& j, double j11, double j12, double j21, double j22) { //матрица якоби постоянная,считаю в main, отправляем ещё ее элемены
	double f1, f2;
	double delta_x, delta_y;
	f1 = -(x - pow(y, 3));
	f2 = -(log(abs(x + 2)) - y);
	j = 1;
	delta_x = (f1 * j22 - f2 * j12) / (j11 * j22 - j12 * j21);
	delta_y = (f2 * j11 - f1 * j21) / (j11 * j22 - j12 * j21);
	while ((sqrt(pow(delta_x, 2) + pow(delta_y, 2)) > e) && (j < 1000)) {
		f1 = -(x - pow(y, 3));
		f2 = -(log(abs(x + 2)) - y);
		delta_x = (f1 * j22 - f2 * j12) / (j11 * j22 - j12 * j21);
		delta_y = (f2 * j11 - f1 * j21) / (j11 * j22 - j12 * j21);
		x += delta_x;
		y += delta_y;
		j += 1;
	}
}
//////для дискр метода с разностной аппроксимацией производных//////
double f1(double x, double y) {//1ое уравнение
	return (x - pow(y, 3));
}
double f2(double x, double y) {//2ое уравнение
	return (log(abs(x + 2)) - y);
}
void J3(double x, double y, double& j11, double& j12, double& j21, double& j22, double h) {//якобиан с помощью разностной схемы
	j11 = (f1(x + h, y) - f1(x - h, y)) / (2 * h);
	j12 = (f1(x, y + h) - f1(x, y - h)) / (2 * h);
	j21 = (f2(x + h, y) - f2(x - h, y)) / (2 * h);
	j22 = (f2(x, y + h) - f2(x, y - h)) / (2 * h);;
}
void Newtone3(double& x, double& y, double e, int& j, double h) {
	double j11, j12, j21, j22;
	double f1, f2;
	double delta_x, delta_y;
	f1 = -(x - pow(y, 3));
	f2 = -(log(abs(x + 2)) - y);
	j11 = 1;   //на первой считаем ручками потом  с помощью разностной апроксимации считаем якобиан
	j12 = -3 * pow(y, 2);
	j21 = 1 / (x + 2);
	j22 = -1;
	j = 1;
	delta_x = (f1 * j22 - f2 * j12) / (j11 * j22 - j12 * j21);
	delta_y = (f2 * j11 - f1 * j21) / (j11 * j22 - j12 * j21);
	while ((sqrt(pow(delta_x, 2) + pow(delta_y, 2)) > e) && (j < 1000)) {
		f1 = -(x - pow(y, 3));
		f2 = -(log(abs(x + 2)) - y);
		J3(x, y, j11, j12, j21, j22, h);
		delta_x = (f1 * j22 - f2 * j12) / (j11 * j22 - j12 * j21);
		delta_y = (f2 * j11 - f1 * j21) / (j11 * j22 - j12 * j21);
		x += delta_x;
		y += delta_y;
		j += 1;
	}
}
int main() {
	double e = 1e-7;//точность вычислений
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!МПИ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//начальные приближения
	double pribl_x1 = 97.0, pribl_y1 = 4.6;//для 1 корня
	double pribl_x2 = -1.7, pribl_y2 = -1.19;//для 2 корня
	double delta_x1 = 1, delta_y1 = 1;//x^(k+1)-x^(k) и y^(k+1)-y^(k)
	double delta_x2 = 1, delta_y2 = 1;
	int kol_iter1_MPI = 0;//счётчик итераций
	int kol_iter2_MPI = 0;
	while((max(abs(delta_x1), abs(delta_y1)) > e)&&(kol_iter1_MPI<1000)) {
		double a = pribl_x1;//запоминаем на предыдущем шаге
		double b = pribl_y1;
		pribl_x1 = pow(b, 3);//вычисояем на новом 
		pribl_y1 = log(a + 2);
		delta_x1 = pribl_x1 - a;//разница между следующей и текущей итерацией
		delta_y1 = pribl_y1 - b;      
		kol_iter1_MPI += 1;
	}
	while((max(abs(delta_x2), abs(delta_y2)) > e) && (kol_iter2_MPI < 1000)) {
		double a = pribl_x2;
		double b = pribl_y2;
		pribl_x2 = exp(b) - 2;
		pribl_y2 = cbrt(a);
		delta_x2 = pribl_x2 - a;        //разница между следующей и текущей итерацией
		delta_y2 = pribl_y2 - b;
		kol_iter2_MPI += 1;
	}
	//cout << setprecision(10) << pribl_x1 << endl << pribl_y1;
	//cout << endl << kol_iter1_MPI << endl;
	//cout << pribl_x2 << endl << pribl_y2;
	//cout << endl << kol_iter2_MPI << endl;


	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!метод Ньютона!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double priblN_x1 = 97.0, priblN_y1 = 4.6;//нач приближение
	double priblN_x2 = -1.7, priblN_y2 = -1.19;
	int kol_iter1N = 0;//интерации
	int kol_iter2N = 0;
	Newtone(priblN_x1, priblN_y1, e, kol_iter1N);
	Newtone(priblN_x2, priblN_y2, e, kol_iter2N);
	//cout << "Newtoun" << endl;
	//cout << setprecision(10) << priblN_x1 << endl << priblN_y1 << endl << kol_iter1N << endl;
	//cout << priblN_x2 << endl << priblN_y2 << endl << kol_iter2N << endl;

	//Модифицированный метод Ньютона
	double priblN2_x1 = 97.0, priblN2_y1 = 4.6;
	double priblN2_x2 = -1.7, priblN2_y2 = -1.19;
	int kol_iter1N2 = 0;
	int kol_iter2N2 = 0;
	double a11, a12, a21, a22;
	a11 = 1;//якобиан для 1 корня
	a12 = -3 * pow(priblN2_y1, 2);
	a21 = 1 / (priblN2_x1 + 2);
	a22 = -1;

	double u11, u12, u21, u22;
	u11 = 1;//якобиан для 2 корня
	u12 = -3 * pow(priblN2_y2, 2);
	u21 = 1 / (priblN2_x2 + 2);
	u22 = -1;

	Newtone2(priblN2_x1, priblN2_y1, e, kol_iter1N2, a11, a12, a21, a22);
	Newtone2(priblN2_x2, priblN2_y2, e, kol_iter2N2, u11, u12, u21, u22);
	//cout << "Newtoun2" << endl;
	//cout << setprecision(10) << priblN2_x1 << endl << priblN2_y1 << endl << kol_iter1N2 << endl;
	//cout << priblN2_x2 << endl << priblN2_y2 << endl << kol_iter2N2 << endl;


	//метод Ньютона c матрицей Якоби посчитанной разностной схемой
	double priblN3_x1 = 97.0, priblN3_y1 = 4.6;
	double priblN3_x2 = -1.7, priblN3_y2 = -1.19;
	int kol_iter1N3 = 0;
	int kol_iter2N3 = 0;
	double h = 1e-4;
	Newtone3(priblN3_x1, priblN3_y1, e, kol_iter1N3, h);
	Newtone3(priblN3_x2, priblN3_y2, e, kol_iter2N3, h);
	//cout << "Newtoun3" << endl;
	//cout << setprecision(10) << priblN3_x1 << endl << priblN3_y1 << endl << kol_iter1N3 << endl;

	ofstream fout8;
	fout8.open("file1.txt");
	fout8 << setw(30) << left << "metod" << setw(30) << left << "x" << setw(30) << left << "y" << setw(30) << left << "kolvo_iter" << endl;
	fout8 << setw(30) << left << "MPI" << setw(30) << endl;
	fout8 << setw(30) << left << " " << setw(30) << left << setprecision(10) << pribl_x1 << setw(30) << left << pribl_y1 << setw(30) << left << kol_iter1_MPI << endl;
	fout8 << setw(30) << left << " " << setw(30) << left << setprecision(10) << pribl_x2 << setw(30) << left << pribl_y2 << setw(30) << left << kol_iter2_MPI << endl;
	fout8 << "---------------------------------------------------------------------------------------------------------" << endl;
	fout8 << setw(30) << left << "Newtoun_1" << setw(30) << endl;
	fout8 << setw(30) << left << " " << setw(30) << left << setprecision(10) << priblN_x1 << setw(30) << left << priblN_y1 << setw(30) << left << kol_iter1N << endl;
	fout8 << setw(30) << left << " " << setw(30) << left << setprecision(10) << priblN_x2 << setw(30) << left << priblN_y2 << setw(30) << left << kol_iter2N << endl;
	fout8 << "---------------------------------------------------------------------------------------------------------" << endl;
	fout8 << setw(30) << left << "Newtoun_2" << setw(30) << endl;
	fout8 << setw(30) << left << " " << setw(30) << left << setprecision(10) << priblN2_x1 << setw(30) << left << priblN2_y1 << setw(30) << left << kol_iter1N2 << endl;
	fout8 << setw(30) << left << " " << setw(30) << left << setprecision(10) << priblN2_x2 << setw(30) << left << priblN2_y2 << setw(30) << left << kol_iter2N2 << endl;
	fout8 << "---------------------------------------------------------------------------------------------------------" << endl;
	fout8 << setw(30) << left << "Newtoun_3" << setw(30) << endl;
	fout8 << setw(30) << left << " " << setw(30) << left << setprecision(10) << priblN3_x1 << setw(30) << left << priblN3_y1 << setw(30) << left << kol_iter1N3 << endl;
	fout8 << setw(30) << left << " " << setw(30) << left << setprecision(10) << priblN3_x2 << setw(30) << left << priblN3_y2 << setw(30) << left << kol_iter2N3 << endl;
	fout8 << "---------------------------------------------------------------------------------------------------------" << endl;
	
	fout8.close();



	std::system("python PythonApplication1.py");
	return 0;
}