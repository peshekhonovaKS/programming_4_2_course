#include <fstream> 
#include <iostream>
#include <time.h>
#include <iomanip>//��� ������ ��������� � �������
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;

void pogreshnost1(double& e1_abs, double& e1_othn, double& e2_abs, double& e2_othn, double& einf_abs, double& einf_othn, int M, double* x_resh, double* x);
void proizv(double** A, double* Ax, double* x, int M);

double func(double x) {
	//return  x+sin(pow(x, 2)) ;
	return sin(x);
}
double fi(int i, double x, double* metr, int N, int k) {//fi_i(x), metr - �����, �� �������� �������, N - ���������� ����� � ���������, k - ��������, � ������� ��������� 
	double znam = 1.0;
	double f = 1.0;
	for (int j = 0; j < N; j++)
		if (i != j + (N - 1) * k) {
			f *= x - metr[j + (N - 1) * k];
			znam *= metr[i] - metr[j + (N - 1) * k];
		}
	f = f / znam;
	return f;
}

int main() {
	srand(time(NULL));
	double a = -3.14;
	double b = 3.14;
	const int K = 3; //���-�� �������������
	const int N = 3; //���-�� ����� �� �������� ��������
	int L = 15;//���-�� ��������� ����� �� ����� ������������
	const int M = (N - 1) * K + 1;//����� ����� ����� �����
	double* arrxM = new double[M];//������ ���� ����� ����� �
	double* arrxL = new double[L * K];//������ ���� ��������� ����� �� ��� ab
	double* arryL = new double[L * K];//y-�� �� ��������� ������

	double h = (b - a) / (M - 1); //����������� ���
	for (int i = 0; i < M; i++)
		arrxM[i] = a + i * (b - a) / (M - 1);//      (b - a) / (M - 1)=h-����������� ���

	int l = 0;
	for (int k = 0; k < K; k++) { //�������� �� �������������
		double k1 = arrxM[k * (N - 1)]; //����,��������������� ������ ������������:����� ���������*���-�� ����� � ���������
		double k2 = arrxM[(k + 1) * (N - 1)];//����� ������������
		//cout << " k1 = " << k1 << " k2 = " << k2 << endl;
		for (; l < L * (k + 1); l++) { //����� �� ����� L c�������� ����� �� ������������, ��������� L*(k+1)=���-�� ����� �� ����� ���������
			//cout << a<<" "<<b<< endl;
			arrxL[l] = k1 + (k2 - k1) * rand() / (float)RAND_MAX; //����� ����� �� ������������ [k1k2]
			//cout << arrxL[l]<<"  ";
		}
	}

	double t;
	for (int i = 1; i < K * L; i++) {
		for (int j = 0; j < K * L - 1; j++)
			if (arrxL[j] > arrxL[j + 1]) {
				t = arrxL[j];
				arrxL[j] = arrxL[j + 1];
				arrxL[j + 1] = t;
			}

	}//��������� ��������� ����� 

	//�-�� ����������� �� ��������� ������ L, ����� �� ������� [ab] �� K*L
	for (int i = 0; i < K * L; i++) {
		arryL[i] = func(arrxL[i]);
	}

	//!!!!!!!!!!!������� ���������!!!!!!!!!!!

	//c����� ������� �
	double** arr = new double* [M];//���-�� �����=����� �� �-��=����� ����� �� �������=M ��� �� ����� ��������� �� ���������
	for (int i = 0; i < M; i++) {
		arr[i] = new double[N];//���-�� �������� �������=����� ����� �� ������������
		for (int j = 0; j < N; j++) {
			arr[i][j] = 0;
		}
	}

	MatrixXd arr1 = MatrixXd::Zero(M, M);//������ �������
	VectorXd arrright = VectorXd::Zero(M);//������ �����
	VectorXd � = VectorXd::Zero(M);//������ �1...�M

	for (int k = 0; k < K; k++) { //�� ����������
		for (int i = k * (N - 1); i <= (k + 1) * (N - 1); i++) {//i ��� ����� ������� ������� ��������� ���������� ��������� ����� ��������� �� ���-�� ����� � ��
			for (int j = i; j <= (k + 1) * (N - 1); j++) {//j-������� ������� �������
				for (int l = k * L; l < (k + 1) * L; l++) {
					arr[i][j - i] += fi(i, arrxL[l], arrxM, N, k) * fi(j, arrxL[l], arrxM, N, k);
					//cout << fi(i, arrxL[l], arrxM, N, k) << endl;
				}
			}

			for (int l = k * L; l < (k + 1) * L; l++)
				arrright(i) += arryL[l] * fi(i, arrxL[l], arrxM, N, k);
		}
	}

	//!!!!!!!!!!!������� ������!!!!!!!!!!!
	for (int j = 0; j < N; j++) {//j-������� i-������ ���� ������� 
		for (int i = 0; i < M - j; i++) {//�� ������ ���� ��������� �� 1 ��-� ������
			arr1(i, i + j) = arr[i][j];
			arr1(i + j, i) = arr[i][j];
		}
	}
	VectorXd C = arr1.lu().solve(arrright);
	//����� ������� �
	cout << "reshenie C" << endl;
	for (int i = 0; i < M; i++) {
		cout << C(i) << " ";
	}
	cout << endl;
	double x[M]; //����� ������� ���������� � ������ x
	//for (int i = 0; i < M; i++) {
		//x[i] = 0.0;
	//}
	for (int i = 0; i < M; i++) {
		x[i] = C(i);
	}
	
///////////////////����� ������� ������/////////////////////
	///������� �*� � ������� ��������� ���������
	double** arr2 = new double* [M];//���-�� �����=����� �� �-��=����� ����� �� �������=M ��� �� ����� ��������� �� ���������
	for (int i = 0; i < M; i++) {
		arr2[i] = new double[M];//���-�� �������� �������=����� ����� �� ������������
		for (int j = 0; j < M; j++) {
			arr2[i][j] = 0;
		}
	}
	for (int j = 0; j < N; j++) {//j-������� i-������ ���� ������� 
		for (int i = 0; i < M - j; i++) {//�� ������ ���� ��������� �� 1 ��-� ������
			arr2[i][i + j] = arr[i][j];
			arr2[i + j][i] = arr[i][j];
		}
	}
	
	
///////////////����� ������/////////////////////////
	double A1[M][M];//����� ������
	double x1[M];//������� ������ ������
	double B1[M];//������ ����� � ������ ������
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			A1[i][j] = arr2[i][j];
		}
	}
	for (int i = 0; i < M; i++) {
		x1[i] = 0.0;
	}
	for (int i = 0; i < M; i++) {
		B1[i] = arrright(i);
	}

	double d=0.0;//��� ���������
	for (int k = 0; k < M; k++){// ������ ���
		for (int j = k + 1; j < M; j++){
			d = A1[j][k] / A1[k][k]; // ������� (�������� �� ������)
			for (int i = k; i < M; i++){
			A1[j][i] = A1[j][i] - d * A1[k][i];
			}
			B1[j] = B1[j] - d * B1[k];
		}
	}
	double s = 0.0;//��� ����� ������������
	for (int k = M-1; k >= 0; k--){ // �������� ���
		d = 0;
		for (int j = k + 1; j < M; j++){
			s = A1[k][j] * x1[j]; 
			d = d + s; 
		}
		x1[k] = (B1[k] - d) / A1[k][k]; 
	}
	cout << endl;
	cout << "reshenie gauss" << endl;
	for (int i = 0; i < M; i++) {
		cout << x1[i] << " ";
	}
/////////////////////////////////LU ����������//////////////////////////////////
	double A2[M][M];//LU
	double x2[M];//������� LU
	double B2[M];//������ ����� � ������ LU
	double y2[M];//��� ��������� ���� LU
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			A2[i][j] = arr2[i][j];
		}
	}
	for (int i = 0; i < M; i++) {
		x2[i] = 0.0;
		y2[i] = 0.0;
	}
	for (int i = 0; i < M; i++) {
		B2[i] = arrright(i);
	}
	double LL[M][M];//������� L � LU ���������� -����������������
	double U[M][M];//������� U � LU ���������� -�����������������
	for (int i = 0; i < M; i++){
		for (int j = 0; j < M; j++){
			LL[i][j] = 0.0;
			U[i][j] = 0.0;
			if (i == j) {
				U[i][j] = 1;///�� ��������� U ����� 1, ��� ��������� ��������� ������ 
			}
		}
	}
	for (int i = 0; i < M; i++){
		LL[i][0] = A2[i][0];//������ �������
		U[0][i] = A2[0][i] / LL[0][0];//1 �������
	}
	s = 0.0;
	for (int i = 1; i < M; i++){
		for (int j = 1; j < M; j++){
			if (i >= j){
				s = 0.0;
				for (int k = 0; k < i; k++) {
					s += LL[i][k] * U[k][j];
				}
				LL[i][j] = A2[i][j] - s;
			}
			else{
			   s = 0.0;
			   for (int k = 0; k < i; k++) {
				   s += LL[i][k] * U[k][j];
			   }
				U[i][j] = (A2[i][j] - s) / LL[i][i];
			}
		}
	}
	for (int i = 0; i < M; i++) {//���������� ������ ����
		y2[i] = B2[i];
		for (int j = 0; j < i; j++) {
			y2[i] -= LL[i][j] * y2[j];
		}
		y2[i] /= LL[i][i];
	}
	
	for (int i = M-1; i >=0; i--) {//���������� ����� ����� � ����� ����� ��� � � ���� ��������
		x2[i] = y2[i];
		for (int j = i+1; j < M; j++) {
			x2[i] -= U[i][j] * x2[j];
		}
		x2[i] /= U[i][i];
	}
	cout << endl;
	cout << endl;
	cout << "reshenie LU" << endl;
	for (int i = 0; i < M; i++) {
		cout << x2[i] << " ";
	}
///////////////////////////////��������////////////////////////////////////////
	double A3[M][M];//��������
	double x3[M];//������� ��������� 
	double B3[M];//������ ����� � 
	double y3[M];//��� ��������� ���� ��� ����������
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			A3[i][j] = arr2[i][j];
		}
	}
	for (int i = 0; i < M; i++) {
		x3[i] = 0.0;
		y3[i] = 0.0;
	}
	for (int i = 0; i < M; i++) {
		B3[i] = arrright(i);
	}
	double F[M][M];//�������  -����������������
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			F[i][j] = 0.0;
		}
	}
	for (int i = 0; i < M; i++) {
		double s = A3[i][i];
		for (int k = 0; k < i; k++) {
			s -= pow(F[i][k], 2);
		}
		F[i][i] = sqrt(s);
		for (int j = i + 1; j < M; j++) {
			s = A3[j][i];
			for (int k = 0; k < i; k++) {
				s -= F[i][k] * F[j][k];
			}
			F[j][i] = s / F[i][i];
		}
	}
	for (int i = 0; i < M; i++) {//����������
		y3[i] = B3[i];
		for (int j = 0; j < i; j++) {
			y3[i] -= F[i][j] * y3[j];
		}
		y3[i] /=F[i][i];
	}

	for (int i = M - 1; i >= 0; i--) {//����������
		x3[i] = y3[i];
		for (int j = i + 1; j < M; j++) {
			x3[i] -= F[j][i] * x3[j];
		}
		x3[i] /= F[i][i];
	}
	cout << endl;
	cout << endl;
	cout << "reshenie Xolezky" << endl;
	for (int i = 0; i < M; i++) {
		cout << x3[i] << " ";
	}
	/////////////////////////////////////////////////////////����� ������� ����������///////////////////////////////////
	double A_R[M][M];
	double B_R[M];
	double x_R[M];//�� i+1
	double x_pred_R[M];//�� i

	for (int i = 0; i < M; i++) {
		x_R[i] = 0.0;
		x_pred_R[i] = 0.0;
	}

	for (int i = 0; i < M; i++) {
		B_R[i] = arrright(i);
	}

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			A_R[i][j] = arr2[i][j];
		}
	}

	double w = 1.05;
	int iter = 0;
	double delta = 1.0;
	double e = 1e-10;

	while ((delta > e) && (iter < 1000)) {
		delta = 0.0;
		for (int i = 0; i < M; i++) {
			double sum1 = 0.0;
			double sum2 = 0.0;
			for (int j = 0; j < i; j++) {
				sum1 += w * A_R[i][j] * x_R[j] / A_R[i][i];
			}
			for (int j = i + 1; j < M; j++) {
				sum2 += w * A_R[i][j] * x_pred_R[j] / A_R[i][i];
			}
			x_R[i] = (1 - w) * x_pred_R[i] - sum1 - sum2 + w * B_R[i] / A_R[i][i];
		}

		for (int i = 0; i < M; i++) {
			delta += pow(abs(x_R[i] - x_pred_R[i]), 2);

		}
		delta = sqrt(delta);

		for (int i = 0; i < M; i++) {
			x_pred_R[i] = x_R[i];//x(i)=x(i+1)
		}
		iter++;
	}

	cout << endl;
	cout << endl;
	cout << "kol-vo iter RELAX=" << iter;
	cout << endl;
	cout << "reshenie RELAX" << endl;
	for (int k = 0; k < M; k++) {
		cout << x_R[k] << ' ';
	}
//////////////////////////////////////���������� ����������///////////////////////////////////
	//�� ��� � pred ��������� � "�", ��� ����- "�+1"
	double A_G[M][M];//������� A ����� ������
	double B_G[M];//������ �����
	double x_G[M];//�� �+1
	double x_pred_G[M];//����������� �� �
	double* p_pred = new double[M];//pk ����� �����������
	double* r_pred = new double[M];// ������� �� k
	double* r = new double[M];//����� ������� �� �+1
	double alfa = 0.0;
	double betta = 0.0;
	double* A_mult_p = new double[M];//������������ ������� � �� ������ p A*pk

	for (int i = 0; i < M; i++){
		x_G[i] = 0.0;
		x_pred_G[i] = 0.0;//��������� ����������� - 0
		r[i] = 0.0;
		p_pred[i] = 0.0;
		r_pred[i] = 0.0;
	}

	for (int i = 0; i < M; i++){//��������� ������ �����
		B_G[i] = arrright(i);
	}

	for (int i = 0; i < M; i++){//��������� ������� �
		for (int j = 0; j < M; j++){
			A_G[i][j] = arr2[i][j];
		}
	}

	for (int i = 0; i < M; i++){
		r_pred[i] = B_G[i];//r0=B-Ax0=B-0
		p_pred[i] = r_pred[i];//po=r0
	}

	delta = 1.0;//����� ����� � � ����
	iter = 0; //����� ��������
	e = 1e-6;

	double norm = 0.0;
	for (int i = 0; i < M; i++) {
		if (abs(r_pred[i]) > norm) {
			norm = abs(r_pred[i]);
		}
	}
	while ((norm > e)&&(iter<1000)) {
		double mult_vect = 0.0;//����������� ��� �������� ����������� �����
	    double mult_r_r_pred = 0.0;//����������� ��� �������� ��������� �����
		double mult_r_r = 0.0;// ������������ ��� ������� ��������� �����
		delta = 0.0;

		for (int i = 0; i < M; i++) {
			A_mult_p[i] = 0;
		}

		//������� ������������ ������� �� ������ A*p� 
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				A_mult_p[i] += A_G[i][j] * p_pred[j];
			}
		}

		for (int i = 0; i < M; i++) {
			mult_vect += p_pred[i] * A_mult_p[i];//pk*(Ap�)
			mult_r_r_pred += r_pred[i] * r_pred[i];//rk*rk
		}


		alfa = mult_r_r_pred / mult_vect;

		for (int i = 0; i < M; i++) {
			x_G[i] = x_pred_G[i] + alfa * p_pred[i];//����� �����������
			r[i] = r_pred[i] - alfa * A_mult_p[i];//����� �������

		}

		for (int i = 0; i < M; i++) {
			mult_r_r += r[i] * r[i];//������������ ��� ������� ����� � ��������� r(k+1)*r(k+1)
		}

		betta = mult_r_r / mult_r_r_pred;

		for (int i = 0; i < M; i++) {//����� ������ �����������
			p_pred[i] = r[i] + betta * p_pred[i];
		}

		norm = 0.0;
		for (int i = 0; i < M; i++) {
			if (abs(r[i]) > norm) {
				norm = abs(r[i]);
			}
		}

		for (int i = 0; i < M; i++) {
			x_pred_G[i] = x_G[i];//x(k)=x(k+1)
			r_pred[i] = r[i];//r(�)=r(�+1)
		}
		iter++;
	}

	cout << endl;
	cout << endl;
	cout <<"kol-vo iter SOPR_G=" << iter;
	cout << endl;
	cout << "reshenie SOPR_G" << endl;
	for (int k = 0; k < M; k++) {
		cout << x_G[k] << ' ';
	}


	
	///////////////////////////////�����������//////////////////////////////////////
	double e1_abs_G = 0.0, e1_otn_G = 0.0, e2_abs_G = 0.0, e2_otn_G = 0.0, einf_abs_G = 0.0, einf_otn_G = 0.0;
	double e1_abs_LU = 0.0, e1_otn_LU = 0.0, e2_abs_LU = 0.0, e2_otn_LU = 0.0, einf_abs_LU = 0.0, einf_otn_LU = 0.0;
	double e1_abs_X = 0.0, e1_otn_X = 0.0, e2_abs_X = 0.0, e2_otn_X = 0.0, einf_abs_X = 0.0, einf_otn_X = 0.0;
	double e1_abs_R = 0.0, e1_otn_R = 0.0, e2_abs_R = 0.0, e2_otn_R = 0.0, einf_abs_R = 0.0, einf_otn_R = 0.0;
	double e1_abs_SG = 0.0, e1_otn_SG = 0.0, e2_abs_SG = 0.0, e2_otn_SG = 0.0, einf_abs_SG = 0.0, einf_otn_SG = 0.0;

	pogreshnost1(e1_abs_G, e1_otn_G, e2_abs_G, e2_otn_G, einf_abs_G, einf_otn_G, M, x, x1);//x-�������� � ���������,x1-������
	pogreshnost1(e1_abs_LU, e1_otn_LU, e2_abs_LU, e2_otn_LU, einf_abs_LU, einf_otn_LU, M, x, x2);
	pogreshnost1(e1_abs_X, e1_otn_X, e2_abs_X, e2_otn_X, einf_abs_X, einf_otn_X, M, x, x3);
	pogreshnost1(e1_abs_R, e1_otn_R, e2_abs_R, e2_otn_R, einf_abs_R, einf_otn_R, M, x, x_R);
	pogreshnost1(e1_abs_SG, e1_otn_SG, e2_abs_SG, e2_otn_SG, einf_abs_SG, einf_otn_SG, M, x, x_G);

	ofstream fout12;
	fout12.open("file12.txt");
	fout12.setf(ios::scientific); 
	fout12 << "POGRESHNOST" << endl;

	fout12 << endl;
	fout12 << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(23) << left << " " << setw(30) << left << "|| . ||1" << setw(30) << left << "|| . || 2" << setw(30) << left << "|| . || inf " << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "metod" << setw(15) << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs " << setw(15) << left << "otn" << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "Gauss" << setw(15) << left << e1_abs_G << setw(15) << left << e1_otn_G << setw(15) << left << e2_abs_G << setw(15) << left << e2_otn_G << setw(15) << left << einf_abs_G << setw(15) << left << einf_otn_G << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "LU" << setw(15) << left << e1_abs_LU << setw(15) << left << e1_otn_LU << setw(15) << left << e2_abs_LU << setw(15) << left << e2_otn_LU << setw(15) << left << einf_abs_LU << setw(15) << left << einf_otn_LU << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "Holezkiy" << setw(15) << left << e1_abs_X << setw(15) << left << e1_otn_X << setw(15) << left << e2_abs_X << setw(15) << left << e2_otn_X << setw(15) << left << einf_abs_X << setw(15) << left << einf_otn_X << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "Ralax" << setw(15) << left << e1_abs_R << setw(15) << left << e1_otn_R << setw(15) << left << e2_abs_R << setw(15) << left << e2_otn_R << setw(15) << left << einf_abs_R << setw(15) << left << einf_otn_R << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout12 << setw(15) << left << "Sopr Grad" << setw(15) << left << e1_abs_SG << setw(15) << left << e1_otn_SG << setw(15) << left << e2_abs_SG << setw(15) << left << e2_otn_SG << setw(15) << left << einf_abs_SG << setw(15) << left << einf_otn_SG << endl;
	fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;

	fout12.close();
	/////////////////////////////////////////////����� �������////////////////////////////////////////////////////
	///Ax-b=r, ����� ������� � ������� ������� ������������ �*x ����� ������� ����� Ax � B
	double Ax_G[M], Ax_LU[M], Ax_X[M], Ax_R[M], Ax_SG[M];
	for (int i = 0; i < M; i++) {
		Ax_G[i] = 0.0;
		Ax_LU[i] = 0.0;
		Ax_X[i] = 0.0;
		Ax_R[i] = 0.0;
		Ax_SG[i] = 0.0;
	}
	
	proizv(arr2, Ax_G, x1, M);
	proizv(arr2, Ax_LU, x2, M);
	proizv(arr2, Ax_X, x3, M); 
	proizv(arr2, Ax_R, x_R,M);
	proizv(arr2, Ax_SG, x_G,M);

	double e1_abs_nev_G = 0.0, e1_otn_nev_G = 0.0, e2_abs_nev_G = 0.0, e2_otn_nev_G = 0.0, einf_abs_nev_G = 0.0, einf_otn_nev_G = 0.0;
	double e1_abs_nev_LU = 0.0, e1_otn_nev_LU = 0.0, e2_abs_nev_LU = 0.0, e2_otn_nev_LU = 0.0, einf_abs_nev_LU = 0.0, einf_otn_nev_LU = 0.0;
	double e1_abs_nev_X = 0.0, e1_otn_nev_X = 0.0, e2_abs_nev_X = 0.0, e2_otn_nev_X = 0.0, einf_abs_nev_X = 0.0, einf_otn_nev_X = 0.0;
	double e1_abs_nev_R = 0.0, e1_otn_nev_R = 0.0, e2_abs_nev_R = 0.0, e2_otn_nev_R = 0.0, einf_abs_nev_R = 0.0, einf_otn_nev_R = 0.0;
	double e1_abs_nev_SG = 0.0, e1_otn_nev_SG = 0.0, e2_abs_nev_SG = 0.0, e2_otn_nev_SG = 0.0, einf_abs_nev_SG = 0.0, einf_otn_nev_SG = 0.0;

	pogreshnost1(e1_abs_nev_G, e1_otn_nev_G, e2_abs_nev_G, e2_otn_nev_G, einf_abs_nev_G, einf_otn_nev_G, M, B3, Ax_G);
	pogreshnost1(e1_abs_nev_LU, e1_otn_nev_LU, e2_abs_nev_LU, e2_otn_nev_LU, einf_abs_nev_LU, einf_otn_nev_LU, M, B3, Ax_LU);
	pogreshnost1(e1_abs_nev_X, e1_otn_nev_X, e2_abs_nev_X, e2_otn_nev_X, einf_abs_nev_X, einf_otn_nev_X, M, B3, Ax_X);
	pogreshnost1(e1_abs_nev_R, e1_otn_nev_R, e2_abs_nev_R, e2_otn_nev_R, einf_abs_nev_R, einf_otn_nev_R, M, B3, Ax_R);
	pogreshnost1(e1_abs_nev_SG, e1_otn_nev_SG, e2_abs_nev_SG, e2_otn_nev_SG, einf_abs_nev_SG, einf_otn_nev_SG, M, B3, Ax_SG);

	cout << endl;
	//cout << e1_abs_nev_G << " " << e1_otn_nev_G << " " << e2_abs_nev_G << " " << e2_otn_nev_G << " " << einf_abs_nev_G << " " << einf_otn_nev_G;
	ofstream fout1;
	fout1.open("file1.txt");
	fout1.setf(ios::scientific);
	fout1 << "NEVYAZKA NORMA"<<endl;
	fout1 << endl;
	fout1 << endl;
	fout1 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout1 << setw(23) << left << " " << setw(30) << left << "|| . ||1" << setw(30) << left << "|| . || 2" << setw(30) << left << "|| . || inf " << endl;
	fout1 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout1 << setw(15) << left << "metod" << setw(15) << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs " << setw(15) << left << "otn" << endl;
	fout1 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout1 << setw(15) << left << "Gauss" << setw(15) << left << e1_abs_nev_G << setw(15) << left << e1_otn_nev_G << setw(15) << left << e2_abs_nev_G << setw(15) << left << e2_otn_nev_G << setw(15) << left << einf_abs_nev_G << setw(15) << left << einf_otn_nev_G << endl;
	fout1 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout1 << setw(15) << left << "LU" << setw(15) << left << e1_abs_nev_LU << setw(15) << left << e1_otn_nev_LU << setw(15) << left << e2_abs_nev_LU << setw(15) << left << e2_otn_nev_LU << setw(15) << left << einf_abs_nev_LU << setw(15) << left << einf_otn_nev_LU << endl;
	fout1 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout1 << setw(15) << left << "Holezkiy" << setw(15) << left << e1_abs_nev_X << setw(15) << left << e1_otn_nev_X << setw(15) << left << e2_abs_nev_X << setw(15) << left << e2_otn_nev_X << setw(15) << left << einf_abs_nev_X << setw(15) << left << einf_otn_nev_X << endl;
	fout1 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout1 << setw(15) << left << "Ralax" << setw(15) << left << e1_abs_nev_R << setw(15) << left << e1_otn_nev_R << setw(15) << left << e2_abs_nev_R << setw(15) << left << e2_otn_nev_R << setw(15) << left << einf_abs_nev_R << setw(15) << left << einf_otn_nev_R << endl;
	fout1 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	fout1 << setw(15) << left << "Sopr Grad" << setw(15) << left << e1_abs_nev_SG << setw(15) << left << e1_otn_nev_SG << setw(15) << left << e2_abs_nev_SG << setw(15) << left << e2_otn_nev_SG << setw(15) << left << einf_abs_nev_SG << setw(15) << left << einf_otn_nev_SG << endl;
	fout1 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
	
	fout1.close();

///////////////////////////////////////////////////�������� ��������/////////////////////////////////////////////////////////
	delete[] arrxL;
	delete[] arrxM;
	delete[] arryL;
	
	for (int i = 0; i < M; i++) {
		delete[] arr[i];
	}
	delete[] arr;

	for (int i = 0; i < M; i++) {
		delete[] arr2[i];
	}
	delete[] arr2;
	delete[] p_pred;//pk ����� �����������
	delete[] r_pred;// ������� �� k
	delete[] r;
	return 0;
}
void pogreshnost1(double& e1_abs, double& e1_othn, double& e2_abs, double& e2_othn, double& einf_abs, double& einf_othn, int M, double* x_resh, double* x) {
	for (int i = 0; i < M; i++){
		e1_abs += fabs(x_resh[i] - x[i]);
		e1_othn += fabs(x_resh[i]);
		e2_abs += pow(fabs(x_resh[i] - x[i]), 2);
		e2_othn += pow(fabs(x_resh[i]), 2);
		einf_abs = max(einf_abs, fabs(x_resh[i] - x[i]));
		einf_othn = max(einf_othn, fabs(x_resh[i]));
	}

	e1_othn = e1_abs / e1_othn;
	e2_abs = sqrt(e2_abs);
	e2_othn = sqrt(e2_othn);
	e2_othn = e2_abs / e2_othn;
	einf_othn = einf_abs / einf_othn;
}

void proizv(double** A, double* Ax, double* x, int M) {//������������ ������� �� �����
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			Ax[i] += A[i][j] * x[j];
		}
	}
}
