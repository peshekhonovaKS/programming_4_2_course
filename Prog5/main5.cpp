#include <fstream> 
#include <iostream>
#include <time.h>
#include <iomanip>//дл€ вывода красивого в файлики
#include <cmath>
using namespace std;
double psi0(double t) { //краевое значение слева
	return t;
}
double psi1(double t) {//краевое значение справа
	return 1 - pow(t, 2);

}
double fi(double x) {//н у
	return pow(x, 2);
}
double f(double t, double x) {
	return pow(t, 2) * cos(3.14 * x);
}
void pogreshnost1(double& e1_abs, double& e1_othn, double& e2_abs, double& e2_othn, double& einf_abs, double& einf_othn, int N, int M, double** U, double** U_05);
void pogreshnost2(double& e1_abs, double& e1_othn, double& e2_abs, double& e2_othn, double& einf_abs, double& einf_othn, int N, int M, double** U, double** U_05);

void progonka(int N, int M, double h, double tau, double** U, double* x, double* t, double* Mat1, double* Mat2, double* Mat3, double* Mat_prav, double* A, double* B) {
	//double sigma = 1.0;
	double sigma = 0.5;
	for (int i = 1; i < N; i++) {
		Mat1[0] = 1, Mat1[M - 1] = 1;//на главную диагональ
		Mat2[0] = 0;//первый эл-т верхней диагонали
		Mat3[M - 2] = 0;//последний эл-т нижней диагонали
		Mat_prav[0] = psi0(t[i]);//1 урав-ие гр условие x=x0 u(t,0)=пси0(t)
		Mat_prav[M - 1] = psi1(t[i]);//последнее x=X u(t,M-1)=пси1(t)

		for (int j = 1; j < M - 1; j++) {
			Mat2[j] = -tau * sigma / pow(h, 2);//верхн€€
			Mat3[j - 1] = -tau * sigma / pow(h, 2);//нижн€€
			Mat1[j] = 2 * tau * sigma / pow(h, 2) + 1;
			Mat_prav[j] = U[i-1][j] + tau * (U[i-1][j + 1] - 2 * U[i-1][j] + U[i-1][j - 1]) * (1 - sigma) / (h * h) + tau * (1 - sigma) * f(t[i-1], x[j]) + tau * sigma * f(t[i], x[j]);
		
		}

		A[0] = -Mat2[0] / Mat1[0];
		B[0] = Mat_prav[0] / Mat1[0];

		
		for (int j = 1; j <M ; j++) {
			A[j] = -Mat2[j] / (Mat1[j] + Mat3[j - 1] * A[j - 1]);
			B[j] = (Mat_prav[j] - Mat3[j - 1] * B[j - 1]) / (Mat1[j] + Mat3[j - 1] * A[j - 1]);
		}

		U[i][M - 1] = B[M - 1];//послежнее уравнение
		for (int j = M - 2; j >= 0; j--) {
			U[i][j] = B[j] +A[j] * U[i][j + 1];
		}
	}
	
}



int main() {

	double x0 = 0.0; //xc(0,1) по усл
	double X = 1.0;
	double t0 = 0.0; //tс[0,1]
	double T = 1.0;

	//////////////////////////////////////////////////////////////////////////////////////////////////////дл€ сетки h///////////////////////////////////////////////////////
	int M = 11;
	double h = (X - x0) / (M - 1);
	int N = 1 + 2 * pow((M - 1), 2); //из услови€ что тау<=h^2/2(тогда €вна€ схема устойчива) (на листочке всЄ расписано)
	double tau = (T - t0) / (N - 1);

    ////////////////////////////////////////////////////////дл€ не€вной схемы
	int N_n = 1 + 2*(M - 1);////дл€ не€вной схемы
	double tau_n = (T - t0) / (N_n - 1);

	//зададим U дл€ €вной схемы
	double** U = new double* [N];//кол-во строк по времени t
	for (int i = 0; i < N; i++) {
		U[i] = new double[M];//кол-во столбцов по x
		for (int j = 0; j < M; j++) {
			U[i][j] = 0.0;   //U[i][j]=U[t][x]=U[N][M]
		}
	}
	//зададим U дл€ не€вной схемы с весами+прогонки
	double** U2 = new double* [N_n];//кол-во строк по времени t
	for (int i = 0; i < N_n; i++) {
		U2[i] = new double[M];//кол-во столбцов по x
		for (int j = 0; j < M; j++) {
			U2[i][j] = 0.0;   //U[i][j]=U[t][x]=U[N][M]
		}
	}
	/////////////////////////////дл€ €вной схемы
	//зададим сетку по x
	double* x = new double[M];
	for (int i = 0; i < M; i++) {//равномерна€ сетка с шагом h
		x[i] = x0 + i * h;
	}
	//зададим сетку по t
	double* t = new double[N];
	for (int i = 0; i < N; i++) {//равномерна€ сетка с шагом tau
		t[i] = t0 + i * tau;
	}
	//зададим граничные услви€ дл€ U
	for (int j = 0; j < M; j++) {
		U[0][j] = fi(x[j]); //ну t=0
		U2[0][j] = fi(x[j]);
	}
	for (int i = 0; i < N; i++) {
		U[i][0] = psi0(t[i]);//x=0
		U[i][M - 1] = psi1(t[i]);//x=X
	}
	/////////////////////////////дл€ не€вной схемы
	//зададим сетку по t
	double* t_n = new double[N_n];
	for (int i = 0; i < N_n; i++) {//равномерна€ сетка с шагом tau
		t_n[i] = t0 + i * tau_n;
	}
	for (int i = 0; i < N_n; i++) {
		U2[i][0] = psi0(t_n[i]);
		U2[i][M - 1] = psi1(t_n[i]);
	}



	
////////////////////////////////////////////////////////////////////////сетка h/2////////////////////////////////////////////////////////////////////
	int M_05 = 2 * M - 1;
	int N_05 = 1 + 2 * pow((M_05 - 1), 2);
	double h_05 = (X - x0) / (M_05 - 1);
	double tau_05 = (T - t0) / (N_05 - 1);

	////////////////////дл€ не€вной   
	int N_05_n = 1 + 2 *(M_05 - 1);
	double tau_05_n = (T - t0) / (N_05_n - 1);


	//зададим U дл€ €вной схемы
	double** U_05 = new double* [N_05];//кол-во строк по времени t
	for (int i = 0; i < N_05; i++) {
		U_05[i] = new double[M_05];//кол-во столбцов по x
		for (int j = 0; j < M_05; j++) {
			U_05[i][j] = 0.0;   //U[i][j]=U[t][x]=U[N][M]
		}
	}
	
	//зададим U дл€ не€вной схемы с весами+прогонки
	double** U2_05 = new double* [N_05_n];//кол-во строк по времени t
	for (int i = 0; i < N_05_n; i++) {
		U2_05[i] = new double[M_05];//кол-во столбцов по x
		for (int j = 0; j < M_05; j++) {
			U2_05[i][j] = 0.0;   //U[i][j]=U[t][x]=U[N][M]
		}
	}

	//зададим сетку по x
	double* x_05 = new double[M_05];
	for (int i = 0; i < M_05; i++) {//равномерна€ сетка с шагом h
		x_05[i] = x0 + i * h_05;
	}
	//зададим сетку по t
	double* t_05 = new double[N_05];
	for (int i = 0; i < N_05; i++) {//равномерна€ сетка с шагом tau
		t_05[i] = t0 + i * tau_05;
	}
	//зададим граничные услви€ дл€ U
	for (int j = 0; j < M_05; j++) {
		U_05[0][j] = fi(x_05[j]);
		U2_05[0][j] = fi(x_05[j]);
	}
	for (int i = 0; i < N_05; i++) {
		U_05[i][0] = psi0(t_05[i]);;
		U_05[i][M_05 - 1] = psi1(t_05[i]);
	}
	//////////////////////////////дл€ не€вной
	//зададим сетку по t
	double* t_05_n = new double[N_05_n];
	for (int i = 0; i < N_05_n; i++) {//равномерна€ сетка с шагом tau
		t_05_n[i] = t0 + i * tau_05_n;
	}
	//зададим граничные услви€ дл€ U
	for (int i = 0; i < N_05_n; i++) {
		U2_05[i][0] = psi0(t_05_n[i]);
		U2_05[i][M_05 - 1] = psi1(t_05_n[i]);
	}

///////////////////////////////////////////////////////////////////////////////я¬Ќјя —’≈ћј////////////////////////////	
		// дл€ сетки h
		for (int i = 0; i < N-1 ; i++){//врем€
			for (int j = 1; j < M - 1; j++){//иксы
				U[i + 1][j] = U[i][j] + (tau * (U[i][j + 1] - 2 * U[i][j] + U[i][j - 1]))/pow(h,2) + tau * f(t[i], x[j]);
			}
		}
		//дл€ сетки h/2
		for (int i = 0; i < N_05 - 1; i++){ 
			for (int j = 1; j < M_05 - 1; j++){
				U_05[i + 1][j] = U_05[i][j] + (tau_05 * (U_05[i][j + 1] - 2 * U_05[i][j] + U_05[i][j - 1])) / pow(h_05, 2) + tau_05 * f(t_05[i], x_05[j]);
			}
		}

	//////////////////////////////////////////////////////////Ќе€вна€ схема с весами////////////////////////////////////////////////////
	///////////////////////сетка h/////////////////////////////////
   double* Mat1 = new double[M]; //средн€€ диагональ
   double* Mat2 = new double[M - 1]; //верхн€€ диагональ
   double* Mat3 = new double[M - 1]; //нижн€€ диагональ
   double* Mat_prav = new double[M]; //правые части
   double* A = new double[M];
   double* B = new double[M];
   for (int i = 0; i < M; i++) {
	   A[i] = 0.0;
	   B[i] = 0.0;
	   Mat1[i] = 0.0;
	   Mat_prav[i] = 0.0;
   }
   for (int i = 0; i < M-1; i++) {
	   Mat3[i] = 0.0;
	   Mat2[i] = 0.0;
   }
  progonka(N_n, M, h, tau_n, U2, x, t_n, Mat1, Mat2, Mat3, Mat_prav, A, B);
   ////////////////////////////cетка h/2////////////////////////////
   double* Mat1_05 = new double[M_05]; //средн€€ диагональ
   double* Mat2_05 = new double[M_05 - 1]; //верхн€€ диагональ
   double* Mat3_05 = new double[M_05 - 1]; //нижн€€ диагональ
   double* Mat_prav_05 = new double[M_05]; //правые части
   double* A_05 = new double[M_05];
   double* B_05 = new double[M_05];
   for (int i = 0; i < M_05; i++) {
	   A_05[i] = 0.0;
	   B_05[i] = 0.0;
	   Mat1_05[i] = 0.0;
	   Mat_prav_05[i] = 0.0;
   }
   for (int i = 0; i < M_05 - 1; i++) {
	   Mat3_05[i] = 0.0;
	   Mat2_05[i] = 0.0;
   }
  progonka(N_05_n, M_05, h_05, tau_05_n, U2_05, x_05, t_05_n, Mat1_05, Mat2_05, Mat3_05, Mat_prav_05, A_05, B_05);


   cout << endl;
   cout << endl;
   cout << endl;

  /* for (int p = 0; p < N; p++) {//строки матрицы ј1
	   for (int i = 0; i < M; i++) {//столбец матрицы ј1
		   cout << U2[p][i] << " ";
	   }
	   cout << endl;
   }*/
 ////////////////////////////////////////////////////ѕќ√–≈ЎЌќ—“»//////////////////////////////////////////////////////////
   double e1_abs_1 = 0.0,  //  e1 -||.||1 и тд, 1-явна€ схема
	   e1_otn_1 = 0.0,
	   e2_abs_1 = 0.0,
	   e2_otn_1 = 0.0,
	   einf_abs_1 = 0.0,
	   einf_otn_1 = 0.0;
   double e1_abs_2 = 0.0,  //  e1 -||.||1 и тд, 2-Ќе€вна€ схема
	   e1_otn_2 = 0.0,
	   e2_abs_2 = 0.0,
	   e2_otn_2 = 0.0,
	   einf_abs_2 = 0.0,
	   einf_otn_2 = 0.0;

  pogreshnost1(e1_abs_1, e1_otn_1, e2_abs_1, e2_otn_1, einf_abs_1, einf_otn_1, N ,M, U, U_05);
  pogreshnost2(e1_abs_2, e1_otn_2, e2_abs_2, e2_otn_2, einf_abs_2, einf_otn_2, N_n, M, U2, U2_05);
  cout << endl;
  cout << endl;
  cout << endl;

  cout << e1_abs_1 << ' ' << e1_otn_1 << ' ' << e2_abs_1 << ' ' << e2_otn_1 << ' ' << einf_abs_1 << ' ' << einf_otn_1 << endl;
  cout << e1_abs_2 << ' ' << e1_otn_2 << ' ' << e2_abs_2 << ' ' << e2_otn_2 << ' ' << einf_abs_2 << ' ' << einf_otn_2 << endl;
  ////////////////////////////////////////////////‘ј…Ћџ//////////////////////////////////////////////////////////////////////
  ofstream fout1;
  fout1.open("file1.txt");
  for (int i = 0; i < M; i++) {
	  fout1 << fixed;
	  fout1.precision(5);
	  fout1.setf(ios::right);
	  fout1.width(10);
	  fout1 << x[i] << "\n";
  }
  fout1.close();
  ofstream fout2;
  fout2.open("file2.txt");
  for (int i = 0; i < N; i++) {
	  fout2 << fixed;
	  fout2.precision(5);
	  fout2.setf(ios::right);
	  fout2.width(10);
	  fout2 << t[i] << "\n";
  }
  fout2.close();

  ofstream fout7;
  fout7.open("file7.txt");
  for (int i = 0; i < N_n; i++) {
	  fout7 << fixed;
	  fout7.precision(5);
	  fout7.setf(ios::right);
	  fout7.width(10);
	  fout7 << t_n[i] << "\n";
  }
  fout7.close();


  ofstream fout3;
  fout3.open("file3.txt");
  for (int i = 0; i < N; i++) {
	  for (int j = 0; j < M; j++) {
		  fout3 << fixed;
		  fout3.precision(6);
		  //fout3.setf(ios::right);
		  //fout3.width(7);
		  fout3 << U[i][j] << " ";
	  }
	  fout3<< endl;
  }
  fout3.close();
  ofstream fout4;
  fout4.open("file4.txt");
  for (int i = 0; i < N_05; i++) {
	  for (int j = 0; j < M_05; j++) {
		  fout4 << fixed;
		  fout4.precision(7);
		  //fout4.setf(ios::right);
		  //fout4.width(7);
		  fout4 << U_05[i][j] << " ";
	  }
	  fout4 << endl;
  }
  fout4.close();

  ofstream fout5;
  fout5.open("file5.txt");
  for (int i = 0; i < N_n; i++) {
	  for (int j = 0; j < M; j++) {
		  fout5 << fixed;
		  fout5.precision(6);
		  //fout4.setf(ios::right);
		  //fout4.width(7);
		  fout5 << U2[i][j] << " ";
	  }
	  fout5 << endl;
  }
  fout5.close();

  ofstream fout6;
  fout6.open("file6.txt");
  for (int i = 0; i < N_05_n; i++) {
	  for (int j = 0; j < M_05; j++) {
		  fout6 << fixed;
		  fout6.precision(6);
		  //fout4.setf(ios::right);
		  //fout4.width(7);
		  fout6 << U2_05[i][j] << " ";
	  }
	  fout6 << endl;
  }
  fout6.close();

  ofstream fout12;
  fout12.open("file12.txt");
  fout12.setf(ios::scientific);
  fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  fout12 << setw(27) << left << " " << setw(30) << left << "|| . ||1" << setw(30) << left << "|| . || 2" << setw(30) << left << "|| . || inf " << endl;
  fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  fout12 << setw(20) << left << "metod" << setw(15) << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs" << setw(15) << left << "otn" << setw(15) << left << "abs " << setw(15) << left << "otn" << endl;
  fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  fout12 << setw(20) << left << "iavnaya shema" << setw(15) << left << e1_abs_1 << setw(15) << left << e1_otn_1 << setw(15) << left << e2_abs_1 << setw(15) << left << e2_otn_1 << setw(15) << left << einf_abs_1 << setw(15) << left << einf_otn_1 << endl;
  fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  fout12 << setw(20) << left << "neiavnaya shema" << setw(15) << left << e1_abs_2 << setw(15) << left << e1_otn_2 << setw(15) << left << e2_abs_2 << setw(15) << left << e2_otn_2 << setw(15) << left << einf_abs_2 << setw(15) << left << einf_otn_2 << endl;
  fout12 << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  fout12.close();















	for (int i = 0; i < N; i++) {
		delete[] U[i];
	}
	delete[] U;
	for (int i = 0; i < N_05; i++) {
		delete[] U_05[i];
	}
	delete[] U_05;
	for (int i = 0; i < N_n; i++) {
		delete[] U2[i];
	}
	delete[] U2;
	for (int i = 0; i < N_05_n; i++) {
		delete[] U2_05[i];
	}
	delete[] U2_05;
	delete[] t;
	delete[] x;
	delete[] t_05;
	delete[] x_05;
	delete[] t_n;
	delete[] t_05_n;
	delete[] Mat1;
	delete[] Mat2;
	delete[] Mat3;
	delete[] Mat_prav;
	delete[] A;
	delete[] B;
	delete[] Mat1_05;
	delete[] Mat2_05;
	delete[] Mat3_05;
	delete[] Mat_prav_05;
	delete[] A_05;
	delete[] B_05;
	


















	std::system("python PythonApplication1.py");

	return 0;
}
void pogreshnost1(double& e1_abs, double& e1_othn, double& e2_abs, double& e2_othn, double& einf_abs, double& einf_othn, int N, int M, double** U, double** U_05)
{
	for (int i = 0; i < N; i++){
		for (int j = 0; j < M; j++){
			e1_abs += fabs(U[i][j] - U_05[4 * i][2 * j]);
			//e1_othn += fabs(U[i][j]);
			e1_othn += fabs(U_05[4*i][2*j]);
			e2_abs += pow(fabs(U[i][j] - U_05[4* i][2 * j]), 2);
			//e2_othn += pow(fabs(U[i][j]), 2);
			e2_othn += pow(fabs(U_05[4*i][2*j]), 2);
			einf_abs = max(einf_abs, fabs(U[i][j] - U_05[4* i][2 * j]));
			//einf_othn = max(einf_othn, fabs(U[i][j]));
			einf_othn = max(einf_othn, fabs(U_05[4*i][2*j]));
		}
	}

	e1_othn = e1_abs / e1_othn;
	e2_abs = sqrt(e2_abs);
	e2_othn = sqrt(e2_othn);
	e2_othn = e2_abs / e2_othn;
	einf_othn = einf_abs / einf_othn;
}
void pogreshnost2(double& e1_abs, double& e1_othn, double& e2_abs, double& e2_othn, double& einf_abs, double& einf_othn, int N, int M, double** U, double** U_05)//дл€ не€вной схемы
{
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			e1_abs += fabs(U[i][j] - U_05[2 * i][2 * j]);
			//e1_othn += fabs(U[i][j]);
			e1_othn += fabs(U_05[2 * i][2 * j]);
			e2_abs += pow(fabs(U[i][j] - U_05[2 * i][2 * j]), 2);
			//e2_othn += pow(fabs(U[i][j]), 2);
			e2_othn += pow(fabs(U_05[2 * i][2 * j]), 2);
			einf_abs = max(einf_abs, fabs(U[i][j] - U_05[2 * i][2 * j]));
			//einf_othn = max(einf_othn, fabs(U[i][j]));
			einf_othn = max(einf_othn, fabs(U_05[2 * i][2 * j]));
		}
	}

	e1_othn = e1_abs / e1_othn;
	e2_abs = sqrt(e2_abs);
	e2_othn = sqrt(e2_othn);
	e2_othn = e2_abs / e2_othn;
	einf_othn = einf_abs / einf_othn;
}