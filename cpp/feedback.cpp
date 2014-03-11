#include <fstream>
#include <Eigen/Dense>
#include <iostream>
#include <math.h>

using namespace Eigen;
using namespace std;


const double epsilon = 1.0e-3;
/** Not Problem Specific.  For Problem Statement Scroll further down.
                  Various Definitions, Functions, and Procedures  **/
typedef struct
{
  VectorXd Reactant;
  VectorXd Product;
  VectorXd Change;
  double (* Coeff) (VectorXd);
} Reaction;

int Array_Size (const int Min[], const int Max[], const int N)
{
  int res = 1;
  for (int i = 0; i < N; i++)
    res = res * (Max [i] - Min [i] + 1);
  return res;
}

void Set_Size (Reaction &X, const int N)
{
  X.Reactant.resize (N);
  X.Product.resize (N);
  X.Change.resize (N);
}

void Create_Rxn (Reaction &X,
		 VectorXd Reactant,
		 VectorXd Product,
		 double func (VectorXd))
{
  X.Reactant = Reactant;
  X.Product  = Product;
  X.Change   = Product - Reactant;
  X.Coeff    = func;
}

int Number (VectorXd X, const int Min[], const int Max[], const int N)
{
  int res = 0;
  for (int i = N - 1; i >= 0; i--)
    {
      res = res * (Max [i] - Min [i] + 1) + X (i) - Min [i];
    }
  return res;
}

VectorXd Inv_Number (int I, const int Min[], const int Max[], const int N)
{
  VectorXd res (N);
  int tmp = I;
  int tmp2, length;

  for (int j = 0; j < N; j++)
    {
      length  = Max [j] - Min [j] + 1;
      tmp2    = tmp / length;
      res (j) = tmp - tmp2 * length + Min [j];
      tmp     = tmp2;
    }
  return res;
}

double Propensity (Reaction Rxn, VectorXd X, const int N)
{
  double res = Rxn.Coeff (X);
  double tmp = 0.0;

  for (int j = 0; j < N; j++)
    {
      if (Rxn.Reactant (j) > 0)
	{
	  tmp = (double) pow (X (j), Rxn.Reactant (j));
	  res = res * tmp;
	}
    }
  return res;
}

bool Valid (VectorXd X, const int N)
{
  for (int i = 0; i < N; i++)
    {
      if (X (i) < 0)
	  return false;
    }
  return true;
}

MatrixXd Propensity_Matrix (Reaction Rxn[],
			    const int N_Rxn,
			    const int Min[],
			    const int Max[],
			    const int N)
{
  VectorXd X (N), Y (N), Z (N);
  double A;
  int K, N_max = Array_Size (Min, Max, N);
  MatrixXd Mat = MatrixXd::Zero (N_max, N_max);

  for (int i = 0; i < N_max; i++)
    {
      X = Inv_Number (i, Min, Max, N);
      for (int r = 0; r < N_Rxn; r++)
	{
	  Z = X - Rxn [r].Reactant;
	  Y = X + Rxn [r].Change;
	  if (Valid (Z, N) and Valid (Y, N))
	    {
	      A          = Propensity (Rxn [r], X, N);
	      K          = Number (Y, Min, Max, N);
	      Mat (i, i) = Mat (i, i) - A;
	      if (K < N_max) Mat (K, i) = Mat (K, i) + A;
	    }
	}
    }
  return Mat;
}


VectorXd Forward_Euler (VectorXd P,
			double Dt,
			MatrixXd Mat)
{
  return P + Dt * (Mat * P);
}


VectorXd Reverse_Euler (VectorXd P,
			MatrixXd Mat)
{
  // MatrixXd Mat2 = Mat;
  // VectorXd res = Mat2.fullPivLu ().solve (P);
  return Mat * P;
}



VectorXd Adaptive_FE (VectorXd X,
		      double time,
		      double t_final,
		      Reaction Rxn[],
		      const int N_Rxn,
		      const int Min[],
		      const int Max[],
		      const int N)
{
  MatrixXd Mat = Propensity_Matrix (Rxn, N_Rxn, Min, Max, N);
  VectorXd X1, X2;
  double Dt1 = 0.001, Dt2 = 0.5 * Dt1;
  double error;
  while (time < t_final)
    {
      if (Dt1 > t_final - time)
  	{
  	  Dt1 = t_final - time;
	  Dt2 = 0.5 * Dt1;
  	  Mat = Propensity_Matrix (Rxn, N_Rxn, Min, Max, N);
  	}
      X1 = Forward_Euler (X, Dt1, Mat);
      X2 = Forward_Euler (X, Dt2, Mat);
      X2 = Forward_Euler (X2, Dt2, Mat);

      error = (X2 - X1).norm ();
      if (error < epsilon)
	{
	  X = X2;
	  time = time + Dt1;
	  // Update time step based
	  Dt1 = Dt1 * sqrt (epsilon / error);
	  Dt2 = 0.5 * Dt1;
	  cout << "time = " << time << " new dt = " << Dt1 << endl;
	}
      else
	{
	  Dt1 = 0.5 * Dt1;
	  Dt2 = 0.5 * Dt2;
	}
    }
  return X;
}


VectorXd Adaptive_BE (VectorXd X,
		      double time,
		      double t_final,
		      Reaction Rxn[],
		      const int N_Rxn,
		      const int Min[],
		      const int Max[],
		      const int N)
{
  const int N_max = Array_Size (Min, Max, N);
  const MatrixXd Id  = MatrixXd::Identity (N_max, N_max);
  const MatrixXd Mat = Propensity_Matrix (Rxn, N_Rxn, Min, Max, N);
  
  MatrixXd Mat1, Mat2;
  VectorXd X1, X2;
  double Dt1 = 1.0e-5;
  double Dt2 = 0.5 * Dt1;
  double error;

  Mat1 = Id - Dt1 * Mat;
  Mat2 = Id - Dt2 * Mat;
  
  while (time < t_final)
    {
      if (Dt1 > t_final - time)
  	{
  	  Dt1 = t_final - time;
	  Dt2 = 0.5 * Dt1;
	  cout << "     Reducing time step to " << Dt1 << endl;
	  Mat1 = Id - Dt1 * Mat;
	  Mat2 = Id - Dt2 * Mat;
  	}
      X1 = Mat1.partialPivLu ().solve (X);
      X2 = Mat2.partialPivLu ().solve (X);
      X2 = Mat2.partialPivLu ().solve (X2);
      error = (X2 - X1).norm ();
      
      if (error < epsilon)
	{ // then accept new values
	  X = X2; time = time + Dt1;
	  // Update time step estimate
	  Dt1 = Dt1 * sqrt (epsilon / error);
	  Dt2 = 0.5 * Dt1;
	  cout << "time = " << time << "  dt = " << Dt1 << endl;
	}
      else
	{ // reduce time step by 1/5
	  Dt1 = 0.2 * Dt1; Dt2 = 0.2 * Dt2;
	  Mat1 = Id - Dt1 * Mat;
	  Mat2 = Id - Dt2 * Mat;

	  cout << "     reducing time step to " << Dt1 << endl;
	}
    }
  return X;
}
/**************************************************************/












// Problem Definition

const int N = 6;       // Number of Species
const int N_Rxns = 10; // Number of Reactions


// Define Reaction Coefficient Functions
const double alpha = 1.0e10; // Just some constant
const double P = 1.0e-7; // Just some constant
double C1 (VectorXd X) { return 1.0 / 6.0; }
double C2 (VectorXd X) { return 1.0e8; }
double C3 (VectorXd X) { return C2 (X) / alpha; }
double C4 (VectorXd X) { return 1.0e-4; }
double C5 (VectorXd X) { return 1.8e10 * P; }
double C6 (VectorXd X) { return 1.0; }
double C7 (VectorXd X) { return 1.0e-2; }
double C8 (VectorXd X) { return 1.0 / 6.0; }
double C9 (VectorXd X) { return 1.0e9; }
double C10 (VectorXd X) { return 1.0; }




int main ()
{

  // Set Boundaries
  int Min[N] = {0, 0, 0, 0, 0, 0};
  int Max[N] = {1, 1, 1, 4, 4, 4};
  const int N_max = Array_Size (Min, Max, N);
  

  // Create data.csv
  // ofstream file;
  // file.open("data.csv"); 

  
  // Define Species
  VectorXd Zero (N), DNA (N), DNA_P (N), DNA_Tf2 (N), mRNA (N), Tf (N), Tf2 (N);
  Zero    << 0, 0, 0, 0, 0, 0;
  DNA     << 1, 0, 0, 0, 0, 0;
  DNA_P   << 0, 1, 0, 0, 0, 0;
  DNA_Tf2 << 0, 0, 1, 0, 0, 0;
  mRNA    << 0, 0, 0, 1, 0, 0;
  Tf      << 0, 0, 0, 0, 1, 0;
  Tf2     << 0, 0, 0, 0, 0, 1;


  cout << "Defining Rxns . . .";
  // Define Reactions
  Reaction Rxns [N_Rxns];
  for (int i = 0; i < N_Rxns; i++)  Set_Size (Rxns [i], N); // Set VectorXd size
  Create_Rxn (Rxns [0], DNA_P, DNA + mRNA, C1);
  Create_Rxn (Rxns [1], DNA + Tf2, DNA_Tf2, C2);
  Create_Rxn (Rxns [2], DNA_Tf2, DNA + Tf2, C3);
  Create_Rxn (Rxns [3], Tf, Zero, C4);
  Create_Rxn (Rxns [4], DNA, DNA_P, C5);
  Create_Rxn (Rxns [5], DNA_P, DNA, C6);
  Create_Rxn (Rxns [6], mRNA, Zero, C7);
  Create_Rxn (Rxns [7], mRNA, mRNA + Tf, C8);
  Create_Rxn (Rxns [8], 2 * Tf, Tf2, C9);
  Create_Rxn (Rxns [9], Tf2, 2 * Tf, C10);
  cout << "Finished" << endl;

  cout << "Setting Initial Values . . .";
  // Initial Value
  VectorXd X = VectorXd::Zero (N_max);
  int num = Number (DNA, Min, Max, N);
  X (num) = 1.0;
  cout << "Finished" << endl;


  
  cout << "Start Integration" << endl;
  // Integrate
  double time = 0.0;
  double t_final = 10.0;

  X = Adaptive_FE (X, time, t_final, Rxns, N_Rxns, Min, Max, N);
  // X = Adaptive_BE (X, time, t_final, Rxns, N_Rxns, Min, Max, N);

  
  // Output Results to file
  // file << "PapI, g1, g2, g3, g4" << endl;
  // for (int I = Min [N - 1]; I < Max [N - 1]; I++)
  //   {
  //     file << I << ", ";
  //     file << X (Number (G1 + I * PapI, Min, Max, N)) << ", ";
  //     file << X (Number (G2 + I * PapI, Min, Max, N)) << ", ";
  //     file << X (Number (G3 + I * PapI, Min, Max, N)) << ", ";
  //     file << X (Number (G4 + I * PapI, Min, Max, N)) << endl;
  //   }



  // Output Probability Sum to Monitor
  double sum = 0.0;
  for (int i = 0; i < N_max; i++) sum = sum + X (i);
  cout << "Probability Sum = " << sum * 100.0d << "%" << endl;



  return 0;
}
