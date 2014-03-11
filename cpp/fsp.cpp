#include <fstream>
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;


const double epsilon = 1.0e-6;
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
			Reaction Rxn[],
			const int N_Rxn,
			const int Min[],
			const int Max[],
			const int N)
{
  MatrixXd Mat = Propensity_Matrix (Rxn, N_Rxn, Min, Max, N);
  VectorXd res (Array_Size (Min, Max, N));
  res = P + Dt * (Mat * P);
  return res;
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
  	}
      X1 = X  + Dt1 * (Mat * X);
      X2 = X  + Dt2 * (Mat * X);
      X2 = X2 + Dt2 * (Mat * X2);

      error = (X2 - X1).norm ();
      if (error < epsilon)
	{
	  X = X2;
	  time = time + Dt1;
	  // Update time step based
	  Dt1 = Dt1 * sqrt (epsilon / error);
	  Dt2 = 0.5 * Dt1;
	  cout << "time = " << time << "  dt = " << Dt1 << endl;
	}
      else
	{
	  Dt1 = 0.2 * Dt1;
	  Dt2 = 0.2 * Dt2;
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

const int N = 5;       // Number of Species
const int N_Rxns = 10; // Number of Reactions


// Define Reaction Coefficient Functions
const double U0 = 100.0d; // Just some constant
double C1 (VectorXd X) { return 1.0 * U0; }
double C2 (VectorXd X) {
  double r = X (N - 1);
  return 2.5 - 2.25 * r / (1 + r);
}
double C3 (VectorXd X) { return 1.0 * U0; }
double C4 (VectorXd X) {
  double r = X (N - 1);
  return 1.2 - 0.20 * r / (1 + r);
}
double C5 (VectorXd X) { return 0.01 * (U0 - 1.0); }
double C6 (VectorXd X) {
  double r = X (N - 1);
  return 1.2 - 0.20 * r / (1 + r);
}
double C7 (VectorXd X) { return 0.01 * (U0 - 1.0); }
double C8 (VectorXd X) {
  double r = X (N - 1);
  return 2.5 - 2.25 * r / (1 + r);
}
double CT (VectorXd X) { return 10.0; }
double CD (VectorXd X) { return 1.0; }




int main ()
{

  // Set Boundaries
  const int Min [N] = {0, 0, 0, 0, 0};
  const int Max [N] = {1, 1, 1, 1, 20};
  const int N_max = Array_Size (Min, Max, N);
  

  // Create data.csv
  ofstream file;
  file.open("data.csv"); 

  
  // Define Species
  VectorXd Zero (N), G1 (N), G2 (N), G3 (N), G4 (N), PapI (N);
  Zero << 0, 0, 0, 0, 0;
  G1   << 1, 0, 0, 0, 0;
  G2   << 0, 1, 0, 0, 0;
  G3   << 0, 0, 1, 0, 0;
  G4   << 0, 0, 0, 1, 0;
  PapI << 0, 0, 0, 0, 1;

  
  // Define Reactions
  Reaction Rxns [N_Rxns];
  for (int i = 0; i < N_Rxns; i++)  Set_Size (Rxns [0], N); // Set VectorXd size
  Create_Rxn (Rxns [0], G1, G2, C1);
  Create_Rxn (Rxns [1], G2, G1, C2);
  Create_Rxn (Rxns [2], G1, G3, C3);
  Create_Rxn (Rxns [3], G3, G1, C4);
  Create_Rxn (Rxns [4], G2, G4, C5);
  Create_Rxn (Rxns [5], G4, G2, C6);
  Create_Rxn (Rxns [6], G3, G4, C7);
  Create_Rxn (Rxns [7], G4, G3, C8);
  Create_Rxn (Rxns [8], G2, G2 + PapI, CT);
  Create_Rxn (Rxns [9], PapI, Zero, CD);


  // Initial Value
  VectorXd X = VectorXd::Zero (N_max);
  X (Number (G1 + 5 * PapI, Min, Max, N)) = 1.0;





  // Integrate
  double time = 0.0;
  double t_final = 10.0;
  X = Adaptive_FE (X, time, t_final, Rxns, N_Rxns, Min, Max, N);
  // X = Adaptive_BE (X, time, t_final, Rxns, N_Rxns, Min, Max, N);

  
  // Output Results to file
  file << "PapI, g1, g2, g3, g4" << endl;
  for (int i = Min [N - 1]; i < Max [N - 1]; i++)
    {
      file << i << ", ";
      file << X (Number (G1 + i * PapI, Min, Max, N)) << ", ";
      file << X (Number (G2 + i * PapI, Min, Max, N)) << ", ";
      file << X (Number (G3 + i * PapI, Min, Max, N)) << ", ";
      file << X (Number (G4 + i * PapI, Min, Max, N)) << endl;
    }

  
  // Output Probability Sum to Monitor
  double sum = 0.0;
  for (int i = 0; i < N_max; i++) sum = sum + X (i);
  cout << "Probability Sum = " << sum * 100.0d << "%" << endl;

  
  return 0;
}
