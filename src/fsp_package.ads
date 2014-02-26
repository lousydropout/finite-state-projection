with Numerics, Numerics.Real_Arrays;
use  Numerics, Numerics.Real_Arrays;

generic
   
   N     : Positive;
   N_Rxn : Positive;
   Min   : Int_Vector;
   Max   : Int_Vector;
   
package FSP_Package is
   
   type Integration_Method is 
     (Forward_Euler, Reverse_Euler, FE_Adaptive);
   subtype Species is Int_Vector (1 .. N);
   type Reaction is private;
   type Reaction_List is array (1 .. N_Rxn) of Reaction;
   type Func is access function (X : in Species;
				 T : in Real) return Real;
   
   
   function Array_Size return Positive;
   function "+" (Left, Right : in Int_Vector) return Int_Vector;
   function "-" (Left, Right : in Int_Vector) return Int_Vector;
   function "*" (Left  : in Integer;
		 Right : in Int_Vector) return Int_Vector;
   
   function Create_Rxn (Start, Finish : in Species;
			C	      : in Func) return Reaction;
   
   function Number (X : in Species) return Positive;
   -- Define Integrators
   function Integrate (X       : in Real_Vector;
		       T       : in Real;
		       T_Final : in Real;
		       Rxns    : in Reaction_List;
		       By      : in Integration_Method := Forward_Euler;
		       Dt      : in Real	       := 0.1)
     return Real_Vector;
   
private
   function Norm (X : in Real_Vector) return Real;
   function Calc_Lengths return Species;
   function Inv_Number (I : in Positive) return Species;
   function Propensity (Rxn : in Reaction;
			X   : in Species;
			T   : in Real) return Real;
   function Valid (X : in Species) return Boolean;
   function Propensity_Matrix (Rxn_List	: in Reaction_List;
			       Time	: in Real) return Real_Matrix;
   function Forward_Euler (P	: in Real_Vector;
			   T	: in Real;
			   Dt	: in Real;
			   Rxns	: in Reaction_List) return Real_Vector;
   
   function Reverse_Euler (P	: in Real_Vector;
			   T	: in Real;
			   Dt	: in Real;
			   Rxns	: in Reaction_List) return Real_Vector;


   type Reaction is
      record
	 Reactant    : Species;
	 Change      : Species;
	 Coefficient : Func;
      end record;
   
   Length : constant Species  := Calc_Lengths;
   N_Max  : constant Positive := Array_Size;
   
end FSP_Package;
