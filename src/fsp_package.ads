with Ada.Text_IO; use Ada.Text_IO;
with Numerics; use Numerics;
with Numerics.Real_Arrays; use Numerics.Real_Arrays;
with Matrix_Exponential; use Matrix_Exponential;

generic
   
   N     : Positive;
   N_Rxn : Positive;
   Min   : Int_Vector;
   Max   : Int_Vector;
   
package FSP_Package is
   
   subtype Species is Int_Vector (1 .. N);
   type Reaction is private;
   type Reaction_List is array (1 .. N_Rxn) of Reaction;
   
   type Func is access function (X : in Species) return Real;
   
   function Array_Size return Positive;
   function "+" (Left, Right : in Int_Vector) return Int_Vector;
   function "-" (Left, Right : in Int_Vector) return Int_Vector;
   function "*" (Left  : in Integer;
		 Right : in Int_Vector) return Int_Vector;
   
   function Create_Rxn (Start, Finish : in Species;
			C	      : in Func) return Reaction;
   
   function Number (X : in Species) return Positive;
   function Inv_Number (I : in Positive) return Species;
   
   function Propensity (Rxn : in Reaction;
   			I   : in Positive) return Real;
   function Propensity (Rxn : in Reaction;
   			X   : in Species) return Real;
   
   
   function Valid (X : in Species) return Boolean;
   function Propensity_Matrix (Rxn_List : in Reaction_List) return Real_Matrix;


private
   function Calc_Lengths return Species;
			  
   type Reaction is
      record
	 Reactant    : Species;
	 Product     : Species;
	 Change      : Species;
	 Coefficient : Func;
      end record;
   Length : constant Species := Calc_Lengths;
   N_Max : constant Positive := Array_Size;
   
end FSP_Package;
