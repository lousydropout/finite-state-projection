
package body FSP_Package is
   
   function "+" (Left, Right : in Int_Vector) return Int_Vector is
      Result : Int_Vector (Left'Range);
   begin
      for I in Result'Range loop
	 Result (I) := Left (I) + Right (I);
      end loop;
      return Result;
   end "+";
   
   function "-" (Left, Right : in Int_Vector) return Int_Vector is
      Result : Int_Vector (Left'Range);
   begin
      for I in Result'Range loop
	 Result (I) := Left (I) - Right (I);
      end loop;
      return Result;
   end "-";
   
   function "*" (Left  : in Integer;
		 Right : in Int_Vector) return Int_Vector is
      Result : Int_Vector (Right'Range);
   begin
      for I in Result'Range loop
	 Result (I) := Left * Right (I);
      end loop;
      return Result;
   end "*";

   function Array_Size return Positive is
      Result : Positive := 1;
   begin
      for X of Length loop
	 Result := Result * X;
      end loop;
      return Result;
   end Array_Size;
   
   
   
   function Create_Rxn (Start, Finish : in Species;
			C : in Func) return Reaction is
      Result : Reaction;
   begin
      Result.Reactant    := Start;
      Result.Product     := Finish;
      Result.Coefficient := C;
      Result.Change      := Finish - Start;
      return Result;
   end Create_Rxn;
   
   
   
   function Number (X : in Species) return Positive is
      Result : Natural := 0;
   begin
      for I in reverse X'Range loop
	 Result := Result * Length (I) + (X (I) - Min (I));
      end loop;
      return Positive (Result + 1);
   end Number;
   
   function Inv_Number (I : in Positive) return Species is
      Result : Species;
      Tmp : Integer := I - 1;
      Tmp2 : Integer;
   begin
      for J in Result'Range loop
	 Tmp2 := Tmp / Length (J);
	 Result (J) := Tmp - Tmp2 * Length (J) + Min (J);
	 Tmp := Tmp2;
      end loop;
      return Result;
   end Inv_Number;
   
   function Calc_Lengths return Species is
      Result : Species;
   begin
      for I in Result'Range loop
	 Result (I) := Max (I) - Min (I) + 1;
      end loop;
      return Result;
   end Calc_Lengths;
   
   
   function Propensity (Rxn : in Reaction;
   			I   : in Positive) return Real is
      X      : Species := Inv_Number (I);
      Result : Real    := Rxn.Coefficient (X);
   begin
      for J in Rxn.Reactant'Range loop
      	 if Rxn.Reactant (J) > 0 then
      	    Result := Result * (Real (X (J)) ** Rxn.Reactant (J));
      	 end if;
      end loop;
      return Result;
   end Propensity;
   
   
   function Propensity (Rxn : in Reaction;
   			X   : in Species) return Real is
      Result : Real := Rxn.Coefficient (X);
   begin
      for J in Rxn.Reactant'Range loop
      	 if Rxn.Reactant (J) > 0 then
      	    Result := Result * Real ((X (J)) ** Rxn.Reactant (J));
      	 end if;
      end loop;
      return Result;
   end Propensity;
   
   
   function Valid (X : in Species) return Boolean is
   begin
      for I in X'Range loop
	 if X (I) < 0 then
	    return False; 
	 end if;
      end loop;
      return True;
   end Valid;
   
   
   function Propensity_Matrix (Rxn_List : in Reaction_List) return Real_Matrix is
      Result : Real_Matrix (1 .. N_Max, 1 .. N_Max) := (others => (others => 0.0));
      X, Y : Species;
      A : Real;
      K : Positive;
   begin
      for I in 1 .. N_Max loop
	 X := Inv_Number (I);
	 for R of Rxn_List loop
	    if Valid (X - R.Reactant) and then Valid (X + R.Change) then
	       A := Propensity (R, X);
	       Result (I, I) := Result (I, I) - A;
	       
	       Y := (X + R.Change);
	       K := Number (Y);
	       if (K < N_Max) and then Valid (Y) then 
		  Result (K, I) := Result (K, I) + A;
	       end if;
	    end if;
	 end loop;
      end loop;
      return Result;
   end Propensity_Matrix;
   
   
end FSP_Package;
   
