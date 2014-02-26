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
   			X   : in Species;
			T   : in Real) return Real is
      Result : Real := Rxn.Coefficient (X, T);
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
   
   
   function Propensity_Matrix (Rxn_List	: in Reaction_List;
			       Time	: in Real) return Real_Matrix is
      Result : Real_Matrix (1 .. N_Max, 1 .. N_Max) 
	:= (others => (others => 0.0));
      X, Y : Species;
      A : Real;
      I : Positive;
   begin
      for J in 1 .. N_Max loop
	 X := Inv_Number (J);
	 
	 for R of Rxn_List loop
	    Y := (X + R.Change);
	    
	    if Valid (X - R.Reactant) and then Valid (Y) then
	       I := Number (Y);
	       A := Propensity (R, X, Time);
	       Result (J, J) := Result (J, J) - A;
	       if (I <= N_Max) then Result (I, J) := Result (I, J) + A; end if;
	    end if;
	 end loop;
      end loop;
      
      return Result;
   end Propensity_Matrix;
   
   
   
   -- Define Integrators
   function Forward_Euler (P	: in Real_Vector;
			   T	: in Real;
			   Dt	: in Real;
			   Rxns	: in Reaction_List) return Real_Vector is
      Mat    : Real_Matrix := Propensity_Matrix (Rxns, T);
      Result : Real_Vector := P + Dt * (Mat * P);
   begin
      return Result;
   end Forward_Euler;
   
   
   function Reverse_Euler (P	: in Real_Vector;
			   T	: in Real;
			   Dt	: in Real;
			   Rxns	: in Reaction_List) return Real_Vector is
      Mat : Real_Matrix 
	:= Unit_Matrix (P'Length) - Dt * Propensity_Matrix (Rxns, T);
   begin
      return Solve (Mat, P);
   end Reverse_Euler;
   
   
   
   function Integrate (X       : in Real_Vector;
		       T       : in Real;
		       T_Final : in Real;
		       Rxns    : in Reaction_List;
		       By      : in Integration_Method := Forward_Euler;
		       Dt      : in Real	       := 0.1)
		      return Real_Vector is
      Time     : Real        := T;
      Δt       : Real        := Dt;
      Dt1, Dt2 : Real        := Dt;
      E1, E2   : Real        := 1.0;
      C1, C2   : Real        := 1.0;
      C, Error : Real        := 100.0;
      Result   : Real_Vector := X;
      Y1, Y2   : Real_Vector := X;
      Iter     : Positive    := 1;
      
      Epsilon  : Real  := 1.0e-6;
   begin
      
      case By is
	 when Forward_Euler =>
	    while Time < T_Final loop
	       Δt     := Real'Min (Δt, T_Final - Time);
	       Result := Forward_Euler (Result, Time, Δt, Rxns);
	       Time   := Time + Δt;
	    end loop;
	    
	 when Reverse_Euler =>
	    while Time < T_Final loop
	       Δt     := Real'Min (Δt, T_Final - Time);
	       Result := Reverse_Euler (Result, Time, Δt, Rxns);
	       Time   := Time + Δt;
	    end loop;
	    
	 when FE_Adaptive =>
	    while Time < T_Final loop
	       Dt1 := Real'Min (Dt1, T_Final - Time);
	       Dt2 := Dt1 / 2.0;
	       
	       Y1 := Forward_Euler (Y1, Time, Dt1, Rxns);
	       Y2 := Forward_Euler (Y2, Time, Dt2, Rxns);
	       Y2 := Forward_Euler (Y2, Time, Dt2, Rxns);
	       
	       -- Incomplete
	       Y1   := Y2;
	       Time := Time + Dt1;
	    end loop;
	    
      end case;
      return Result;
   end Integrate;
   
   
   function Norm (X : in Real_Vector) return Real is
      Result : Real := 0.0;
   begin
      for Y of X loop
	 Result := Real'Max (Result, abs (Y));
      end loop;
      return Result;
   end Norm;
   
   function C_New (C_Current, C_Previous : in Real) return Real is
   begin
      return C_Current ** 2 / C_Previous;
   end C_New;
   
end FSP_Package;
