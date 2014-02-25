with Ada.Text_IO;
package body Matrix_Exponential is
   
   function Eye (N : in Positive;
		 X : in Real	 := 1.0) return Real_Matrix is
      Result : Real_Matrix (1 .. N, 1 .. N) 
	:= (others => (others => 0.0));
   begin
      for I in 1 .. N loop
	 Result (I, I) := X;
      end loop;
      return Result;
   end Eye;
   
   
   -- P = Q = 15
   function Exp (A : in Real_Matrix) return Real_Matrix is
      N, D : Real_Matrix (A'Range (1), A'Range (2))
      	:= (others => (others => 0.0));
      Z : Real_Matrix (A'Range (1), A'Range (2))
	:= Eye (A'Length (1));
      Tmp : Real := 1.0;
      M : Integer := 1;
      P : constant Natural := 50;
   begin
      N := Z;
      for K in 1 .. P loop
      	 Tmp := Tmp / Real (K);
      	 Z := Z * A;
      	 N := N + Tmp * Z;
      end loop;
      return N;
      
      --  for K in 0 .. P loop
      --  	 --  Put ("iter: "); Put (K'img); New_Line;
      --  	 Tmp := Apqk (K);
      --  	 N := N + (Tmp * Real (M)) * Z;
      --  	 D := D + Tmp * Z;
      --  	 M := -M;
      --  	 Z := Z * A;
      --  end loop;
      --  return Solve (N, D);
   end Exp;
   
   
   
   
   function Apqk (K : in Natural) return Real is
      Tmp : Real := 1.0;
      P, Q : constant Natural := 50;
      Result : Real := 1.0;
   begin
      for J in 0 .. K - 1 loop
	 Tmp    := Real (P - J) / Real (P + Q - J);
	 Result := Result * Tmp;
      end loop;
      for I in 1 .. K loop
	 Result := Result / Real (I);
      end loop;
      
      return Result;
   end Apqk;
   

end Matrix_Exponential;
