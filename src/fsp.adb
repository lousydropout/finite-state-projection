with Ada.Text_IO; use Ada.Text_IO;
with Numerics; use Numerics;
with Numerics.Real_Arrays; use Numerics.Real_Arrays;
with Matrix_Exponential; use Matrix_Exponential;
with FSP_Package;
procedure Fsp is
   
   package Int_IO  is new Ada.Text_IO.Integer_IO (Integer); use Int_IO;
   package Real_IO is new Ada.Text_IO.Float_IO   (Real);    use Real_IO;
   
   N     : constant Positive := 5;
   N_Rxn : constant Positive := 10;
   
   Min : constant Int_Vector := (0, 0, 0, 0, 0);
   Max : constant Int_Vector := (1, 1, 1, 1, 15);
   
   package FSP_Pkg is new FSP_Package (N, N_Rxn, Min, Max); use FSP_Pkg;
   
   X    : Real_Vector (1 .. Array_Size) := (others => 0.0);
   Rxns : Reaction_List;
   
   Zero : constant Species := (0, 0, 0, 0, 0);
   G1   : constant Species := (1, 0, 0, 0, 0);
   G2   : constant Species := (0, 1, 0, 0, 0);
   G3   : constant Species := (0, 0, 1, 0, 0);
   G4   : constant Species := (0, 0, 0, 1, 0);
   PapI : constant Species := (0, 0, 0, 0, 1);
   

   
   U0 : constant Real := 100.0;
   
   function C1 (X : in Species) return Real is
   begin
      return 1.0 * U0;
   end C1;
   function C2 (X : in Species) return Real is
      R : Real := Real (X (5));
   begin
      return 2.5 - 2.25 * R / (1.0 + R);
   end C2;
   function C3 (X : in Species) return Real renames C1; 
   function C4 (X : in Species) return Real is
      R : Real := Real (X (5));
   begin
      return 1.2 - 0.2 * R / (1.0 + R);
   end C4;
   function C5 (X : in Species) return Real is
   begin
      return 0.01 * (U0 - 1.0);
   end C5;
   function C6 (X : in Species) return Real renames C4;
   function C7 (X : in Species) return Real renames C5;
   function C8 (X : in Species) return Real renames C2;
   function CT (X : in Species) return Real is
   begin
      if X (2) > 0 then
	 return 10.0;
      else
	 return 0.0;
      end if;
   end CT;
   function CD (X : in Species) return Real is
   begin
      return 1.0;
   end CD;
   
   
   function Forward_Euler (P  : Real_Vector;
			   Dt : Real) return Real_Vector is
      Mat : Real_Matrix := Propensity_Matrix (Rxns);
      Result : Real_Vector := P + Dt * (Mat * P);
   begin
      return Result;
   end Forward_Euler;
   
   K    : Real := 0.0;
   Time : Real := 0.0;
   Dt   : constant Real := 0.001;
   F    : File_Type;
begin
   -- Create file
   Create (File => F, Name => "data.csv");
   
   Rxns (1) := Create_Rxn (G1, G2, C1'Access);
   Rxns (2) := Create_Rxn (G2, G1, C2'Access);
   Rxns (3) := Create_Rxn (G1, G3, C3'Access);
   Rxns (4) := Create_Rxn (G3, G1, C4'Access);
   Rxns (5) := Create_Rxn (G2, G4, C5'Access);
   Rxns (6) := Create_Rxn (G4, G2, C6'Access);
   Rxns (7) := Create_Rxn (G3, G4, C7'Access);
   Rxns (8) := Create_Rxn (G4, G3, C8'Access);
   Rxns (9) := Create_Rxn (G2, G2 + PapI, CT'Access);
   Rxns (10) := Create_Rxn (PapI, Zero, CD'Access);
   
   
   X (Number (G1 + 5 * PapI)) := 1.0;
   
   while Time < 10.0 loop
      X := Forward_Euler (X, Dt);
      Time := Time + Dt;
   end loop;
   
   Put_Line (F, "PapI, g1, g2, g3, g4");
   for I in Min (5) .. Max (5) loop
      Put (F, I, 0); Put (F, ", ");
      Put (F, X (Number (G1 + I * PapI)), Exp => 3, Aft => 4); Put (F, ", ");
      Put (F, X (Number (G2 + I * PapI)), Exp => 3, Aft => 4); Put (F, ", ");
      Put (F, X (Number (G3 + I * PapI)), Exp => 3, Aft => 4); Put (F, ", ");
      Put (F, X (Number (G4 + I * PapI)), Exp => 3, Aft => 4); New_Line (F);
   end loop;
   
   for I of X loop
      K := K + I;
   end loop;
   Put (K); New_Line;
end Fsp;
