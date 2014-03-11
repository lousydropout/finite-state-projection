with Ada.Text_IO, Numerics, Numerics.Real_Arrays;
use  Ada.Text_IO, Numerics, Numerics.Real_Arrays;
with FSP_Package;

procedure Feedback is
   
   package Int_IO  is new Ada.Text_IO.Integer_IO (Integer); use Int_IO;
   package Real_IO is new Ada.Text_IO.Float_IO   (Real);    use Real_IO;
   
   Test : Real_Matrix (1 .. 10_000, 1 .. 10_000) := (others => (others => 0.0));
   
   -- Define Boundaries and FSP_Package
   N     : constant Positive := 6;
   N_Rxn : constant Positive := 10;
   
   Min : constant Int_Vector := (0, 0, 0,  0,  0,  0);
   Max : constant Int_Vector := (1, 1, 1, 3, 3, 3);
   package FSP_Pkg is new FSP_Package (N, N_Rxn, Min, Max); use FSP_Pkg;
   -----------------------------
   
   
   -- Define Species
   Zero   : constant Species := (others => 0);
   DNA    : constant Species := (1 => 1, others => 0);
   DNAP   : constant Species := (2 => 1, others => 0);
   DNATF2 : constant Species := (3 => 1, others => 0);
   MRNA   : constant Species := (4 => 1, others => 0);
   TF     : constant Species := (5 => 1, others => 0);
   TF2    : constant Species := (6 => 1, others => 0);
   
   
   -- Define Reaction Coefficient Functions
   Alpha : constant Real := 1.0e12;
   P     : constant Real := 1.0e-7; -- only affects K5
   
   function K1 (X : in Species; T : in Real) return Real is
   begin
      return 1.0 / 6.0;
   end K1;
   
   function K2 (X : in Species; T : in Real) return Real is
   begin
      return 1.0e8;
   end K2;
   
   function K3 (X : in Species; T : in Real) return Real is
   begin
      return K2 (X, T) / Alpha;
   end K3;
   
   function K4 (X : in Species; T : in Real) return Real is
   begin
      return 1.0e-4;
   end K4;
   
   function K5 (X : in Species; T : in Real) return Real is
   begin
      return 1.8e10 * P;
   end K5;
   
   function K6 (X : in Species; T : in Real) return Real is
   begin
      return 1.0;
   end K6;
   
   function K7 (X : in Species; T : in Real) return Real is
   begin
      return 1.0e-2;
   end K7;
   
   function K8 (X : in Species; T : in Real) return Real is
   begin
      return 1.0 / 6.0;
   end K8;
   
   function K9 (X : in Species; T : in Real) return Real is
   begin
      return 1.0e9;
   end K9;
   
   function K10 (X : in Species; T : in Real) return Real is
   begin
      return 1.0;
   end K10;
   
   
   procedure Print_Results (F	 : in File_Type;
			    Prob : in Real_Vector) is
      Previous : File_Type := Ada.Text_IO.Current_Output;
      Values   : Species;
   begin
      Set_Output (F);

      for I in Prob'Range loop
	 Values := Inv_Number (I);
	 for J in Values'Range loop
	    Put (Values (J)); Put (", ");
	 end loop;
	 Put (Prob (I)); New_Line;
      end loop;
      
      Set_Output (Previous);
   end Print_Results;
   
   
   -- Various Variables
   Prob : Real_Vector (1 .. Array_Size) := (others => 0.0); -- Initial Value
   Rxns : Reaction_List; -- Saves information on the reactions
   Sum  : Real := 0.0;   -- For summing up probabilities at the end
   F    : File_Type;     -- Output file
   
   Mat : Real_Matrix := ((0.0, 0.0), (1.0, 0.0));
   Vec : Real_Vector := (1.1, 2.3);
begin
   Vec := Solve (Mat, Vec);
   Vec := Mat * Vec;
   for V of Vec loop
      Put_Line (V'Img);
   end loop;
   -- Create File
   Create (File => F, Name => "data.csv");
   
   -- Define Reactions 
   Rxns (1)  := Create_Rxn (DNAP,        DNA + MRNA,   K1'Access);
   Rxns (2)  := Create_Rxn (DNA + TF2,   DNATF2,       K2'Access);
   Rxns (3)  := Create_Rxn (DNATF2,      DNA + TF2,    K3'Access);
   Rxns (4)  := Create_Rxn (TF,          Zero,         K4'Access);
   Rxns (5)  := Create_Rxn (DNA,         DNAP,         K5'Access);
   Rxns (6)  := Create_Rxn (DNAP,        DNA,          K6'Access);
   Rxns (7)  := Create_Rxn (MRNA,        Zero,         K7'Access);
   Rxns (8)  := Create_Rxn (MRNA,        MRNA + TF,    K8'Access);
   Rxns (9)  := Create_Rxn (2 * TF,      TF2,          K9'Access);
   Rxns (10) := Create_Rxn (TF2,         2 * TF,      K10'Access);
   
   -- Initial Value
   -- Starts out with only 1 DNA molecule and background RNA Polymerase P
   Prob (Number (DNA)) := 1.0; 
   
   -- Print Header
   Put_Line (F, "Time, 0mRNA, 1mRNA, 2mRNA");
   Put_Line ("Begin Integrating");
   -- Integrate
   --  --  Prob := Integrate (Prob, 0.0, 1.0, Rxns, FE_Adaptive);
   --  --  Prob := Integrate (Prob, 0.0, 1.0, Rxns, FE_Adaptive);
   Prob := Integrate (Prob, 0.0, 1.0, Rxns, Reverse_Euler, 0.001);
   --  --  Prob := Integrate (Prob, 0.0, 10.0, Rxns, Forward_Euler, 0.001);
   --  Put_Line ("Finished Integrating");
   --  Print_Results (F, Prob);
   
   -- Output Probability Sum to Monitor
   Set_Output (Standard_Output);
   for P of Prob loop
      Sum := Sum + P;
   end loop;
   Put ("Probability Sum = ");
   Put (Sum * 100.0, Exp => 0, Aft => 3); 
   Put_Line (" %");
   
exception
   
   when Storage_Error =>
      Put_Line ("From Feedback");
      raise Storage_Error;   
   
end Feedback;
