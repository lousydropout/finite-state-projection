with Ada.Text_IO, Numerics, Numerics.Real_Arrays;
use  Ada.Text_IO, Numerics, Numerics.Real_Arrays;
with FSP_Package;

procedure Fsp is
   
   package Int_IO  is new Ada.Text_IO.Integer_IO (Integer); use Int_IO;
   package Real_IO is new Ada.Text_IO.Float_IO   (Real);    use Real_IO;
   
   
   -- Define Boundaries and FSP_Package
   N     : constant Positive := 5;
   N_Rxn : constant Positive := 10;
   
   Min : constant Int_Vector := (0, 0, 0, 0, 0);
   Max : constant Int_Vector := (1, 1, 1, 1, 15);
   package FSP_Pkg is new FSP_Package (N, N_Rxn, Min, Max); use FSP_Pkg;
   -----------------------------
   
   
   -- Define Species
   Zero : constant Species := (0, 0, 0, 0, 0);
   G1   : constant Species := (1, 0, 0, 0, 0);
   G2   : constant Species := (0, 1, 0, 0, 0);
   G3   : constant Species := (0, 0, 1, 0, 0);
   G4   : constant Species := (0, 0, 0, 1, 0);
   PapI : constant Species := (0, 0, 0, 0, 1);
   
   
   -- Define Reaction Coefficient Functions
   U0 : constant Real := 100.0;
   
   function C1 (X : in Species; T : in Real) return Real is
   begin
      return 1.0 * U0;
   end C1;
   
   function C2 (X : in Species; T : in Real) return Real is
      R : Real := Real (X (5));
   begin
      return 2.5 - 2.25 * R / (1.0 + R);
   end C2;
   
   function C3 (X : in Species; T : in Real) return Real renames C1; 
   
   function C4 (X : in Species; T : in Real) return Real is
      R : Real := Real (X (5));
   begin
      return 1.2 - 0.2 * R / (1.0 + R);
   end C4;
   
   function C5 (X : in Species; T : in Real) return Real is
   begin
      return 0.01 * (U0 - 1.0);
   end C5;
   
   function C6 (X : in Species; T : in Real) return Real renames C4;
   
   function C7 (X : in Species; T : in Real) return Real renames C5;
   
   function C8 (X : in Species; T : in Real) return Real renames C2;
   
   function CT (X : in Species; T : in Real) return Real is
   begin
      return 10.0;
   end CT;
   
   function CD (X : in Species; T : in Real) return Real is
   begin
      return 1.0;
   end CD;
   
   
   
   
   -- Various Variables
   Prob : Real_Vector (1 .. Array_Size) := (others => 0.0); -- Initial Value
   Rxns : Reaction_List; -- Saves information on the reactions
   Sum  : Real := 0.0;   -- For summing up probabilities at the end
   F    : File_Type;     -- Output file
   
begin
   -- Create File
   Create (File => F, Name => "data.csv");
   
   -- Define Reactions
   Rxns (1)  := Create_Rxn (G1,   G2,        C1'Access);
   Rxns (2)  := Create_Rxn (G2,   G1,        C2'Access);
   Rxns (3)  := Create_Rxn (G1,   G3,        C3'Access);
   Rxns (4)  := Create_Rxn (G3,   G1,        C4'Access);
   Rxns (5)  := Create_Rxn (G2,   G4,        C5'Access);
   Rxns (6)  := Create_Rxn (G4,   G2,        C6'Access);
   Rxns (7)  := Create_Rxn (G3,   G4,        C7'Access);
   Rxns (8)  := Create_Rxn (G4,   G3,        C8'Access);
   Rxns (9)  := Create_Rxn (G2,   G2 + PapI, CT'Access);
   Rxns (10) := Create_Rxn (PapI, Zero,      CD'Access);
   
   -- Initial Value
   Prob (Number (G1 + 5 * PapI)) := 1.0;
   
   -- Integrate
   Prob := Integrate (Prob, 0.0, 10.0, Rxns, Reverse_Euler);
   
   -- Output Results To File
   Set_Output (F);
   Put_Line ("PapI, g1, g2, g3, g4");
   for I in Min (5) .. Max (5) loop
      Put (I, 0); Put (", ");
      Put (Prob (Number (G1 + I * PapI)), Exp => 3, Aft => 4); Put (", ");
      Put (Prob (Number (G2 + I * PapI)), Exp => 3, Aft => 4); Put (", ");
      Put (Prob (Number (G3 + I * PapI)), Exp => 3, Aft => 4); Put (", ");
      Put (Prob (Number (G4 + I * PapI)), Exp => 3, Aft => 4); New_Line;
   end loop;
   
   
   -- Output Probability Sum to Monitor
   Set_Output (Standard_Output);
   for P of Prob loop
      Sum := Sum + P;
   end loop;
   Put ("Probability Sum = ");
   Put (Sum * 100.0, Exp => 0, Aft => 3); 
   Put_Line (" %");
   
end Fsp;
