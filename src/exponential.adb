with Ada.Text_IO; use Ada.Text_IO;
with Numerics; use Numerics;
with Numerics.Real_Arrays; use Numerics.Real_Arrays;
with Matrix_Exponential; use Matrix_Exponential;
with Ada.Containers; use Ada.Containers;
with Ada.Containers.Vectors;
with Ada.Containers.Indefinite_Hashed_Maps; -- , Ada.Strings.Hash;
procedure Exponential is
   package Int_IO is new Ada.Text_IO.Integer_IO (Integer); use Int_IO;
   package Real_IO is new Ada.Text_IO.Float_IO (Real); use Real_IO;
   
   --  package Int_Vector_Package is 
   --     new Ada.Containers.Vectors (Positive, Integer, "=");
   --  package Real_Vector_Package is 
   --     new Ada.Containers.Vectors (Positive, Real   , "=");
   --  use Int_Vector_Package, Real_Vector_Package;
   
   
   --  subtype Int_Vector  is Int_Vector_Package.Vector;
   --  subtype Vector is Real_Vector_Package.Vector;

   
   function Hash (N : Natural) return Hash_Type is
   begin
      return Hash_Type (N);
   end Hash;
   
   package Int_Mat is new Ada.Containers.Indefinite_Hashed_Maps
     (Key_Type        => Natural,
      Element_Type    => Natural,
      Hash            => Hash,
      Equivalent_Keys => "=");
   
   
   
   N : constant Natural := 6;
   type Int_Array is array (1 .. N) of Integer;
   N_Array : constant Int_Array := (99, 104, 203, 7, 87, 103);
   
   function Number (Item : in Int_Array) return Integer is
      Result : Integer := 0;
   begin
      for I in reverse Item'Range loop
	 Result := Result * N_Array (I) + Item (I);
      end loop;
      return Result;
   end Number;
   
   function Inv_Number (I : in Integer) return Int_Array is
      Result : Int_Array;
      Tmp    : Integer := I;
      Tmp2   : Integer;
   begin
      for J in Result'Range loop
	 Tmp2 := Tmp / N_Array (J);
	 Result (J) := Tmp - Tmp2 * N_Array (J);
	 Tmp := Tmp2;
      end loop;
      return Result;
   end Inv_Number;
   
   function "+" (Left, Right : in Int_Array) return Int_Array is
      Result : Int_Array;
   begin
      for I in Result'Range loop
	 Result (I) := Left (I) + Right (I);
      end loop;
      return Result;
   end "+";
   
   function "=" (Left, Right : in Int_Array) return Boolean is
   begin
      for I in Left'Range loop
	 if Left (I) /= Right (I) then return False; end if;
      end loop;
      return True;
   end "=";
   
   Diff : Int_Array := (1, 0, -1, 0, 0, 0);
   X    : Int_Array := (3, 0,  0, 1, 2, 0);
   Z    : Int_Array := X + Diff;
   
   Num : Integer := Number (X) + Number (Diff);
   Y : Int_Array := Inv_Number (Num);
   
   A : Real_Matrix (1 .. 2, 1 .. 2) 
     := ((1.0, 2.0),
	 (-3.0, 4.0));
   
   
   package Vector_Package is 
      new Ada.Containers.Vectors (Positive, Int_Array, "=");
   
   --  List : Vector := To_Vector (100);
begin
   Put (Num); New_Line;
   for I in 1 .. N loop
      Put (Z (I)); Put (Y (I)); New_Line;
   end loop;
   New_Line;
   
   A := Exp (-A);
   for I in 1 .. 2 loop
      for J in 1 .. 2 loop
	 Put (A (I, J), Exp => 0, Aft => 5); Put ("   ");
      end loop;
      New_Line;
   end loop;
   New_Line;
   null;
end Exponential;
