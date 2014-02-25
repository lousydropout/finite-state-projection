package body Numerics is
   
   function "+" (Left  : in Natural_Vector;
		 Right : in Int_Vector) return Natural_Vector is
      Result : Natural_Vector (Left'Range);
   begin
      for I in Left'Range loop
	 Result (I) := Natural (Left (I) + Right (I + Right'First - Left'First));
      end loop;
      return Result;
   end "+";
   
   function "+" (Left, Right : in String) return String is
      Result : String (1 .. Left'Length + Right'Length);
   begin
      for I in Left'Range loop
	 Result (I - Left'First + Result'First) := Left (I);
      end loop;
      for I in Right'Range loop
	 Result (I - Right'First + Result'First + Left'Length) := Right (I);
      end loop;
      return Result;
   end "+";

end Numerics;
