with Numerics, Numerics.Real_Arrays;
use Numerics, Numerics.Real_Arrays;
package Matrix_Exponential is
   
   function Eye (N : in Positive;
		 X : in Real	 := 1.0) return Real_Matrix;
   function Exp (A : in Real_Matrix) return Real_Matrix
     with Pre => A'Length (1) = A'Length (2);
   
   
private
   
   function Apqk (K : in Natural) return Real
     with Inline => True;
   
end Matrix_Exponential;
