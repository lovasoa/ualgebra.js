ualgebra.js is a tiny javascript linear algebra library.
It contains numerical implementations of basic matric operations.

Keep in mind that this this contains only numerical implementations. There's no type checking of any kind nor suport for rational numbers.
This is _not_ a symbolic calculus tool either.

Usage
=====

The functions accept matrices as two dimentional arrays. A matrix is represented as an array of lines, a line is simply an array.

examples

```javascript
M = new linalg.Matrix([[1,2],[3,4]]);

// transpose
T = M.transpose(); //returns the matrix [[1,3],[2,4]]

// inverse
M.inverse(); // returns [[-2,-1],[-1.5,-0.5]]

//Not implemented yet:
// determinant
determinant(M); //returns -2

//Not implemented yet:
// cofactors
cofactorMatrix(M); //returns [[4,3],[2,1]]


// element-by-element operations
//first parameter must be a function that accepts two paremeters
//this can be used for example to sum matrices:
N = linalg.Matrix([[1,1],[1,1]]);
dotOp(function(a,b){return a+b;},M,N); //returns [[2,3],[4,5]]


// basic aritmetic matrix operations

// matrix product
M.multiply(N); //returns [[3,3],[7,7]]

// scalar-matrix product
M.scalarMultiply(2) //returns [[2,4],[6,8]]


// Matrix generation utilities

linalg.NullMatrix(3,5); //returns [[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]

//Not implemented yet:
ones(3,3); //returns [[1,1,1],[1,1,1],[1,1,1]]

linalg.Identity(4); //returns [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]

```
