// The MIT License
// Copyright (c) 2011 Pedro http://lamehacks.net
// Copyright (c) 2012 Ophir LOJKINE
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

"use strict";

var linalg = {

	Matrix : function (elems) {
		if (elems === undefined) {return linalg.Identity(3);}
		else if (! elems instanceof Array) {throw "InvalidType";}
		else if (typeof(elems[0]) !== "object") {return new linalg.Vector(elems);}

		this.rows = elems.length;
		this.columns = elems[0].length;
		for (var i=0; i<this.rows; i++) {
			if (elems[i].length != this.columns) throw "ImpossibleSize";
		}
		this.elems = elems;
		
		this.toVector = function() {
			if (this.rows = 1) return linalg.Vector(this.elems[0]);
			else if (this.columns = 1) return linalg.Vector(this.transpose()[0]);
			else throw "NotAVector";
		}
		
		/* /!\ Add another matrix in place. It alters the matrix */
		this.addInPlace = function (m2) {
			if (m2.rows!==this.rows || m2.columns!==this.columns) {
				throw "SizesDoNotMatch";
			}
			for (var i = 0; i < this.rows; i++) {
				for (var j = 0; j < this.columns; j++) {
				    this.elems[i][j] += m2.elems[i][j]
				}
			}
			return this;
		};

		this.add = function (m2) {
			var m1 = this.clone();
			return m1.addInPlace(m2);
		};

		this.multiply = function (m2) {
			if ( ! m2 instanceof linalg.Matrix) throw "IncompatibleTypes";
			if (m2.rows != this.columns) throw "SizesDoNotMatch";
			var result = [];
			for (var i=0; i<this.rows; i++) {
				var row = new Array(m2.columns);
				for (var j=0; j<m2.columns; j++) {
					row[j] = 0;
					for (var k=0; k<this.columns; k++) {
						row[j] += this.elems[i][k]*m2.elems[k][j];
					}
				}
				result.push(row);
			}
			return new linalg.Matrix(result);
		};

		/*Warning : In-place operations*/
		this.scalarMultiplyInPlace = function (s) {
			for (var k = this.rows*this.columns-1; k>=0; k--) {
				this.elems[parseInt(k/this.columns)][k%this.columns] *= s;
			}
			return this;
		}

		this.scalarMultiply = function (m2) {
			var m1 = this.clone();
			return m1.scalarMultiplyInPlace(m2);
		};

		/* /!\ Transpose a square matrix in place. It alters the matrix */
		this.transposeFast = function () {
			if ( this.rows !== this.columns) throw "NonSquareMatrix";
			for (var i=0; i<this.rows; i++) {
				for (var j=i+1; j<this.columns; j++) {
					var t = this.elems[i][j];
					this.elems[i][j] = this.elems[j][i];
					this.elems[j][i] = t;
				}
			}
			return this;
		};

		this.transpose = function() {
			var t = [];
			for (var i = 0; i < this.columns; i++) {
				var row = [];
				for (var j = 0; j < this.rows; j++) {
				    row.push(this.elems[j][i]);
				}
				t.push(row);
			}
			return (new linalg.Matrix(t));
		};

		
		this.toString = function () {
			var str="";
			for (var i=0;i<this.rows;i++){
				str+="(";
				for (var j=0;j<this.columns;j++){
					str += this.elems[i][j];
					if(j!==this.columns-1) str+="\t";
				}
				str+= ")\n";
			}
			return str;
		};

		this.clone = function () {
			var elems = this.elems;
			var func = function(i, j) {return elems[i][j];};
			return linalg.generateMatrix(this.rows, this.columns, func);
		};

	},
	
	
	generateMatrix : function (nlines, ncols, elements) {
		var func;
		if (typeof(elements) === "number"){
		    func = function(){return elements;};
		}else{
		    func = elements;
		}
		var elems = [];
		for (var i = 0; i < nlines; i++) {
		    var row = [];
		    for (var j = 0; j < ncols; j++) {
		        row.push(func(i, j));
		    }
		    elems.push(row);
		}
		return (new linalg.Matrix(elems));
	},

	Vector : function (elems, rowVector) {
		if (elems === undefined) return new linalg.Vector([0]);
		if (typeof(elems) !== "object") throw "InvalidType";
		if ( elems instanceof linalg.Matrix) return elems.toVector();
		var m = new linalg.Matrix(new Array(elems));
		if (rowVector) return m;
		else return m.transpose();
	},

};

function matrixMultiply(a, b) {
    var btrans = transposeMatrix(b);
    var result = [];
    for (var i = 0; i < a.length; i++) {
        var row = [];
        for (var j = 0; j < btrans.length; j++) {
            var value = internalProduct(a[i], btrans[j]);
            row.push(value);
        }
        result.push(row);
    }
    return result;
}

function matrixScalarMultiply(m, s) {
    var result = [];
    for (var i = 0; i < m.length; i++) {
        var row = [];
        for (var j = 0; j < m[0].length; j++) {
            row.push(s * m[i][j]);
        }
        result.push(row);
    }
    return result;
}


//unnecessary - use dotOp instead


function matrixAdd(a, b) {
    var result = [];
    for (var i = 0; i < a.length; i++) {
        var row = [];
        for (var j = 0; j < a[0].length; j++) {
            row.push(a[i][j] + b[i][j]);
        }
        result.push(row);
    }
    return result;
}

//for internal usage only


function internalProduct(u, v) {
    if (u.length != v.length) throw "SizesDoNotMatch";
    var sum = 0;
    for (var i = 0; i < u.length; i++) {
        sum += u[i] * v[i];
    }
    return sum;
}

function transposeMatrix(m) {
    var t = [];
    for (var i = 0; i < m[0].length; i++) {
        var row = [];
        for (var j = 0; j < m.length; j++) {
            row.push(m[j][i]);
        }
        t.push(row);
    }
    return t;
}

function minorMatrix(m, k, l) {
    var reduced = [];
    for (var i = 0; i < m.length; i++) {
        if (i == k) continue;
        var row = [];
        for (var j = 0; j < m.length; j++) {
            if (j == l) continue;
            row.push(m[i][j]);
        }
        reduced.push(row);
    }
    return reduced;
}


/*
Recursive implementation using laplace expansion
http://www.webcitation.org/61AGedZlm
*/

function determinant(m) {
    var size = m.length;
    if (size == 1) return m[0][0];
    if (size == 2) return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    var det = 0;
    for (var i = 0; i < size; i++) {
        var minor = minorMatrix(m, 0, i);
        var signal = (i % 2 > 0) ? -1 : 1;
        det += signal * m[0][i] * determinant(minor);
    }
    return det;
}


/*
http://en.wikipedia.org/wiki/Cofactor_(linear_algebra)
*/

function cofactor(m, k, l) {
    var minor = minorMatrix(m, k, l);
    return determinant(minor);
}

/*
http://en.wikipedia.org/wiki/Cofactor_(linear_algebra)#Matrix_of_cofactors
*/

function cofactorMatrix(m) {
    var cofactors = [];
    for (var i = 0; i < m.length; i++) {
        var row = [];
        for (var j = 0; j < m.length; j++) {
            var cofactorval = cofactor(m, i, j);
            row.push(cofactorval);
        }
        cofactors.push(row);
    }
    return cofactors;
}


/*
Used the adjoint method
http://www.webcitation.org/61BTRqAoZ
*/

function inverseMatrix(m) {
    var det = determinant(m);
    if (det === 0) throw "SingularMatrix";
    var deti = 1 / det;
    var cof = cofactorMatrix(m);
    var adj = transposeMatrix(cof);
    var result = matrixScalarMultiply(adj, deti);
    return result;
}

/*
performs operation element by element between to matrices
*/

function dotOp(func, m, n) {
    var result = [];
    for (var i = 0; i < m.length; i++) {
        var row = [];
        for (var j = 0; j < m[0].length; j++) {
            row.push(func(m[i][j], n[i][j]));
        }
        result.push(row);
    }
    return result;
}

/*Returns the matrix obtained by applying the function func
on each couple of factors m(i,k) and m(j,k).
For instance, with i=1, j=2, and func=f
(1 2)    ( f(1,3)[0] f(3,4)[0] )
(3 4) -> ( f(1,3)[1] f(3,4)[1] )
(5 6)    (    5         6      )
*/
function rowOperation (m, i, j, func) {
    var h = m.length, w = m[0].length;
    if ( ! ( 0<=i<h || 0<=j<h) ) {
        throw "IncompatibleIndexes";
    }
    for (var k=0; k<w; k++){
        var r = func(m[i][k], m[j][k]);
        m[i][k] = r[0];
        m[j][k] = r[1];
    }
    return m;
}

/*Swap rows i and j in the matrix m*/
function swapRows(m, i, j) {
    return rowOperation (m, i, j, function(factor1, factor2) {
        return [factor2, factor1];
    });
}

/*Substract lambda * the jth row to the ith row*/
function substractAndMultiply(m, i, j, lambda) {
    return rowOperation(m, i, j, function (factor1, factor2) {
        return [factor1-lambda*factor2, factor2];
    });
}

/*Multiply the ith row by a factor of lambda*/
function rowMultiply (m, i, lambda) {
    return rowOperation(m, i, 0, function (factor1, factor2) {
        return [lambda*factor1, factor2];
    });
}

/*Use the Gauss-Jordan elimination on m.
If the matrix n is provided, the same row operations will be applied to n.
So if n is the identity matrix, it will become the inverse of m*/
function gauss_jordan_elimination (m, n) {
    var curRow, curColumn, h=m.length, w=m[0].length;
    for (curColumn=0; curColumn<Math.min(w,h); curColumn++) {
        //Find the row on which we will start the elimination
        for (var i=curColumn; i<h; i++) {
            if (m[i][curColumn]!==0) break;
        }
        if (i<h) curRow=i;
        else continue;
        swapRows(m, curRow, curColumn);
        if (n) swapRows(n, curRow, curColumn);
        curRow = curColumn;

    	var lambda =  1/m[curRow][curColumn];
		rowMultiply(m, curRow, lambda);
        if (n) rowMultiply(n, curRow, lambda);

        for (var i=0; i<h; i++) {
        	if (i==curRow) continue;
            var lambda = m[i][curColumn]/m[curRow][curColumn];
            substractAndMultiply(m, i, curRow, lambda); 
            if (n) substractAndMultiply(n, i, curRow, lambda);
        }
    }
    return m;
}

function fasterInverse (m) {
    var size=0;
    if (m.length != m[0].length) throw "SingularMatrix";
    else size = m.length;
    var mcopy = clone(m);
    var inverse = identity(size);
    gauss_jordan_elimination(mcopy, inverse);
    return inverse;    
}

/*Solve A.x = B*/
function linearSolve(A, B) {
    var Acopy = clone(A);

    var Bcopy;
	if (typeof(B[0])==="number") Bcopy = vector(B);
    else if (B) Bcopy = clone(B);
    else Bcopy = zeros(A.length, 1);

    gauss_jordan_elimination(Acopy, Bcopy);
    return Bcopy;
}

function generateMatrix(nlines, ncols, elements) {
    var func;
    if (typeof(elements) === "number"){
        func = function(){return elements;};
    }else{
        func = elements;
    }
    var m = [];
    for (var i = 0; i < nlines; i++) {
        var row = [];
        for (var j = 0; j < ncols; j++) {
            row.push(func(i, j));
        }
        m.push(row);
    }
    return m;
}

function clone(m) {
    var h=m.length, w=m[0].length;
    return generateMatrix(h, w, function(i, j) {
        return m[i][j];
    });
}

function vector(arr) {
    var h=arr.length;
    return generateMatrix(h, 1, function(i) {
        return arr[i];
    });
}

function zeros(nlines, ncols) {
    return generateMatrix(nlines, ncols, 0);
}

function ones(nlines, ncols) {
    return generateMatrix(nlines, ncols, 1);
}

function identity(size) {
    return generateMatrix(size, size, function(i, j) {
        if (i === j) return 1;
        else return 0;
    });
}

function randomNumberInRange (min, max) {
    return min + (max-min)*Math.random();
}

function randomMatrix (nlines, ncols, min, max) {
    if(!ncols) ncols = nlines;
    if (!min) min = 0;
    if(!max) max=1;

    return generateMatrix(nlines, ncols, function(){
        return randomNumberInRange(min, max);
    });
}

function matrix2str (m) {
	var w=m[0].length, h=m.length;
	var str="";
	for (var i=0;i<h;i++){
		str+="(";
		for (var j=0;j<w;j++){
			str += m[i][j];
			if(j!==w-1) str+="\t";
		}
		str+= ")\n";
	}
	return str;
}
