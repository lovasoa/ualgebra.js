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

	Matrix : function (elems, nocheck) {
		if (elems === undefined) {return linalg.Identity(3);}
		else if (! elems instanceof Array) {throw "InvalidType";}
		else if (typeof(elems[0]) !== "object") {return new linalg.Vector(elems);}

		this.rows = elems.length;
		this.columns = elems[0].length;
		
		if (!nocheck) {
			for (var i=0; i<this.rows; i++) {
				if (elems[i].length != this.columns) throw "ImpossibleSize";
			}
		}

		//slice(0) to copy the array
		this.elems = elems.slice(0);
		
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
		this.plus = this.add;

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
		this.dot = this.multiply;

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
		this.transposeFastSquare = function () {
			if ( this.rows !== this.columns) 	throw "NonSquareMatrix";
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

		/*Returns the matrix obtained by applying the function func
		on each couple of elements (i,k) and (j,k).
		For instance, with i=1, j=2, and func=f
		(1 2)    ( f(1,3)[0] f(3,4)[0] )
		(3 4) -> ( f(1,3)[1] f(3,4)[1] )
		(5 6)    (    5         6      )
		*/
		this.rowOperationInPlace = function (i, j, func) {
				if ( ! ( 0<=i<this.rows || 0<=j<this.rows) ) {
				    throw "IncompatibleIndexes";
				}
				for (var k=0; k<this.columns; k++){
				    var r = func(this.elems[i][k], this.elems[j][k]);
				    if (r[0]!==undefined) this.elems[i][k] = r[0];
				    if (r[1]!==undefined) this.elems[j][k] = r[1];
				}
		}
				//All the sub-functions here operate on the matrix directly
						/*Swap rows i and j in the matrix m*/
				this.swapRows =  function(i, j) {
						this.rowOperationInPlace (i, j, function(factor1, factor2) {
								return [factor2, factor1];
						});
				}

				/*Substract lambda * the jth row to the ith row*/
				this.substractAndMultiply = function(i, j, lambda) {
						this.rowOperationInPlace (i, j, function (factor1, factor2) {
								return [factor1-lambda*factor2, undefined];
						});
				}

				/*Multiply the ith row by a factor of lambda*/
				this.rowMultiply = function (i, lambda) {
						this.rowOperationInPlace (i, 0, function (factor1, factor2) {
								return [lambda*factor1, undefined];
						});
				}


		//Apply the Gauss-Jordan elimination on m.
		this.GaussJordanEliminationInPlace = function (nothrow) {

			var curRow, curColumn, h=this.rows, w=this.columns;
			for (curColumn=0; curColumn<Math.min(w,h); curColumn++) {
			    //Find the row on which we will start the elimination
			    for (var i=curColumn; i<h; i++) {
			        if (this.elems[i][curColumn]!==0) break;
			    }
			    if (i<h) curRow=i;
			    else if (nothrow) continue;
			    else throw "SingularMatrix";

			    this.swapRows(curRow, curColumn);
			    curRow = curColumn;
			    
					this.rowMultiply(curRow, 1/this.elems[curRow][curColumn]);

			    for (var i=0; i<h; i++) {
			    	if (i==curRow) continue;
			        var lambda = this.elems[i][curColumn]/this.elems[curRow][curColumn];
			        this.substractAndMultiply(i, curRow, lambda);
			    }
			}
		};
		
		this.GaussJordanElimination = function () {
			var m = this.clone();
			m.GaussJordanEliminationInPlace();
			return m;
		}

		this.inverse = function () {
				var size = this.rows;
				if (this.rows != this.columns) throw "SingularMatrix";
				var largeMat = this.concat(linalg.Identity(size));
				largeMat.GaussJordanEliminationInPlace();
				return largeMat.submatrix(0,size, size-1, 2*size-1);
		};

		//Submatrix composed of all coeffs between (l1,c1) and (l2,c2)
		this.submatrix = function (l1, c1, l2, c2) {
			if ( c2==undefined || c2 >= this.columns ) c2 = this.columns-1;
			if ( l2==undefined || l2 >= this.rows ) l2 = this.rows-1;

			var newColumns = c2-c1+1;
			var newRows = l2-l1+1;
			var newElems = elems.slice(l1,l2+1);
			if (c1!=0 || c2!=this.columns-1) {
				for (var l=0; l < newRows; l++){
					newElems[l] = newElems[l].slice(c1, c2+1);
				}
			}
			return (new linalg.Matrix(newElems, true));
		};
		
		//Concatenate the matrix Mat to the left of current matrix (or to the bottom if bottom is set)
		this.concat = function (Mat, bottom) {
			if ( (!bottom && Mat.rows!=this.rows) || (bottom && Mat.columns!=this.columns) ){
				throw "SizesDoNotMatch"
			}
			var newElems;
			if (bottom) {
				newElems = this.elems.concat(Mat.elems);
			} else {
				newElems = new Array(this.rows);
				for (var i=0; i<this.rows; i++) {
					newElems[i] = this.elems[i].concat(Mat.elems[i]);
				}
			}
			return (new linalg.Matrix(newElems, true));
		}
		
		this.clone = function () {
			var newElems = new Array(this.rows);
			for (var i=0; i<this.rows; i++) {
				newElems[i] = this.elems[i].slice(0);
			}
			return linalg.Matrix(newElems, true);
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
		if (elems === undefined) return new linalg.Vector([0,0,0]);
		if (typeof(elems) !== "object") throw "InvalidType";
		if ( elems instanceof linalg.Matrix) return elems.toVector();
		var m = new linalg.Matrix(new Array(elems));
		if (rowVector) return m;
		else return m.transpose();
	},
	
	
	/*Solve A.x = B*/
	linearSolve : function (A, B) {
		if (B === undefined) B = Vector(A.rows);
    var M = A.concat(B);
		M.GaussJordanEliminationInPlace();
		return M.submatrix(0,A.columns, A.rows, A.columns);
    return Bcopy;
	},

	Identity : function (size) {
		if ( ! size ) size = 3;
		return linalg.generateMatrix (size, size, function (i, j) {
			return ( (i==j) ? 1 : 0 )
		});
	},
	
	NullMatrix : function (nlines, ncols) {
		if ( ! nlines || ! ncols ) nlines = ncols = 3;
		return linalg.generateMatrix (nlines, ncols, 0);
	},
	
	NullVector : function (size) {
		return linalg.NullMatrix(size, 1);
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
