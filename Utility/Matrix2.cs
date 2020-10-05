using System;
//using System.Dynamic;
//using System.Linq.Expressions;
//using System.Linq;
using PUtil.Operators;
//using PUtil.Interfaces;
using PUtil.GenericMath;

[Serializable]
public class Matrix<T> : Tensor<T>, IEquatable<Matrix<T>>, IFormattable, IDisposable where T:new() //matrix is a rank 2 tensor
{
    public Matrix() //create scalar
    {
        dim = new int[2];
        dim[0] = 1;
        dim[1] = 1;
        values = new T[1];
        values[0] = new T();
    }

    public Matrix(int row, int column, bool identity = false)
    {
        dim = new int[2];
        dim[0] = row;
        dim[1] = column;

        values = new T[row * column];
        if(identity)
        {
            int n = 0;
            int m = 0;
            if (rows >= columns) { n = rows; m = columns; }
            else { n = columns; m = rows; }

            for (int i = 0; i < n; i++)
            {
                values[i * m + i] = Operator<int, T>.convert(1); //we try to fill it up with ones
            }
        }
    }
    public Matrix(int row, int col, T[] vals)
    {
        if (row <= 0 || col <= 0) throw new System.Exception("Can't have negative rows or columns");
        if (vals.Length != row * col) throw new System.Exception("Length of array doesn't match dimensions");
        dim = new int[2];
        rows = row;
        columns = col;
        values = new T[row * columns];
        for (int i = 0; i < vals.Length; i++)
        {
            values[i] = vals[i];
        }
    }
    public Matrix(Matrix<T> m)
    {
        dim = new int[2];
        rows = m.rows;
        columns = m.columns;
        values = new T[m.values.Length];
        for(int i=0; i<m.values.Length; i++)
        {
            values[i] = m.values[i];
        }
    }
    public Matrix(Vector<T>[] v, bool columnVectors = true)
    {
        //we obviosly require all v to be the same length

        dim = new int[2];
        if(columnVectors)
        {
            rows = v[0].rows;
            columns = v.Length;
            for (int j = 0; j < v.Length; j++) //iterate of "columns"
            {
                for (int i = 0; i < rows; i++)
                {
                    this[i, j] = v[j][i]; //the jth vector's ith component
                }
            }
        }
        else
        {
            rows = v.Length;
            columns = v[0].rows;
            for(int i=0; i<v.Length; i++)
            {
                for(int j=0; j<columns; j++)
                {
                    this[i, j] = v[i][j];
                }
            }
        }
    }

    public T this[int i, int j]
    {
        get { return values[(i * dim[1]) + j]; }
        set { values[(i * dim[1]) + j] = value; }
    }

    public int rows
    {
        get { return dim[0]; }
        set { dim[0] = value; }
    }

    public int columns
    {
        get { return dim[1]; }
        set { dim[1] = value; }
    }

    public static bool operator== (Matrix<T> a, Matrix<T> b)
    {
        return a.Equals(b);
    }

    public static bool operator!= (Matrix<T> a, Matrix<T> b)
    {
        return !a.Equals(b);
    }

    public static explicit operator T (Matrix<T> a)
    {
        if (a.rows != 1 || a.columns != 1) throw new System.Exception("a is not scalar dimension");
        return a.values[0];
    }

    public static Matrix<T> operator -(Matrix<T> a)
    {
        Matrix<T> m = new Matrix<T>(a.rows, a.columns);
        for(int i=0; i<a.rows*a.columns; i++)
        {
            m.values[i] = Operator<T>.negate(a.values[i]);
        }
        return m;
    }

    public static Matrix<T> operator+ (Matrix<T> a, Matrix<T> b)
    {
        if (a.rows != b.rows || a.columns != b.columns) throw new System.Exception("Matrix dimensions don't match");
        Matrix<T> m = new Matrix<T>(a.rows, a.columns);
        for (int i = 0; i < a.rows * a.columns; i++)
        {
            //dynamic _a = a.values[i]; //we need to convert to dynamic to be able to add
            //dynamic _b = b.values[i];
            //dynamic c = _a + _b;
            m.values[i] = Operator<T>.add(a.values[i], b.values[i]); //Add type T together
        }
        return m;
    }
    public static Matrix<T> operator -(Matrix<T> a, Matrix<T> b)
    {
        if (a.rows != b.rows || a.columns != b.columns) throw new System.Exception("Matrix dimensions don't match");
        Matrix<T> m = new Matrix<T>(a.rows, a.columns);
        for (int i = 0; i < a.rows * a.columns; i++)
        {
            //dynamic _a = a.values[i]; //we need to convert to dynamic to be able to add
            //dynamic _b = b.values[i];
            //dynamic c = _a + _b;
            m.values[i] = Operator<T>.substract(a.values[i], b.values[i]); //Add type T together
        }
        return m;
    }

    public static Matrix<T> operator*(Matrix<T> a, Matrix<T> b)
    {
        if (a.columns != b.rows) throw new System.Exception("Matrix dimensions are not compatible");
        Matrix<T> m = new Matrix<T>(a.rows, b.columns);
        //we do row times columns
        for (int i = 0; i < m.rows; i++)
        {
            for (int j = 0; j < m.columns; j++) //row 
            {
                for (int k = 0; k < a.columns; k++) //columns
                {
                    //m[i, j] += (dynamic)a[i, k] * (dynamic)b[k, j];
                    m[i, j] = Operator<T>.add(m[i, j], Operator<T>.multiply(a[i, k], b[k, j]));
                }
            }
        }
        return m;
    }

    public static Matrix<T> operator* (Matrix<T> a, T f)
    {
        Matrix<T> m = new Matrix<T>(a.rows, a.columns);
        for (int i = 0; i < a.rows * a.columns; i++)
        {
            m.values[i] = Operator<T>.multiply(a.values[i], f);
        }
        return m;
    }

    //for testing
    public static Matrix<T> Mult<U> (Matrix<U> a, Matrix<T> b) where U:new()
    {
        if (a.columns != b.rows) throw new System.Exception("Matrix dimensions are not compatible");
        Matrix<T> m = new Matrix<T>(a.rows, b.columns);
        //we do row times columns
        for (int i = 0; i < m.rows; i++)
        {
            for (int j = 0; j < m.columns; j++) //row 
            {
                for (int k = 0; k < a.columns; k++) //columns
                {
                    //m[i, j] += (dynamic)a[i, k] * (dynamic)b[k, j];
                    m[i, j] = Operator<T>.add(m[i, j], Operator<U,T>.multiply(b[k, j], a[i, k]));
                }
            }
        }
        return m;
    }
    public static Matrix<T> Mult<U>(Matrix<T> a, Matrix<U> b) where U : new()
    {
        if (a.columns != b.rows) throw new System.Exception("Matrix dimensions are not compatible");
        Matrix<T> m = new Matrix<T>(a.rows, b.columns);
        //we do row times columns
        for (int i = 0; i < m.rows; i++)
        {
            for (int j = 0; j < m.columns; j++) //row 
            {
                for (int k = 0; k < a.columns; k++) //columns
                {
                    //m[i, j] += (dynamic)a[i, k] * (dynamic)b[k, j];
                    m[i, j] = Operator<T>.add(m[i, j], Operator<U, T>.multiply(a[i, k], b[k, j]));
                }
            }
        }
        return m;
    }
    //end of testing

    public static Matrix<T> operator/(Matrix<T> a, T f) //only works if thing can be divided by float
    {
        return a * Operator<double, T>.divide(f, 1.0);
    }

    public static Matrix<T> operator^ (Matrix<T> a, uint p)
    {
        //we require p is positive
        if (p <= 0) return a;

        Matrix<T> m = a;
        for (uint i = 1; i < p; i++)
        {
            m *= a; //we multiply them one by one;
            //we could optimise it to power law
            //as p=p*p then do remainder
        }
        return m;
    }

    //public static Matrix<T> HadamartProduct(Matrix<T> a, Matrix<T> b)
    public static Matrix<T> EWiseProduct(Matrix<T> a, Matrix<T> b) //element wise product
    {
        if (a.rows != b.rows || a.columns != b.columns) throw new Exception("dimensions don't match, rows and columns need to match");
        Matrix<T> m = new Matrix<T>(a.rows, a.columns);
        for(int i=0; i<a.values.Length; i++)
        {
            m.values[i] = Operator<T>.multiply(a.values[i], b.values[i]);
        }
        return m;
    }
    //End of operators

    public static Matrix<T> Transpose(Matrix<T> a)
    {
        Matrix<T> m = new Matrix<T>(a.columns, a.rows);
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.columns; j++)
            {
                m[j, i] = a[i, j];
            }
        }
        return m;
    }

    public static Matrix<T> Kronecker(int row, int col)
    {
        return new Matrix<T>(row, col, true);
    }

    public static T Trace(Matrix<T> a)
    {
        if (a.rows != a.columns) throw new System.Exception("Invalid matrix, require square matrix");
        T trace = new T();
        for (int i = 0; i < a.rows; i++)
        {
            trace = Operator<T>.add(trace, a[i, i]);
        }
        return trace;
    }


    //inverse related stuff
    public static T Det(Matrix<T> a) //we require T to be able to add and multiply
    {
        if (a.columns != a.rows) throw new System.Exception("Can't find determinant, square Matrix required");
        if (a.rows == 1) return a[0, 0];
        //if (a.rows == 2) return Operator<T>.add(Operator<T>.multiply(a[0, 0],a[1, 1]), Operator<T>.negate(Operator<T>.multiply(a[1, 0],a[0, 1]))); //complicated
        if (a.rows == 2) return Operator<T>.substract(Operator<T>.multiply(a[0, 0], a[1, 1]), Operator<T>.multiply(a[1, 0], a[0, 1]));
        T sum = new T();
        Matrix<T> minor = new Matrix<T>(a.rows - 1, a.columns - 1);
        for (int k = 0; k < a.columns; k++)
        {
            for (int i = 0; i < a.rows - 1; i++)
            {
                int n = 0;
                for (int j = 0; j < a.columns; j++)
                {
                    if (j == k) continue;
                    minor[i, n] = a[i + 1, j];
                    n++;
                }
            }
            /*
            sum = Operator<T>.add(sum,
                Operator<T>.negate(
                    Operator<T>.multiply(a[0, k], Det(minor))
                    )
                    );*/ //sum += (-1f).Pow(k) * a[0, k] * Det(minor);
            sum = Operator<T>.substract(sum, Operator<T>.multiply(a[0, k], Det(minor)));
        }
        return sum;
    }

    public static Matrix<T> LowerTriangleInverse(Matrix<T> a)
    {
        if (a.rows != a.columns) throw new System.Exception("not square matrix");
        Matrix<T> inv = new Matrix<T>(a.rows, a.columns, true);
        int d = (int)Math.Ceiling(a.rows / (double)2);

        for (int k = 1; k < a.columns; k++)
        {
            for (int i = 0; i < a.rows - k; i++)
            {
                T sum = new T();
                for (int j = 0; j < i + k; j++)
                {
                    sum = Operator<T>.add(sum, Operator<T>.multiply(a[i + k, j], inv[j, i])); // sum = sum + (l[i + k, j] * invL[j, i]);
                }
                inv[i + k, i] = Operator<T>.negate(sum);
            }
        }

        return inv;
    }

    public static void LUDecomposition(Matrix<T> a, out Matrix<T> l, out Matrix<T> u)
    {
        //we create the matrices

        l = new Matrix<T>(a.rows, a.columns, true); //identity
        u = new Matrix<T>(a.rows, a.columns, a.values);

        for (int k = 0; k < a.rows - 1; k++)
        {
            //this is just simply taking away the above row from the previous one
            T max = u[k, k];
            for (int i = k + 1; i < a.rows; i++)
            {
                T c = Operator<T>.divide(u[i, k], max);
                l[i, k] = c;
                for (int j = k; j < a.columns; j++)
                {
                    if (j == k) { u[i, k] = new T(); }
                    else { u[i, j] = Operator<T>.substract(u[i, j], Operator<T>.multiply(c, u[k, j])); } //u[i, j] -= c * u[k, j];
                }
            }
        }
    }

    public static Matrix<T> Inverse(Matrix<T> a)
    {
        if (a.rows != a.columns) throw new System.Exception("Can't invert non-square matrix");
        if (a.rows == 1) return new Matrix<T>(1, 1, new T[] { Operator<T>.divide(Operator<int, T>.convert(1), a[0, 0]) }); //trying to do 1/a[0,0]
        if (a.rows == 2) //explicit inverse
        {
            //T det = Det(a);
            T invDet = Operator<T>.divide(Operator<int, T>.convert(1), Det(a));
            T[] v = {   a[1,1], Operator<T>.negate(a[0,1]),
                        Operator<T>.negate(a[1,0]), a[0,0]};
            Matrix<T> m = new Matrix<T>(2, 2, v);
            return m * invDet;
        }
        //throw new System.Exception("not implemented");

        Matrix<T> l = new Matrix<T>(0, 0);
        Matrix<T> u = new Matrix<T>(0, 0);
        LUDecomposition(a, out l, out u);

        Matrix<T> d = new Matrix<T>(a.rows, a.rows);
        for (int i = 0; i < a.rows; i++)
        {
            d[i, i] = Operator<T>.divide(Operator<int, T>.convert(1), u[i, i]);
        }
        u = Transpose(d * u);
        u = Transpose(LowerTriangleInverse(u)) * d;
        l = LowerTriangleInverse(l);

        return u * l;
    }
    //end of inverse stuff


    //NORMS
    public static T OneNorm(Matrix<T> a) //columns sum
    {
        //column sum
        T max = new T();
        for (int j = 0; j < a.columns; j++)
        {
            T sum = new T();
            for (int i = 0; i < a.rows; i++)
            {
                sum = Operator<T>.add(sum, GenericMath.Abs(a[i, j])); //sum += abs(a[i,j])
            }
            if (j == 0) { max = sum; continue; }
            if (Operator<T>.greaterThan(sum, max)) max = sum;
        }
        return max;
        //throw new SystemException("Not implemented YET");
    }

    public static T InfNorm(Matrix<T> a) //infinity norm
    {
        //row sum
        T max = new T();
        for (int i = 0; i < a.columns; i++)
        {
            T sum = new T();
            for (int j = 0; j < a.columns; j++)
            {
                sum = Operator<T>.add(sum, GenericMath.Abs(a[i, j]));
            }
            if (i == 0) { max = sum; continue; }
            if (Operator<T>.greaterThan(sum, max)) max = sum;
        }
        return max;
    }

    public bool Equals(Matrix<T> other)
    {
        //throw new NotImplementedException();
        if (other.rows != rows) return false; //same dimensions
        if (other.columns != columns) return false;
        if (!(values[0].GetType() is T)) return false; //to make sure they're the same type

        //compare values
        for (int i = 0; i < values.Length; i++)
        {
            if (Operator<T>.notEqual(values[i], other.values[i])) return false;
        }
        return true; //should only get here if every values is the same
    }

    public override bool Equals(object obj) //this is for the comparison stuff
    {
        if (obj == null || !(obj is Matrix<T>)) return false;
        //return base.Equals(obj);
        return this.Equals((Matrix<T>)obj);
    }

    public string ToString(string format, IFormatProvider formatProvider)
    {
        //throw new NotImplementedException();
        string str = "";

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                str += this[i, j].ToString();
                if (j < columns - 1) str += ", ";
            }
            str += "@";
        }
        //str = str.Replace("@", Environment.NewLine);
        return str.Replace("@", Environment.NewLine);
    }

    public override string ToString()
    {
        //ToString("RC", );
        return ToString("RC", System.Globalization.CultureInfo.CurrentCulture); //we need the long thing just for this
    }

    public void Dispose()
    {
        //throw new NotImplementedException();
        values = null; //for garbage collector
        rows = new int(); //for garbage collector
        columns = new int(); //for garbage collector
    }
}
