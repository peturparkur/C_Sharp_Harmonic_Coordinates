using System;
//using System.Dynamic;
using System.Linq.Expressions;
//using System.Linq;
using PUtil.Operators;
using PUtil.Interfaces;
using PUtil.GenericMath;

//for testing purposes only
using Unity.Entities;
using Unity.Collections;

public struct Matrix1<T> : IEquatable<Matrix1<T>>, IFormattable, IDisposable where T:new()
{
    public T[] values;
    public int rows;
    public int columns;

    public Matrix1(int row, int col, T[] vals)
    {
        if (row <= 0 || col <= 0) throw new System.Exception("Can't have negative rows or columns");
        if (vals.Length != row * col) throw new System.Exception("Length of array doesn't match dimensions");
        rows = row;
        columns = col;
        values = new T[row * columns];
        for (int i = 0; i < vals.Length; i++)
        {
            values[i] = vals[i];
        }
    }

    public Matrix1(int row, int col, bool identity = false)
    {
        if (row <= 0 || col <= 0) throw new System.Exception("Can't have negative rows or columns");
        rows = row;
        columns = col;
        values = new T[row * col];
        if(identity)
        {
            int n = 0;
            int m = 0;
            if (rows >= columns) { n = rows; m = columns; }
            else { n = columns; m = rows; }

            for(int i=0; i<n; i++)
            {
                values[i*m + i] = Operator<int,T>.convert(1); //we try to fill it up with ones
            }
        }
    }

    public T this[int i, int j]
    {
        get { return values[(i * columns) + j]; }
        set { values[(i * columns) + j] = value; }
    }

    public static bool operator ==(Matrix1<T> a, Matrix1<T> b)
    {
        return a.Equals(b);
    }

    public static bool operator !=(Matrix1<T> a, Matrix1<T> b)
    {
        return !a.Equals(b);
    }
    
    public static Matrix1<T> Add(Matrix1<T> a, Matrix1<T> b)
    {
        if (a.rows != b.rows || a.columns != b.columns) throw new System.Exception("Matrix dimensions don't match");
        Matrix1<T> m = new Matrix1<T>(a.rows, a.columns);
        for (int i = 0; i < a.rows * a.columns; i++)
        {
            //dynamic _a = a.values[i]; //we need to convert to dynamic to be able to add
            //dynamic _b = b.values[i];
            //dynamic c = _a + _b;
            m.values[i] = Operator<T>.add(a.values[i], b.values[i]); //Add type T together
        }

        return m;
    }
    /*public static bool Equal<T>(T a, T b) where T : class
    {
        return a==b;
    }*/

    public static Matrix1<T> operator+ (Matrix1<T> a, Matrix1<T> b)
    {
        //return Add(a, b); //don't know if this is good practice
        return Add(a, b);
    }

    public static Matrix1<T> operator* (Matrix1<T> a, T f)
    {
        return Multiply(a, f);
    }

    /*public static Matrix<T> operator*= (Matrix<T> b)
    {
        return Multiply(this, b);
    }*/

    public static Matrix1<T> operator/ (Matrix1<T> a, T f) //only works if thing can be divided by float
    {
        return Multiply(a, Operator<double, T>.divide(f, 1.0));
    }

    public static Matrix1<T> operator* (Matrix1<T> a, Matrix1<T> b)
    {
        return Multiply(a, b);
    }

    public static Matrix1<T> Multiply(Matrix1<T> a, Matrix1<T> b) //Matrix multiplication
    {
        if (a.rows != b.columns) throw new System.Exception("Matrix dimensions are not compatible");
        Matrix1<T> m = new Matrix1<T>(a.rows, a.columns);
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

    public static Matrix1<T> Multiply(Matrix1<T> a, T f)
    {
        Matrix1<T> m = new Matrix1<T>(a.rows, a.columns);
        for(int i=0; i<a.rows*a.columns; i++)
        {
            m.values[i] = Operator<T>.multiply(a.values[i], f);
        }
        return m;
    }

    public static Matrix1<T> Power(Matrix1<T> a, uint p)
    {
        //we require p is positive
        if (p <= 0) return a;

        Matrix1<T> m = a;
        for(uint i=1; i<p; i++)
        {
            m *= a; //we multiply them one by one;
            //we could optimise it to power law
            //as p=p*p then do remainder
        }
        return m;
    }

    public static Matrix1<T> Transpose(Matrix1<T> a)
    {
        Matrix1<T> m = new Matrix1<T>(a.columns, a.rows);
        for(int i=0; i<a.rows; i++)
        {
            for(int j=0; j<a.columns; j++)
            {
                m[i, j] = a[j, i];
            }
        }
        return m;
    }

    public static Matrix1<T> kronecker(int row, int col)
    {
        return new Matrix1<T>(row, col, true);
    }

    public static T Trace(Matrix1<T> a)
    {
        if (a.rows != a.columns) throw new System.Exception("Invalid matrix, require square matrix");
        T trace = new T();
        for(int i=0; i<a.rows; i++)
        {
            trace = Operator<T>.add(trace, a[i, i]);
        }
        return trace;
    }

    public static T Det(Matrix1<T> a) //we require T to be able to add and multiply
    {
        if (a.columns != a.rows) throw new System.Exception("Can't find determinant, square Matrix required");
        if (a.rows == 1) return a[0, 0];
        //if (a.rows == 2) return Operator<T>.add(Operator<T>.multiply(a[0, 0],a[1, 1]), Operator<T>.negate(Operator<T>.multiply(a[1, 0],a[0, 1]))); //complicated
        if (a.rows == 2) return Operator<T>.substract(Operator<T>.multiply(a[0, 0], a[1, 1]), Operator<T>.multiply(a[1, 0], a[0, 1]));
        T sum = new T();
        Matrix1<T> minor = new Matrix1<T>(a.rows - 1, a.columns - 1);
        for(int k=0; k<a.columns; k++)
        {
            for(int i=0; i<a.rows-1; i++)
            {
                int n = 0;
                for(int j=0; j<a.columns; j++)
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

    public static Matrix1<T> Inverse(Matrix1<T> a)
    {
        if (a.rows != a.columns) throw new System.Exception("Can't invert non-square matrix");
        if (a.rows == 1) return new Matrix1<T>(1, 1, new T[] { Operator<T>.divide(Operator<int, T>.convert(1), a[0, 0]) }); //trying to do 1/a[0,0]
        if(a.rows == 2) //explicit inverse
        {
            //T det = Det(a);
            T invDet = Operator<T>.divide(Operator<int, T>.convert(1), Det(a));
            Matrix1<T> m = new Matrix1<T>(2 , 2, new T[]
                {
                    a[1,1], Operator<T>.negate(a[0,1]),
                    Operator<T>.negate(a[1,0]), a[0,0]
                });
            return m*invDet;
        }
        //throw new System.Exception("not implemented");

        Matrix1<T> l = new Matrix1<T>(0, 0);
        Matrix1<T> u = new Matrix1<T>(0, 0);
        LUDecomposition(a, out l, out u);

        Matrix1<T> d = new Matrix1<T>(a.rows, a.rows);
        for(int i=0; i<a.rows; i++)
        {
            d[i, i] = Operator<T>.divide(Operator<int, T>.convert(1), u[i, i]);
        }
        u = Transpose(d * u);
        u = Transpose(LowerTriangleInverse(u)) * d;
        l = LowerTriangleInverse(l);

        return u * l;
    }

    public static T OneNorm(Matrix1<T> a) //columns sum
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
            if(j==0) { max = sum; continue; }
            if (Operator<T>.greaterThan(sum, max)) max = sum;
        }
        return max;
        //throw new SystemException("Not implemented YET");
    }

    public static T InfNorm(Matrix1<T> a) //infinity norm
    {
        //row sum
        T max = new T();
        for(int i=0; i<a.columns; i++)
        {
            T sum = new T();
            for(int j=0; j<a.columns; j++)
            {
                sum = Operator<T>.add(sum, GenericMath.Abs(a[i, j]));
            }
            if (i == 0) { max = sum; continue; }
            if (Operator<T>.greaterThan(sum, max)) max = sum;
        }
        return max;
    }

    public static void LUDecomposition(Matrix1<T> a, out Matrix1<T> l, out Matrix1<T> u)
    {
        //we create the matrices

        l = new Matrix1<T>(a.rows, a.columns, true); //identity
        u = new Matrix1<T>(a.rows, a.columns, a.values);

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

    public static Matrix1<T> LowerTriangleInverse(Matrix1<T> a)
    {
        if (a.rows != a.columns) throw new System.Exception("not square matrix");
        Matrix1<T> inv = new Matrix1<T>(a.rows, a.columns, true);
        int d = (int)Math.Ceiling(a.rows / (double)2);

        for(int k=1; k<a.columns; k++)
        {
            for(int i=0; i<a.rows - k; i++)
            {
                T sum = new T();
                for(int j=0; j<i+k; j++)
                {
                    sum = Operator<T>.add(sum, Operator<T>.multiply(a[i + k, j], inv[j, i])); // sum + l[i + k, j] * invL[j, i];
                }
                inv[i + k, i] = Operator<T>.negate(sum);
            }
        }

        return inv;
    }

    public bool Equals(Matrix1<T> other) //for IEquatable
    {
        //throw new NotImplementedException();
        if (other.rows != rows) return false; //same dimensions
        if (other.columns != columns) return false;
        if (!(values[0].GetType() is T)) return false; //to make sure they're the same type

        //compare values
        for (int i=0; i<values.Length; i++)
        {
            if (Operator<T>.notEqual(values[i], other.values[i])) return false;
        }
        return true; //should only get here if every values is the same
    }

    public override bool Equals(object obj) //this is for the comparison stuff
    {
        if (obj == null || !(obj is Matrix1<T>)) return false;
        //return base.Equals(obj);
        return this.Equals(obj);
    }

    public override string ToString()
    {
        //ToString("RC", );
        return ToString("RC", System.Globalization.CultureInfo.CurrentCulture); //we need the long thing just for this
    }
    
    public string ToString(string format, IFormatProvider formatProvider)
    {
        //throw new NotImplementedException();
        string str = "";

        for(int i=0; i<rows; i++)
        {
            for(int j=0; j<columns; j++)
            {
                str += this[i,j].ToString();
                if (j < columns - 1) str += ", ";
            }
            str += "@";
        }
        //str = str.Replace("@", Environment.NewLine);
        return str.Replace("@", Environment.NewLine);
    }

    public void Dispose() //For IDisposable
    {
        //throw new NotImplementedException();
        values = null; //for garbage collector
        rows = new int(); //for garbage collector
        columns = new int(); //for garbage collector
    }
}

//Testing my own interface for arithmetic operators
public interface IArithmetic<T>
{
    T Addition(T a, T b);
    T Substraction(T a, T b);
    T Multiplication(T a, T b);
    T Division(T a, T b);

    T ScalarMultiplication<U>(T a, U f);
    T ScalarDivision<U>(T a, U f);
}

/* //This has been moved to Putil_Operators, Accessed using Putil.Operators namespace
public static class PMathUtil //Just for implementing the generic addition and multiplication methods
{

    //public static Func<Tfrom, Tto> convert<Tfrom, Tto>() = CreateExpression<Tfrom, Tto>(body => Expression.Convert(body, typeof(Tto)));
    
    public static T Negate<T>(T a)
    {
        ParameterExpression paramA = Expression.Parameter(typeof(T), "a");
        UnaryExpression body = Expression.Negate(paramA);
        Func<T, T> negate = Expression.Lambda<Func<T, T>>(body, paramA).Compile();
        return negate(a);
    }

    //These obviously return errors if it isn't defined for the specified type
    public static T Add<T>(T a, T b) //Creating an addition of type T together
    {
        ParameterExpression paramA = Expression.Parameter(typeof(T), "a");
        ParameterExpression paramB = Expression.Parameter(typeof(T), "b");
        BinaryExpression body = Expression.Add(paramA, paramB);
        //compilation
        //Func<Tin, Tin2, Tout>
        Func<T, T, T> add = Expression.Lambda<Func<T, T, T>>(body, paramA, paramB).Compile();
        return add(a, b);
    }

    public static T Multiply<T>(T a, T b) //Creating a Multiplication of type T together
    {
        ParameterExpression paramA = Expression.Parameter(typeof(T), "a");
        ParameterExpression paramB = Expression.Parameter(typeof(T), "b");
        BinaryExpression body = Expression.Multiply(paramA, paramB);
        //compilation
        //Func<Tin, Tin2, Tout>
        Func<T, T, T> multiply = Expression.Lambda<Func<T, T, T>>(body, paramA, paramB).Compile();
        return multiply(a, b);
    }

    public static T Multiply2<T, U>(T a, U b) //Creating an multiplication of type T and U together
    {
        ParameterExpression paramA = Expression.Parameter(typeof(T), "a");
        ParameterExpression paramB = Expression.Parameter(typeof(U), "b");
        BinaryExpression body = Expression.Multiply(paramA, paramB);
        //compilation
        //Func<Tin, Tin2, Tout>
        Func<T, U, T> multiply = Expression.Lambda<Func<T, U, T>>(body, paramA, paramB).Compile();
        return multiply(a, b);
    }

    public static T Divide<T>(T a, T b) //Creating a Multiplication of type T together
    {
        ParameterExpression paramA = Expression.Parameter(typeof(T), "a");
        ParameterExpression paramB = Expression.Parameter(typeof(T), "b");
        BinaryExpression body = Expression.Divide(paramA, paramB);
        //compilation
        //Func<Tin, Tin2, Tout>
        Func<T, T, T> divide = Expression.Lambda<Func<T, T, T>>(body, paramA, paramB).Compile();
        return divide(a, b);
    }

    public static T Divide<T,U>(U a, T b) //Creating a Multiplication of type T together
    {
        ParameterExpression paramA = Expression.Parameter(typeof(U), "a");
        ParameterExpression paramB = Expression.Parameter(typeof(T), "b");
        BinaryExpression body = Expression.Divide(paramA, paramB);
        //compilation
        //Func<Tin, Tin2, Tout>
        Func<U, T, T> divide = Expression.Lambda<Func<U, T, T>>(body, paramA, paramB).Compile();
        return divide(a, b);
    }
}
*/
public unsafe struct MatrixNxN : IEquatable<MatrixNxN>, IFormattable //unsafe is when we don't check memory space, can use pointers
{
    int row;
    int column;
    public fixed float val[9]; //fixed array, *pointer --> it's unmanaged code

    public bool Equals(MatrixNxN other)
    {
        throw new NotImplementedException();
    }

    public string ToString(string format, IFormatProvider formatProvider)
    {
        throw new NotImplementedException();
    }
}

public struct MatrixComponent : IComponentData
{
    MatrixNxN matrix;
}
