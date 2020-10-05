using System;
//using System.Dynamic;
//using System.Linq.Expressions;
//using System.Linq;
using PUtil.Operators;
//using PUtil.Interfaces;
using PUtil.GenericMath;

[System.Serializable]
public class Vector<T> : Matrix<T>, IEquatable<Vector<T>>, IFormattable, IDisposable where T : new() //matrix is a rank 2 tensor
{
    public Vector()
    {
        dim = new int[2];
        dim[0] = 1;
        dim[1] = 1;
        values = new T[1];
        values[0] = new T();
    }

    public Vector(int row, T[] vals)
    {
        dim = new int[2];
        dim[0] = row;
        dim[1] = 1;
        values = new T[row];
        for(int i=0; i<row; i++)
        {
            values[i] = vals[i];
        }
    }

    public Vector(int row, bool one=false)
    {
        dim = new int[2];
        dim[0] = row;
        dim[1] = 1;
        values = new T[row];
        if(one)
        {
            for(int i=0; i<row; i++)
            {
                values[i] = Operator<int, T>.convert(1);
            }
        }
    }

    public Vector(Matrix<T> m, bool truncate = false)
    {
        if (!truncate)
        {
            if (m.columns != 1) throw new System.Exception("Matrix dimension doesn't match"); //maybe we could use the first row of values then we don't need that it should match
        }

        dim = new int[2];
        dim[0] = m.rows;
        dim[1] = 1;
        values = new T[m.rows];
        for(int i=0; i<m.rows; i++)
        {
            values[i] = m[i,0];
        }
    }

    public T this[int i]
    {
        get { return values[i]; }
        set { values[i] = value; }
    }

    public bool Equals(Vector<T> other)
    {
        //throw new NotImplementedException();
        if (rows != other.rows) return false;
        for(int i=0; i<rows; i++)
        {
            if (Operator<T>.notEqual(values[i], other.values[i])) return false;
        }
        return true;

    }
    public override bool Equals(object obj)
    {
        //throw new NotImplementedException();
        if (obj == null || !(obj is Vector<T>)) return false;
        //return base.Equals(obj);
        return this.Equals((Vector<T>)obj);
    }

    public static bool operator ==(Vector<T> a, Vector<T> b)
    {
        return a.Equals(b);
    }
    public static bool operator !=(Vector<T> a, Vector<T> b)
    {
        return !a.Equals(b);
    }

    public static Vector<T> operator -(Vector<T> a)
    {
        Vector<T> v = new Vector<T>(a.rows);
        for (int i = 0; i < a.rows; i++)
        {
            v[i] = Operator<T>.negate(a[i]);
        }
        return v;
    }

    public static Vector<T> operator+ (Vector<T> a, Vector<T> b)
    {
        Vector<T> c = new Vector<T>(a.rows);
        for(int i=0; i<a.rows; i++)
        {
            c[i] = Operator<T>.add(a[i], b[i]);
        }
        return c;
    }

    public static Vector<T> operator -(Vector<T> a, Vector<T> b)
    {
        Vector<T> c = new Vector<T>(a.rows);
        for (int i = 0; i < a.rows; i++)
        {
            c[i] = Operator<T>.substract(a[i], b[i]);
        }
        return c;
    }

    public static Vector<T> operator* (Vector<T> a, T f)
    {
        Vector<T> c = new Vector<T>(a.rows);
        for (int i = 0; i < a.rows; i++)
        {
            c[i] = Operator<T>.multiply(a[i], f);
        }
        return c;
    }
    public static Vector<T> operator* (T f, Vector<T> a)
    {
        return a * f;
    }

    public static Vector<T> operator/ (Vector<T> a, T f)
    {
        Vector<T> c = new Vector<T>(a.rows);
        for (int i = 0; i < a.rows; i++)
        {
            c[i] = Operator<T>.divide(a[i], f);
        }
        return c;
    }

    public static Vector<T> operator* (Matrix<T> m, Vector<T> a)
    {
        if (m.columns != a.rows) throw new System.Exception("Matrix dimension doesn't match");
        return new Vector<T>(m * (Matrix<T>)a);
        //Vector<T> c = new Vector<T>(m * a);
        //return c;
    }

    /*public static explicit operator Vector<T>(Matrix<T> a)
    {
        if (a.columns != 1) throw new System.Exception("a must have only 1 column");
        return new Vector<T>(a.rows, a.values);
    }*/

    public static T Dot(Vector<T> a, Vector<T> b) // transpose(a) * b
    {
        //this should give us a 1,1 matrix
        return (T)Operator<Matrix<T>>.multiply(Transpose(a), b);
    }

    public static Vector<T> EWiseProduct(Vector<T> a, Vector<T> b)
    {
        return new Vector<T>(Matrix<T>.EWiseProduct((Matrix<T>)a, (Matrix<T>)b));
    }

    public static Matrix<T> OuterProduct(Vector<T> a, Vector<T> b)
    {
        Matrix<T> m = new Matrix<T>(a.rows, b.rows); //each row of the matrix corresponds to a row in a
        for(int i=0; i<a.rows; i++)
        {
            for(int j=0; j<b.rows; j++)
            {
                m[i, j] = Operator<T>.multiply(a[i], b[j]); //m[i,j] = a[i] * b[j]
            }
        }
        return m;
    }

    public static Vector<T> Abs(Vector<T> x)
    {
        Vector<T> y = new Vector<T>(x.rows);
        for(int i=0; i<x.rows; i++)
        {
            x[i] = Operator<T>.greaterThanOrEqual(x[i], new T()) ? x[i] : Operator<T>.negate(x[i]);
        }
        return y;
    }

    public static T[] ToArray(Vector<T> x)
    {
        return x.values;
    }
    public T[] ToArray()
    {
        return values;
    }
}
