using System;
using PUtil.Operators;
using PUtil.GenericMath;

public class Complex<T> : Vector<T>, IEquatable<Complex<T>>, IFormattable, IDisposable where T : new() //matrix is a rank 2 tensor
{
    //complex numbers
    public Complex()
    {
        dim = new int[2];
        dim[0] = 2;
        dim[1] = 1;
        values = new T[2];
        values[0] = new T();
        values[1] = new T();
    }
    public Complex(T val)
    {
        dim = new int[2];
        dim[0] = 2;
        dim[1] = 1;
        values = new T[2];
        values[0] = val;
        values[1] = new T();
    }

    public Complex(T x, T y)
    {
        dim = new int[2];
        dim[0] = 2;
        dim[1] = 1;
        values = new T[2];
        values[0] = x;
        values[1] = y;
    }
    public Complex(Vector<T> v)
    {
        if (v.rows != 2) throw new Exception("Vector is not dimension 2!");
        dim = new int[2];
        dim[0] = 2;
        dim[1] = 1;
        values = new T[2];
        values[0] = v.values[0];
        values[1] = v.values[1];
    }

    public T real //real part //could be x
    {
        get { return values[0]; }
        set { values[0] = value; }
    }
    public T imaginary //imaginary part //could be y
    {
        get { return values[1]; }
        set { values[1] = value; }
    }
    public Complex<T> conjugate
    {
        get { return new Complex<T>(real, Operator<T>.negate(imaginary)); }
    }

    public bool Equals(Complex<T> other)
    {
        if (Operator<T>.notEqual(real, other.real)) return false;
        if (Operator<T>.notEqual(imaginary, other.imaginary)) return false;
        return true;
    }

    public override bool Equals(object obj)
    {
        //throw new NotImplementedException();
        if (obj == null || !(obj is Complex<T>)) return false;
        //return base.Equals(obj);
        return this.Equals((Complex<T>)obj);
    }

    public static bool operator ==(Complex<T> a, Complex<T> b)
    {
        return a.Equals(b);
    }
    public static bool operator !=(Complex<T> a, Complex<T> b)
    {
        return !a.Equals(b);
    }

    public static Complex<T> operator+ (Complex<T> a, Complex<T> b)
    {
        return new Complex<T>(Operator<T>.add(a.real, b.real), Operator<T>.add(a.imaginary, b.imaginary));
    }
    public static Complex<T> operator +(Complex<T> a, T b)
    {
        return new Complex<T>(Operator<T>.add(a.real, b), a.imaginary);
    }

    public static Complex<T> operator -(Complex<T> a)
    {
        return new Complex<T>(Operator<T>.negate(a.real), Operator<T>.negate(a.imaginary));
    }

    public static Complex<T> operator* (Complex<T> a, Complex<T> b)
    {
        T x = Operator<T>.substract(Operator<T>.multiply(a.real, b.real), Operator<T>.multiply(a.imaginary, b.imaginary));
        T y = Operator<T>.add(Operator<T>.multiply(a.real, b.imaginary), Operator<T>.multiply(b.real, a.imaginary));
        return new Complex<T>(x, y);
    }

    public static Complex<T> operator *(Complex<T> a, T f)
    {
        return new Complex<T>(Operator<T>.multiply(a.real, f), Operator<T>.multiply(a.imaginary, f));
    }
    public static Complex<T> operator /(Complex<T> a, T f)
    {
        return new Complex<T>(Operator<T>.divide(a.real, f), Operator<T>.divide(a.imaginary, f));
    }

    //Formatting
    new public string ToString(string format, IFormatProvider formatProvider)
    {
        //throw new NotImplementedException();
        return (real + " + i" + imaginary);
    }

    public override string ToString()
    {
        //ToString("RC", );
        return ToString("RC", System.Globalization.CultureInfo.CurrentCulture); //we need the long thing just for this
    }
}
