//using System.Collections;
//using System.Collections.Generic;
using UnityEngine;
using PUtil.Operators;
using PUtil.Interfaces;

[System.Serializable]
public class Tensor<T> //: IVectorSpace<Tensor<T>, T>
{
    public T[] values;
    public int[] dim; //length of each dimension
               //int rank; //how many dimension

    int Rank
    {
        get { return dim.Length; }
    }

    public Tensor() //we create a scalar
    {
        dim = new int[1];
        dim[0] = 1; //scalar
        values = new T[1];
    }

    public Tensor(Tensor<T> a)
    {
        dim = new int[a.dim.Length];
        for (int i = 0; i < a.dim.Length; i++)
        {
            dim[i] = a.dim[i];
            //multiple *= dims[i];
        }
        //length of vals should be the multiple of dims
        values = new T[a.values.Length];
        for (int i = 0; i < a.values.Length; i++)
        {
            values[i] = a.values[i];
        }
    }

    public Tensor(T[] vals, int[] dims)
    {
        dim = new int[dims.Length];
        //int multiple = 1;
        for (int i = 0; i < dims.Length; i++)
        {
            dim[i] = dims[i];
            //multiple *= dims[i];
        }
        //length of vals should be the multiple of dims
        values = new T[vals.Length];
        for (int i = 0; i < vals.Length; i++)
        {
            values[i] = vals[i];
        }
    }
    public Tensor(int[] dims)
    {
        dim = new int[dims.Length];
        int multiple = 1;
        for(int i=0; i<dims.Length; i++)
        {
            dim[i] = dims[i];
            multiple *= dims[i];
        }
        values = new T[multiple]; //created enough empty things
    }

    public static bool CheckValid(Tensor<T> x)
    {
        int multiple = 1;
        for(int i=0; i<x.dim.Length; i++)
        {
            multiple *= x.dim[i];
        }

        if (x.values.Length != multiple) return false; //we only require that the length matches
        return true;
    }

    public int Index(int[] a)
    {
        int p = a[a.Length - 1]; //last value of A
        for (int i = a.Length - 2; i >= 0; i--)
        {
            int m = a[i];
            for (int j = i + 1; j < a.Length; j++)
            {
                m *= dim[j];
            }
            p += m;
        }
        return p;
    }
    public int[] IndexToAccessor(int index)
    {
        int[] access = new int[dim.Length];
        for(int i=dim.Length-1; i>=0; i--)
        {
            int mult = 1;
            for(int j=1; j<=i; j++)
            {
                mult *= dim[dim.Length-j];
            }
            int div = index / mult; //this should round down as far as I am aware
            access[dim.Length-i-1] = div;
            index -= div * mult; //updated index
        }
        return access;
    }

    public T this[int[] a]
    {
        get
        {
            int p = a[a.Length-1]; //last value of A
            for(int i=a.Length-2; i>=0; i--)
            {
                int m = a[i];
                for(int j=i+1; j<a.Length; j++)
                {
                    m *= dim[j];
                }
                p += m;
            }
            return values[p];
        }

        set
        {
            int p = a[a.Length - 1]; //last value of A
            for (int i = a.Length - 2; i >= 0; i--)
            {
                int m = a[i];
                for (int j = i+1; j < a.Length; j++)
                {
                    m *= dim[j];
                }
                p += m;
            }
            values[p] = value;
        }
    }

    public static Tensor<T> operator+ (Tensor<T> a, Tensor<T> b)
    {
        Tensor<T> c = new Tensor<T>(a.dim);
        for (int i = 0; i < c.values.Length; i++)
        {
            c.values[i] = Operator<T>.add(a.values[i], b.values[i]);
        }
        return c;
    }


    public static Tensor<T> ScalarMultiply(Tensor<T> a, T s)
    {
        throw new System.NotImplementedException();
    }

    public static T ScalarIdentity()
    {
        throw new System.NotImplementedException();
    }

    public static Tensor<T> Add(Tensor<T> a, Tensor<T> b)
    {
        Tensor<T> c = new Tensor<T>(a.dim);
        for (int i = 0; i < c.values.Length; i++)
        {
            c.values[i] = Operator<T>.add(a.values[i], b.values[i]);
        }
        return c;
    }

    /*public Tensor<T> AdditiveIdentity() //should be tensor of zeros but we need some dimension
    {

    }*/

    public static Tensor<T> Negate(Tensor<T> a)
    {
        Tensor<T> c = new Tensor<T>(a.dim);
        for(int i=0; i<c.values.Length; i++)
        {
            //c.values[i] = Operator<T>.negate(a.values[i]);
        }
        return c;
    }

    public static Tensor<T> Substract(Tensor<T> a, Tensor<T> b)
    {
        Tensor<T> c = new Tensor<T>(a.dim);
        for (int i = 0; i < c.values.Length; i++)
        {
            c.values[i] = Operator<T>.substract(a.values[i], b.values[i]);
        }
        return c;
    }

    public static double Length(Tensor<T> a) //not sure how to do tensor norm
    {
        throw new System.NotImplementedException();
    }

    public static double Distance(Tensor<T> a, Tensor<T> b)
    {
        throw new System.NotImplementedException();
    }
}
