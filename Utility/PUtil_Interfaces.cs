using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Linq.Expressions;
//using System.Text;
//using System.Threading.Tasks;

namespace PUtil.Interfaces //useful for guidelines for creating new classes with certain structure
{
    ///<summary>
    /// IMetricSpace defines classes where we have a function for distance between two points, where we don't need an origin
    ///</summary>
    public interface IMetricSpace<T> //defines something that can measure distances
    {
        double Distance(T a, T b); //this implements the abstract definition of length
    }

    ///<summary>
    /// INormSpace defines classes where we have a function for distance between two points, where we have an origin
    /// This also provides the Distance from IMetricSpace
    ///</summary>
    public interface INormSpace<T> : IMetricSpace<T> //defines something that can define distance from origin
    //we require that the origin exists
    {
        //this can be applied to vector spaces
        double Length(T a); //distance from the origin
        //here
        //double Distance(T a, T b); //==> should be Length(a-b)
    }

    /*
    public interface IVectorNorm<T, U> : INormSpace<T> , IVectorSpace<T,U>
    {
        double Norm(T a); //this should be Length(a);
    }
    */

    ///<summary>
    /// Group defines classes where we have an operation (eg.: Add) which has an identity (eg.: 0) and has an inverse (eg.: -x), 
    /// We can use the inverse to have Substraction operator
    ///</summary>
    public interface IGroup<T> //we define the operation as addition but could be multiplication
    {
        //Group is defined as 
        //-the Set T is closed under the operation
        //-The operation is associative ==> (a*b)*c = a*(b*c)
        //-Has identity e ==> a*e = e*a = a
        //Every element has inverse ==> a*a^(-1) = a^(-1)*a = identity e
        T Add(T a, T b); //group's operation
        //T AdditiveIdentity(); //return the identity
        T Negate(T a); // returns the inverse element of the operation, this technically doesn't need to be unique
        T Substract(T a, T b); //==> return a + Negate(b);
    }

    ///<summary>
    /// Ring which are groups, and they have another set of operations (eg.: multiplication), 
    /// multiplication is only required to have an inverse (eg.: 1)
    ///</summary>
    public interface IRing<T> : IGroup<T>
    {
        //defines a Set T that is a ring
        //it needs the following
        //-Addition + is
        //  -Associative ==> (a+b)+c = a+(b+c)
        //  -Commutative ==> a+b = b+a
        //  -Has Identity
        //  -Has Inverse
        //-Multiplication * is
        //  -Associative ==> same as above
        //  -Has Identity
        //-Addition and multiplication is distributive ==> a*(b+c) = a*b + a*c or (b+c)*a = b*a + c*a

        //Addition stuff
        //T Add(T a, T b);
        //T Negate(T a); //inverse of a
        //T AdditiveIdentity(); //return identity
        //Substraction is just adding the negation of b ==>
        //T Substract(T a, T b); //==> return a + Negate(b);

        //Multiplication stuff
        T Multiply(T a, T b);
        //T MultiplicativeIdentity();

    }

    ///<summary>
    /// Field which are Rings, where multiplication has an inverse (eg.: 1/x), so we have a division operator.
    ///</summary>
    public interface IField<T> : IRing<T>
    {
        //Same as ring but it also has Inverse for multiplication
        //and maybe Multiplication is commutative ==> a*b = b*a
        T Inverse(T a); //return multiplicative inverse
        T Divide(T a, T b); //==> return a * Inverse(b);
    }

    ///<summary>
    /// VectorSpace defines classes which are Groups with respect to themselves and have an operation defined with respect to another Type T (eg.: float * Vector3(float)),
    /// this operation returns the Type of the group, and has an identity with respect to it. It also defines a NormSpace (Length).
    ///</summary>
    public interface IVectorSpace<T, U> : IGroup<T>, INormSpace<T>
        //Every vector space has the maximum norm defined at least
    {
        //to be linear we require that the set T is closed under addition
        //this is satisfied if T is a group
        //We require that T is linear relative to U

        T ScalarMultiply(T a, U s);
        U ScalarIdentity(); //return the identity where a*s = s*a = a
    }

    /*
     //Just for testing
    public class Vectors<U> : IField<Vectors<U>>, IVectorSpace<Vectors<U>, U>
    {
        public Vectors<U> Add(Vectors<U> a, Vectors<U> b)
        {
            throw new NotImplementedException();
        }

        public Vectors<U> AdditiveIdentity()
        {
            throw new NotImplementedException();
        }

        public double Distance(Vectors<U> a, Vectors<U> b)
        {
            throw new NotImplementedException();
        }

        public Vectors<U> Divide(Vectors<U> a, Vectors<U> b)
        {
            throw new NotImplementedException();
        }

        public Vectors<U> Inverse(Vectors<U> a)
        {
            throw new NotImplementedException();
        }

        public double Length(Vectors<U> a)
        {
            throw new NotImplementedException();
        }

        public Vectors<U> MultiplicativeIdentity()
        {
            throw new NotImplementedException();
        }

        public Vectors<U> Multiply(Vectors<U> a, Vectors<U> b)
        {
            throw new NotImplementedException();
        }

        public Vectors<U> Negate(Vectors<U> a)
        {
            throw new NotImplementedException();
        }

        public U ScalarIdentity()
        {
            throw new NotImplementedException();
        }

        public Vectors<U> ScalarMultiply(Vectors<U> a, U s)
        {
            throw new NotImplementedException();
        }

        public Vectors<U> Substract(Vectors<U> a, Vectors<U> b)
        {
            throw new NotImplementedException();
        }
    }
    */
    
}
