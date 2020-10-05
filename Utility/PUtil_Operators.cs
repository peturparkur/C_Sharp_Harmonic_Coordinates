using System;
//using System.Collections.Generic;
using System.Linq.Expressions;
//using System.Text;
//using System.Threading.Tasks;


namespace PUtil.Operators
{
    ///<summary>
    /// The Expression Util class provides a faster way to creating generic operators
    ///</summary>
    public static class ExpressionUtil
    {
        public static Func<Targ1, Tresult> CreateExpression<Targ1, Tresult>(Func<Expression, UnaryExpression> body)
        //Create a function with 1 argument that follows the lambda expression
        {
            ParameterExpression pm1 = Expression.Parameter(typeof(Targ1), "pm1");
            try
            {
                return Expression.Lambda<Func<Targ1, Tresult>>(body(pm1), pm1).Compile();
            }
            catch (Exception ex)
            {
                string msg = ex.Message;
                return delegate { throw new InvalidOperationException(msg); };
            }
        }

        public static Func<Targ1, Targ2, Tresult> CreateExpression<Targ1, Targ2, Tresult>(Func<Expression, Expression, BinaryExpression> body)
        //Create a function with 2 argument that follows the lambda expression
        {
            ParameterExpression pm1 = Expression.Parameter(typeof(Targ1), "pm1");
            ParameterExpression pm2 = Expression.Parameter(typeof(Targ2), "pm2");
            try
            {
                try
                {
                    return Expression.Lambda<Func<Targ1, Targ2, Tresult>>(body(pm1, pm2), pm1, pm2).Compile();
                }
                catch(InvalidOperationException)
                {
                    if (!(typeof(Targ1) == typeof(Tresult)) && (typeof(Targ2) == typeof(Tresult)))
                    {
                        Expression castL = typeof(Targ1) == typeof(Tresult) ? (Expression)pm1 : (Expression)Expression.Convert(pm1, typeof(Tresult)); //cast left argument to result
                        Expression castR = typeof(Targ2) == typeof(Tresult) ? (Expression)pm2 : (Expression)Expression.Convert(pm2, typeof(Tresult));

                        return Expression.Lambda<Func<Targ1, Targ2, Tresult>>(body(castL, castR), pm1, pm2).Compile();
                    }
                    else throw;
                }
            }
            catch (Exception ex)
            {
                string msg = ex.Message;
                return delegate { throw new InvalidOperationException(msg); };
            }
        }
    }

    ///<summary>
    /// The Operator class provides access to standard operators (eg.: Addition, Multiplication) access to generic types
    /// Targ1 is the Type for to the argument, Tresult is the Type of the return
    ///</summary>
    public static class Operator<Targ1, Tresult>
    {
        public static readonly Func<Targ1, Tresult> convert;
        public static readonly Func<Tresult, Targ1, Tresult> add, substract, multiply, divide, power; //eg.: float * int = float or double * float = double

        static Operator()
        {
            convert = ExpressionUtil.CreateExpression<Targ1, Tresult>(body => Expression.Convert(body, typeof(Tresult)));
            add = ExpressionUtil.CreateExpression<Tresult, Targ1, Tresult>(Expression.Add);
            substract = ExpressionUtil.CreateExpression<Tresult, Targ1, Tresult>(Expression.Subtract);
            multiply = ExpressionUtil.CreateExpression<Tresult, Targ1, Tresult>(Expression.Multiply);
            divide = ExpressionUtil.CreateExpression<Tresult, Targ1, Tresult>(Expression.Divide);
            power = ExpressionUtil.CreateExpression<Tresult, Targ1, Tresult>(Expression.Power);
        }
    }

    ///<summary>
    /// The Operator class provides access to standard operators (eg.: Addition, Multiplication) access to generic types
    /// Type gives the type of the argument and result
    ///</summary>
    public static class Operator<Type>
    {
        public static readonly Func<Type, Type> negate, not;
        public static readonly Func<Type, Type, bool> equal, notEqual, lessThan, greaterThan, greaterThanOrEqual, lessThanOrEqual;
        public static readonly Func<Type, Type, Type> add, substract, multiply, divide;

        static Operator()
        {
            negate = ExpressionUtil.CreateExpression<Type, Type>(Expression.Negate);
            not = ExpressionUtil.CreateExpression<Type, Type>(Expression.Not);

            equal = ExpressionUtil.CreateExpression<Type, Type, bool>(Expression.Equal);
            notEqual = ExpressionUtil.CreateExpression<Type, Type, bool>(Expression.NotEqual);
            lessThan = ExpressionUtil.CreateExpression<Type, Type, bool>(Expression.LessThan);
            greaterThan = ExpressionUtil.CreateExpression<Type, Type, bool>(Expression.GreaterThan);
            lessThanOrEqual = ExpressionUtil.CreateExpression<Type, Type, bool>(Expression.LessThanOrEqual);
            greaterThanOrEqual = ExpressionUtil.CreateExpression<Type, Type, bool>(Expression.GreaterThanOrEqual);

            add = ExpressionUtil.CreateExpression<Type, Type, Type>(Expression.Add);
            substract = ExpressionUtil.CreateExpression<Type, Type, Type>(Expression.Subtract);
            multiply = ExpressionUtil.CreateExpression<Type, Type, Type>(Expression.Multiply);
            divide = ExpressionUtil.CreateExpression<Type, Type, Type>(Expression.Divide);
        }
    }
}
