using System;
using UnityEngine;

namespace PUtil.Plot
{
    ///<summary>
    /// The Plot Util class provides an easy way to create 2D and 3D plots like matlab
    ///</summary>
    public static class PlotUtil
    {
        public static void DrawArrow(Vector3 start, Vector3 end, Color color, float duration, float arrowLength = 0.2f)
        {
            float length = (end - start).magnitude;
            Debug.DrawLine(start, end, color, duration);

            Debug.DrawLine(end, start + (end - start) * (1f - arrowLength) + new Vector3(length * arrowLength, 0f, 0f), color, duration);
            Debug.DrawLine(end, start + (end - start) * (1f - arrowLength) + new Vector3(-length * arrowLength, 0f, 0f), color, duration);
            Debug.DrawLine(end, start + (end - start) * (1f - arrowLength) + new Vector3(0f, length * arrowLength, 0f), color, duration);
            Debug.DrawLine(end, start + (end - start) * (1f - arrowLength) + new Vector3(0f, -length * arrowLength, 0f), color, duration);
            Debug.DrawLine(end, start + (end - start) * (1f - arrowLength) + new Vector3(0f, 0f, length * arrowLength), color, duration);
            Debug.DrawLine(end, start + (end - start) * (1f - arrowLength) + new Vector3(0f, 0f, -length * arrowLength), color, duration);
        }
        public static void DrawArrow(Vector3 start, Vector3 end, Color color, float arrowLength = 0.2f)
        {
            DrawArrow(start, end, color, 0.01f, arrowLength);
        }
        public static void DrawArrow(Vector3 start, Vector3 end, float arrowLength = 0.2f)
        {
            DrawArrow(start, end, Color.white, 0.01f, arrowLength);
        }
        public static float GetMax(float[] x)
        {
            if (x.Length == 1) return x[0];
            float max = x[0];
            for(int i=1; i<x.Length; i++)
            {
                if (x[i] > max) max = x[i];
            }
            return max;
        }
        public static float GetMin(float[] x)
        {
            if (x.Length == 1) return x[0];
            float min = x[0];
            for (int i = 1; i < x.Length; i++)
            {
                if (x[i] < min) min = x[i];
            }
            return min;
        }
        public static void GetMinMax(float[] x, out float min, out float max)
        {
            if(x.Length == 1) { min = x[0]; max = x[0]; }
            min = x[0];
            max = x[0];
            for (int i = 1; i < x.Length; i++)
            {
                if (x[i] < min) min = x[i];
                if (x[i] > max) max = x[i];
            }

        }


        ///<summary>
        /// returns n number of points, evenly spaced between x1 and x2
        ///</summary>
        ///<param name="x1"> starting point</param>
        ///<param name="x2"> end point</param>
        ///<param name="n"> number of points</param>
        public static float[] Linspace(float x1, float x2, int n)
        {
            if (n == 1) return new float[]{ (x2 - x1) / 2f };
            if (n == 2) return new float[] { x1, x2 };

            float d = (x2 - x1) / (float)(n - 1); //delta
            float[] space = new float[n]; //the thing we return
            space[0] = x1;
            for(int i=1; i<n; i++)
            {
                space[i] = space[i - 1] + d; 
            }
            return space;
        }
        ///<summary>
        /// returns 100 points, evenly spaced between x1 and x2
        ///</summary>
        ///<param name="x1"> starting point</param>
        ///<param name="x2"> end point</param>
        public static float[] Linspace(float x1, float x2)
        {
            return Linspace(x1, x2, 100);
        }
        ///<summary>
        /// returns 100 points, evenly spaced between 0 and x2
        ///</summary>
        ///<param name="x2"> end point</param>
        public static float[] Linspace(float x2)
        {
            return Linspace(0f, x2, 100);
        }
        ///<summary>
        /// returns 100 points, evenly spaced between 0 and 1
        ///</summary>
        ///<param name="x2"> end point</param>
        public static float[] Linspace()
        {
            return Linspace(0f, 1f, 100);
        }


        ///<summary>
        /// returns n number of points, evenly spaced between x1 and x2
        ///</summary>
        ///<param name="x1"> starting point</param>
        ///<param name="x2"> end point</param>
        ///<param name="n"> number of points</param>
        public static double[] Linspace(double x1, double x2, int n)
        {
            if (n == 1) return new double[] { (x2 - x1) / 2f };
            if (n == 2) return new double[] { x1, x2 };

            double d = (x2 - x1) / (double)(n - 1); //delta
            double[] space = new double[n]; //the thing we return
            space[0] = x1;
            for (int i = 1; i < n; i++)
            {
                space[i] = space[i - 1] + d;
            }
            return space;
        }
        ///<summary>
        /// returns 100 points, evenly spaced between x1 and x2
        ///</summary>
        ///<param name="x1"> starting point</param>
        ///<param name="x2"> end point</param>
        public static double[] Linspace(double x1, double x2)
        {
            return Linspace(x1, x2, 100);
        }
        ///<summary>
        /// returns 100 points, evenly spaced between 0 and x2
        ///</summary>
        ///<param name="x2"> end point</param>
        public static double[] Linspace(double x2)
        {
            return Linspace(0.0, x2, 100);
        }



        ///<summary>
        /// returns n number of points, logarithmically spaced between a and b
        ///</summary>
        ///<param name="x"> base</param>
        ///<param name="a"> starting point x^a </param>
        ///<param name="b"> end point x^b </param>
        ///<param name="n"> number of points</param>
        public static float[] Logspace(float x, float a, float b, int n)
        {
            float[] ls = Linspace(a, b, n);
            float[] logspace = new float[n];
            for(int i=0; i<n; i++)
            {
                logspace[i] = (float)Math.Pow(x, ls[i]);
            }
            return logspace;
        }



        ///<summary>
        /// Draws the lines connecting the x and y coordinates.
        ///</summary>
        ///<param name="c"> Colour of the plot</param>
        ///<param name="duration"> How long should the plot persist</param>
        ///<param name="x"> x coordinates</param>
        ///<param name="y"> y coordinates</param>
        ///<param name="scale"> scale of the graph, by default it's unit sized</param>
        ///<param name="drawAxis"> Should we draw the axis</param>
        ///<param name="scale"> should the graph be scaled</param>
        public static void Plot2D(float[] x, float[] y, Color c, float duration, float scale = 1f, bool drawAxis = true, bool normalise = true)
        {
            if (x.Length != y.Length) throw new Exception("arrays don't match length");
            //we assume x[0] is min and x[n] is max
            float xmin = new float();
            float xmax = new float();
            float ymin = new float();
            float ymax = new float();

            if (drawAxis || normalise) //we only need to find min and max if we want to draw additional things
            {
                GetMinMax(x, out xmin, out xmax);
                GetMinMax(y, out ymin, out ymax);
            }

            float dx = 1f;
            float dy = 1f;
            if (normalise)
            {
                dx = xmax - xmin;
                dy = ymax - ymin;
            }

            for (int i = 1; i < x.Length; i++)
            {
                Vector3 v0 = new Vector3(x[i - 1]/dx, y[i - 1]/dy, 0f);
                Vector3 v1 = new Vector3(x[i]/dx, y[i]/dy, 0f);
                Debug.DrawLine(scale * v0, scale * v1, c, duration);
            }

            if(drawAxis)
            {
                DrawArrow(scale * new Vector3(xmin / dx, 0f, 0f), scale * new Vector3(xmax / dx, 0f, 0f), Color.red); //x Axis
                DrawArrow(scale * new Vector3(0f, ymin/dy, 0f), scale * new Vector3(0f, ymax/dy, 0f), Color.green); //Y axis
            }
        }

        public static void Plot2D(float[] x, float[] y, Color c, float scale = 1f, bool drawAxis = true, bool normalise = true)
        {
            Plot2D(x, y, c, 0.01f, scale, drawAxis, normalise);
        }

        public static void Plot2D(float[] x, float[] y,float scale = 1f, bool drawAxis = true, bool normalise = true)
        {
            Plot2D(x, y, Color.white, 0.01f, scale, drawAxis, normalise);
        }

        ///<summary>
        /// Draws the lines connecting the x and y coordinates with the height, It draws the height map as a grid
        ///</summary>
        ///<param name="f"> height function in terms of f(x,y) = z; has 2 float parameters and returns float</param>
        ///<param name="x"> x coordinates</param>
        ///<param name="y"> y coordinates</param>
        public static void Plot3D(float[] x, float[] y, Func<float, float, float> f, bool drawAxis = true, bool normalise = true)
        {
            for(int i=0; i<x.Length; i++)
            {
                for(int j=0; j<y.Length; j++)
                {
                    //here this is the i,j th coordinate
                    //float z0 = f(x[i], y[j]);
                    //float z1 = f(x[i + 1], y[j]);
                    //float z2 = f(x[i], y[j+1]);
                    float z0, z1, z2 = 0f;
                    z0 = f(x[i], y[j]); //current point
                    if (i < x.Length - 1)
                    {
                        z1 = f(x[i + 1], y[j]); //point on my right
                        Debug.DrawLine(new Vector3(x[i], z0, y[j]), new Vector3(x[i + 1], z1, y[j]));
                    }
                    if (j < y.Length - 1)
                    {
                        z2 = f(x[i], y[j + 1]); //point front of me
                        Debug.DrawLine(new Vector3(x[i], z0, y[j]), new Vector3(x[i], z2, y[j+1]));
                    }
                }
            }

            if(drawAxis)
            {
                DrawArrow(Vector3.zero, new Vector3(1f, 0f, 0f), Color.red);
                DrawArrow(Vector3.zero, new Vector3(0f, 1f, 0f), Color.blue);
                DrawArrow(Vector3.zero, new Vector3(0f, 0f, 1f), Color.green);
            }
        }


    }
}
