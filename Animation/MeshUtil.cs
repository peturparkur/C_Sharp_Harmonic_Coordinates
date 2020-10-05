using System;
using UnityEngine;
using System.Collections.Generic;
using PUtil.Operators;

namespace PUtil.MeshUtil
{
    public static class MeshUtil
    {
        public static Vector3 PointOnTriangle(Vector3 v1, Vector3 v2, Vector3 v3, float r1, float r2)
        {
            float sqrt1 = Mathf.Sqrt(r1);
            return (1 - sqrt1) * v1 + (sqrt1 * (1 - r2)) * v2 + (r2 * sqrt1) * v3; //this should be a random point on a triangle
        }

        public static Vector<float> PointOnTriangle(Vector<float> v1, Vector<float> v2, Vector<float> v3, float r1, float r2)
        {
            float sqrt1 = Mathf.Sqrt(r1);
            return (1 - sqrt1) * v1 + (sqrt1 * (1 - r2)) * v2 + (r2 * sqrt1) * v3; //this should be a random point on a triangle
        }

        public static Mesh GenerateTetrahedron() //we generate this around the origin
        {
            Mesh mesh = new Mesh();
            mesh.vertices = new Vector3[] { Vector3.one, new Vector3(1, -1, -1), new Vector3(-1, -1, 1), new Vector3(-1, 1, -1) };
            mesh.triangles = new int[] { 0, 1, 3 };
            mesh.RecalculateBounds();
            mesh.RecalculateNormals();
            mesh.RecalculateTangents();
            return mesh;
        }

        public static Mesh GenerateTetrahedron2()
        {
            Mesh mesh = new Mesh();
            float sqrt2 = (float)Math.Sqrt(2);
            mesh.vertices = new Vector3[] { new Vector3(1, -1 / sqrt2, 0), new Vector3(-1, -1f / sqrt2, 0), new Vector3(0, 1 / sqrt2, 1), new Vector3(0, 1 / sqrt2, -1) };
            mesh.triangles = new int[] { 0, 1, 3 };
            mesh.RecalculateBounds();
            mesh.RecalculateNormals();
            mesh.RecalculateTangents();
            return mesh;
        }

        public static Mesh GenerateTetrahedron3(float scale = 1f) //by default this generates a tetrahedron which fits into a unit sphere
        {
            Mesh mesh = new Mesh();
            //float sqrt2 = (float)Math.Sqrt(2);
            mesh.vertices = new Vector3[] 
            {
                new Vector3(Mathf.Sqrt(8f/9f), -1f / 3f, 0f) * scale,
                new Vector3(-Mathf.Sqrt(2f/9f), -1f / 3f, Mathf.Sqrt(2f / 3f)) * scale,
                new Vector3(-Mathf.Sqrt(2f/9f), -1f / 3f, -Mathf.Sqrt(2f / 3f)) * scale,
                new Vector3(0,1,0) * scale
            };

            mesh.triangles = new int[] {
                0, 2, 3,
                0, 3, 1,
                1, 3, 2,
                0, 1, 2};

            mesh.RecalculateBounds();
            mesh.RecalculateNormals();
            mesh.RecalculateTangents();
            return mesh;
        }

        public static Vector3 BarycentricCoordinates(Vector3 v1, Vector3 v2, Vector3 v3, Vector3 p) //described in clockwise order
        {
            //https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/barycentric-coordinates#:~:text=Barycentric%20coordinates%20can%20be%20used%20to%20express%20the%20position%20of,the%20three%20triangle%27s%20vertices%20themselves.
            //p is the point we're trying to describe
            //we assume p is flattened on the triangle
            //u,v,w => (v1,v2,v3)

            //u = Area(1,2,p) / Area(1,2,3)

            float totalArea = Vector3.Cross(v2 - v1, v3 - v1).magnitude / 2f; //this is the total area of the triangle
            float u = Vector3.Cross(v2 - v1, p - v1).magnitude / 2f;
            u = u / totalArea;

            float v = Vector3.Cross(v3 - v2, p - v2).magnitude / 2f;
            v = v / totalArea;

            float w = 1f - u - v;

            return new Vector3(u, v, w); //these are associated in the order of the triangle vectors
            //p = u*v1 + v*v2 + w*v3
        }

        public static bool RaycastTriangle(Vector3 origin, Vector3 rayVector, float length, Vector3 v0, Vector3 v1, Vector3 v2, out Vector3 intersection)
        {
            //we assume that the rayvector is normalised
            //v1,v2,v3 are the triangle vertices
            intersection = new Vector3();
            Vector3 edge1 = v1 - v0;
            Vector3 edge2 = v2 - v0;
            Vector3 h = Vector3.Cross(rayVector, edge2);
            float a = Vector3.Dot(edge1, h);
            if (Mathf.Abs(a) < Mathf.Epsilon) return false; //ray is parallel to triangle
            float f = 1f / a;

            Vector3 s = origin - v0;
            float u = f * Vector3.Dot(s, h);
            if (u < 0 || u > 1) return false;

            Vector3 q = Vector3.Cross(s, edge1);
            float v = f * Vector3.Dot(rayVector, q);
            if (v < 0 || u + v > 1) return false;

            float t = f * Vector3.Dot(edge2, q);
            if(t>Mathf.Epsilon && t<= length)
            {
                intersection = origin + rayVector * t;
                return true;
            }
            else
            {
                return false;
            }
        }

        public static Vector3[] RaycastMesh(Vector3 origin, Vector3 rayVector, float length, Mesh mesh)
        {
            List<Vector3> array = new List<Vector3>();
            for(int i=0; i<mesh.triangles.Length; i+=3)
            {
                Vector3 v = new Vector3();
                bool h = RaycastTriangle(origin, rayVector, length, mesh.vertices[mesh.triangles[i]], mesh.vertices[mesh.triangles[i + 1]], mesh.vertices[mesh.triangles[i + 2]], out v);
                if(h)
                {
                    if(array.Contains(v))
                    {
                        continue;
                    }
                    //now we also need to check if it's too close to an already existing point
                    bool c = false;
                    for(int j=0; j<array.Count; j++)
                    {
                        float d = (array[j] - v).sqrMagnitude;
                        if(d<0.0001f)
                        {
                            //it's too near
                            c=true;
                            break;
                        }
                    }
                    if (c) continue;

                    array.Add(v);
                }
            }
            //we should order the points based on the distance
            return array.ToArray();
        }

        public static void OrderVectorArray(Vector3[] array, Vector3 origin) //sorts it in ascending order
        {
            //for now super quick bubble sort
            float d0 = 0f;
            float d1 = 0f;
            bool ordered = false;
            int n = 0;
            while(!ordered)
            {
                ordered = true;
                for(int i=0; i<array.Length-1; i++)
                {
                    d0 = (array[i] - origin).sqrMagnitude;
                    d1 = (array[i+1] - origin).sqrMagnitude;

                    if(d0 > d1) //swap them
                    {
                        ordered = false;
                        Vector3 tmp = array[i];
                        array[i] = array[i + 1];
                        array[i + 1] = tmp;
                    }
                }
                if (n > array.Length) return;
                n++;
            }
        }
    }
}
