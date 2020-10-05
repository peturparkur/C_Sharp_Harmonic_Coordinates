using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Runtime.InteropServices;
//using System.IO;

public class VertexCalculation
{

    public ComputeShader cs; //this will be the matrixVectormult
    public int kernel;

    public string weightPath = "Assets/Animation/Weights.txt";

    ComputeBuffer matrixBuffer;
    ComputeBuffer vectorBuffer;
    ComputeBuffer vertOutBuffer;
    //ComputeBuffer output;

    public Vector3[] output;


    // Start is called before the first frame update
    void Start()
    {
        
    }

    public void SetBuffers(Matrix<float> matrix, Vector3[] verts, Vector3[] targetVerts)
    {
        matrixBuffer = new ComputeBuffer(matrix.rows * matrix.columns, sizeof(float)); //weights
        vectorBuffer = new ComputeBuffer(verts.Length, Marshal.SizeOf(typeof(Vector3))); //boundary mesh vertices
        vertOutBuffer = new ComputeBuffer(targetVerts.Length, Marshal.SizeOf(typeof(Vector3))); //output vertices
        matrixBuffer.SetData(matrix.values); //setting weights
        vectorBuffer.SetData(verts); //setting boundary verts
    }

    // Update is called once per frame
    public void RunShader(Matrix<float> matrix, Vector3[] verts, Vector3[] targetVerts)
    {
        Debug.Log("Run shader");
        //SetBuffers(matrix, verts, targetVerts);
        kernel = cs.FindKernel("Calculate");
        cs.SetBuffer(kernel, "matrixA", matrixBuffer);
        cs.SetBuffer(kernel, "vectorB", vectorBuffer);
        cs.SetBuffer(kernel, "Result", vertOutBuffer);
        cs.SetInt("columns", matrix.columns);
        cs.Dispatch(kernel, matrix.rows / 64, 1, 1);

        vertOutBuffer.GetData(targetVerts);
        Vector<Vector3> outverts = new Vector<Vector3>(targetVerts.Length, targetVerts);
        Debug.Log("targetVerts = " + outverts.ToString());
    }
}
