using System;
using UnityEngine;
using PUtil.MeshUtil;

public class Harmonic
{
    public int numCells;
    public Tensor<Cell> grid;//= new Tensor<Cell>(new int[]{100,100,100});
    public Tensor<float> gridValues; //voxel values for Laplace

    //we require this so we can easily calculate the positions
    public Matrix<float> cellSize = new Matrix<float>(3, 2); //[0,0] means minX, [0,1] maxX; [1,0] means minZ and so on
    public Vector<float> d = new Vector<float>(3,true); //this might not be necessary

    public void Create3DGrid(int _numCells)
    {
        numCells = _numCells;
        grid = new Tensor<Cell>(new int[] { numCells, numCells, numCells }); //creates a numCells x numCells x numCells space
        gridValues = new Tensor<float>(new int[] { numCells, numCells, numCells });
        for (int i = 0; i < grid.values.Length; i++)
        {
            grid.values[i] = new Cell(CType.Undefined);//, 0); //untyped with 0 value
            gridValues.values[i] = 0f;
        }
    }

    public void SetGridValues(float x)
    {
        for(int i=0; i<gridValues.values.Length; i++)
        {
            gridValues.values[i] = x;
        }
    }

    public void Create2DGrid(int numCells)
    {
        grid = new Tensor<Cell>(new int[] { numCells, numCells });
        for (int i = 0; i < grid.values.Length; i++)
        {
            grid.values[i] = new Cell(CType.Undefined);//, 0); //untyped with 0 value
            gridValues.values[i] = 0f;
        }
    }

    public void Init()
    {
        d = new Vector<float>(3, false);
        for (int i = 0; i < d.rows; i++)
        {
            d[i] = (cellSize[i, 1] - cellSize[i,0]) / (numCells); //calculating dx this is the smallest difference between voxels
            //we need the samples to be atleast this close to each other
        }
    }

    public int[] VectorToGrid(Vector<float> v) //returns the voxel coordinates
    {
        float x = v[0];
        float y = v[0];
        float z = v[0];

        int ix, iy, iz = 0;
        if (x > cellSize[0, 1]) throw new Exception("vector.x is outside the voxels");//ix = numCells - 1; //last x value in voxel
        if (y > cellSize[1, 1]) throw new Exception("vector.y is outside the voxels");//iy = numCells - 1; //last y value in voxels
        if (z > cellSize[2, 1]) throw new Exception("vector.z is outside the voxels");//iz = numCells - 1; //last z value in voxels

        //we assume that v is inside the voxel range between the min and max
        //ix = (int)((x - cellSize[0, 0]) * (numCells - 1) / (cellSize[0, 1] - cellSize[0, 0]));
        ix = (int) ((x - cellSize[0, 0]) / d[0]);
        iy = (int) ((y - cellSize[1, 0]) / d[1]);
        iz = (int) ((z - cellSize[2, 0]) / d[2]);

        return new int[] { ix, iy, iz }; //x,y,z voxel coordinates
    }

    //this is unity specific
    public int[] VectorToGrid(Vector3 v) //returns the voxel coordinates
    {
        float x = v.x;
        float y = v.y;
        float z = v.z;

        int ix, iy, iz = 0;
        if (x > cellSize[0, 1]) throw new Exception("vector.x is outside the voxels");//ix = numCells - 1; //last x value in voxel
        if (y > cellSize[1, 1]) throw new Exception("vector.y is outside the voxels");//iy = numCells - 1; //last y value in voxels
        if (z > cellSize[2, 1]) throw new Exception("vector.z is outside the voxels");//iz = numCells - 1; //last z value in voxels

        //we assume that v is inside the voxel range between the min and max
        //ix = (int)((x - cellSize[0, 0]) * (numCells - 1) / (cellSize[0, 1] - cellSize[0, 0]));
        ix = (int)((x - cellSize[0, 0]) / d[0]);
        //Debug.Log("ix = " + ((x - cellSize[0, 0]) / d[0]));
        //Debug.Log("ix as int = " + ix);
        iy = (int)((y - cellSize[1, 0]) / d[1]);
        iz = (int)((z - cellSize[2, 0]) / d[2]);

        return new int[] { ix, iy, iz }; //x,y,z voxel coordinates
    }


    public Vector<float> GridToVector(int[] c) //returns the starting point of the voxel
    {
        //c is th x,y,z coordinates
        Vector<float> v = new Vector<float>(3);
        v[0] = cellSize[0, 0] + c[0] * d[0]; //this would mark the left hand corner of the grid
        v[1] = cellSize[1, 0] + c[1] * d[1];
        v[2] = cellSize[2, 0] + c[2] * d[2];

        return v;
    }

    public Vector<float> GridToVectorCenter(int[] c) //returns the center of the voxel
    {
        Vector<float> v = new Vector<float>(3);
        v[0] = cellSize[0, 0] + (c[0] - 0.5f) * d[0]; //this would mark the center of the voxel
        v[1] = cellSize[1, 0] + (c[1] - 0.5f) * d[1];
        v[2] = cellSize[2, 0] + (c[2] - 0.5f) * d[2];
        return v;
    }
    public Vector3 GridToVector3Center(int[] c) //returns the center of the voxel
    {
        Vector3 v = new Vector3();
        v[0] = cellSize[0, 0] + (c[0] + 0.5f) * d[0]; //this would mark the center of the voxel
        v[1] = cellSize[1, 0] + (c[1] + 0.5f) * d[1];
        v[2] = cellSize[2, 0] + (c[2] + 0.5f) * d[2];
        return v;
    }
}

public enum CType //this is to store the information if the voxel is interior/boundary and such
{
    Undefined = 0,
    Boundary,
    Exterior,
    Interior
}

/*public struct Cell //the voxel cells
{
    public CType type; //-1 untyped, 0 Boundary, 1 Exterior, 2 Interior
    //public double value; //this represents the weight

    public Cell(CType _type)//, double _val)
    {
        type = _type;
        //value = _val;
    }
}*/
