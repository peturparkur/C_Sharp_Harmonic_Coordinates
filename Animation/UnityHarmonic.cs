using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using PUtil.MeshUtil;

public class UnityHarmonic : MonoBehaviour
{
    public Mesh mesh;
    public Mesh targetMesh;
    Harmonic harmonic;
    public uint gridSize = 2 * 2 * 2 * 2;

    public Matrix<float> weightMatrix;

    public List<int> Interiors; //array for all interior cells
    public List<int> Boundaries; //array of all boundary cells
    //Might need the array of interior points for optimised looping for smoothing

    private void Start()
    {
        Init();
        Interiors = new List<int>();
        harmonic = new Harmonic();
        harmonic.cellSize = new Matrix<float>(3, 2);

        harmonic.cellSize[0, 0] = -1.2f; //min x
        harmonic.cellSize[0, 1] = 1.2f; //max x

        harmonic.cellSize[1, 0] = -1.2f; //min y
        harmonic.cellSize[1, 1] = 1.2f; //max y

        harmonic.cellSize[2, 0] = -1.2f; //min z
        harmonic.cellSize[2, 1] = 1.2f; //max z
        harmonic.Create3DGrid((int)gridSize);
        harmonic.Init();

        /*
        Matrix<float> bounding = TriangleBoundingBox(mesh.vertices[0], mesh.vertices[1], mesh.vertices[2]);
        Debug.Log("bounding = " + bounding.ToString());
        Vector<int> count = new Vector<int>(3);
        int max = 0;
        for(int i=0; i<bounding.rows; i++)
        {
            count[i] = (int)((bounding[i, 1] - bounding[i, 0]) / harmonic.d[i])*2;
            if (count[i] > max) max = count[i];
        }
        */

        //we not create the algorithm for all meshes
        AddRegularTrianglePointsForMesh(mesh); //this adds the boundaries for each triangle of the mesh
        AddExternalPoints(); //setting the voxels for external points - This relies on the boundaries being already set
        DebugVoxelBounds(false, false, false); //debugging

        //AddTrianglePoints(25, mesh.vertices[0], mesh.vertices[1], mesh.vertices[2]);
        //AddRegularTrianglePoints(max, mesh.vertices[0], mesh.vertices[1], mesh.vertices[2]);


    }

    private void Update()
    {
        //Testing();
        MeshTesting();

        targetMesh = MeshUtil.GenerateTetrahedron3(0.5f); //centered at origin, fits in half a radius sphere
        weightMatrix = new Matrix<float>(targetMesh.vertices.Length, mesh.vertices.Length); //if we want to change the targetmeshvertex i get weights as weights[i,1] => 2nd weight of the ith targetmesh vertex
        StartDiffusion(mesh.vertices[0], 100);
        for (int i = 0; i < weightMatrix.rows; i++)
        {
            weightMatrix[i, 0] = harmonic.gridValues[harmonic.VectorToGrid(targetMesh.vertices[i])];
        }
        Debug.Log("WeightMatrix = ");
        Debug.Log(weightMatrix.ToString());
    }

    void StartDiffusion(Vector3 vertex, int numSteps=1)
    {
        harmonic.SetGridValues(0f); //to reset the grid values
        Vector<int> pos = new Vector<int>(3, harmonic.VectorToGrid(vertex));
        //harmonic.grid[pos.ToArray()] = new Cell(CType.Interior); //so the other interiors copy from it
        harmonic.gridValues[pos.ToArray()] = 1f; //the flooding value

        Debug.Log("vertex position = ");
        Debug.Log(pos.ToString());

        Tensor<float> tempGrid = new Tensor<float>(harmonic.gridValues); //copy values from voxel to temp storage for updating


        //find neightbours to the boundary pos
        bool foundN = false;
        pos[0] += 1;
        CType neig = harmonic.grid[new int[] { pos[0], pos[1], pos[2] }].type;
        if (neig == CType.Interior) { Debug.Log("Neighbour found at "); Debug.Log(pos.ToString()); foundN = true; }

        pos[0] -= 2;
        neig = harmonic.grid[new int[] { pos[0], pos[1], pos[2] }].type;
        if (neig == CType.Interior) { Debug.Log("Neighbour found at "); Debug.Log(pos.ToString()); foundN = true; }

        pos[0] += 1;
        pos[1] -= 1;
        neig = harmonic.grid[new int[] { pos[0], pos[1], pos[2] }].type;
        if (neig == CType.Interior) { Debug.Log("Neighbour found at "); Debug.Log(pos.ToString()); foundN = true; }

        pos[1] += 2;
        neig = harmonic.grid[new int[] { pos[0], pos[1], pos[2] }].type;
        if (neig == CType.Interior) { Debug.Log("Neighbour found at "); Debug.Log(pos.ToString()); foundN = true; }

        pos[1] -= 1;
        pos[2] -= 1;
        neig = harmonic.grid[new int[] { pos[0], pos[1], pos[2] }].type;
        if (neig == CType.Interior) { Debug.Log("Neighbour found at "); Debug.Log(pos.ToString()); foundN = true; }

        pos[1] += 2;
        neig = harmonic.grid[new int[] { pos[0], pos[1], pos[2] }].type;
        if (neig == CType.Interior) { Debug.Log("Neighbour found at "); Debug.Log(pos.ToString()); foundN = true; }

        Debug.Log("foundN = " + foundN.ToString());
        if(foundN == false) { Debug.Log("We didn't find a neighbour!!!"); }

        //int[] pos = harmonic.VectorToGrid(vertex); //grid position of the vertex

        bool noValid = false; //to make sure it's not surrounded by boundaries

        //we start looping across the interior points of the grid to smooth out values
        //harmonic.grid[pos.ToArray()] = new Cell(CType.Interior);
        for (int t = 0; t < numSteps; t++)
        {
            for (int i = 0; i < Interiors.Count; i++) //looping through interior points
            {
                //we want to check the neightbours in 3D
                //for simplicity let's just check the voxel neightbours
                int[] accessor = harmonic.grid.IndexToAccessor(Interiors[i]);
                //this is a 3d coordinate

                //we should only take into account the case when it's not a boundary
                CType type = CType.Undefined;
                int count = 0;
                float f0 = 0f;//harmonic.gridValues[new int[] { accessor[0], accessor[1] - 1, accessor[2] }]; //down
                if (harmonic.grid[new int[] { accessor[0], accessor[1] - 1, accessor[2] }].type == CType.Interior)
                {
                    count++;
                    f0 = harmonic.gridValues[new int[] { accessor[0], accessor[1] - 1, accessor[2] }];
                }

                float f1 = 0f;//harmonic.gridValues[new int[] { accessor[0], accessor[1] - 1, accessor[2] }]; //down
                if (harmonic.grid[new int[] { accessor[0] - 1, accessor[1], accessor[2] }].type == CType.Interior)
                {
                    count++;
                    f1 = harmonic.gridValues[new int[] { accessor[0] - 1, accessor[1], accessor[2] }]; //left
                }

                float f2 = 0f;//harmonic.gridValues[new int[] { accessor[0], accessor[1] - 1, accessor[2] }]; //down
                if (harmonic.grid[new int[] { accessor[0] + 1, accessor[1], accessor[2] }].type == CType.Interior)
                {
                    count++;
                    f2 = harmonic.gridValues[new int[] { accessor[0] + 1, accessor[1], accessor[2] }]; //left
                }

                float f3 = 0f;//harmonic.gridValues[new int[] { accessor[0], accessor[1] - 1, accessor[2] }]; //down
                if (harmonic.grid[new int[] { accessor[0], accessor[1], accessor[2] - 1 }].type == CType.Interior)
                {
                    count++;
                    f3 = harmonic.gridValues[new int[] { accessor[0], accessor[1], accessor[2] - 1 }]; //left
                }

                float f4 = 0f;//harmonic.gridValues[new int[] { accessor[0], accessor[1] - 1, accessor[2] }]; //down
                if (harmonic.grid[new int[] { accessor[0], accessor[1], accessor[2] + 1 }].type == CType.Interior)
                {
                    count++;
                    f4 = harmonic.gridValues[new int[] { accessor[0], accessor[1], accessor[2] + 1 }]; //left
                }

                float f5 = 0f;//harmonic.gridValues[new int[] { accessor[0], accessor[1] - 1, accessor[2] }]; //down
                if (harmonic.grid[new int[] { accessor[0], accessor[1] + 1, accessor[2] }].type == CType.Interior)
                {
                    count++;
                    f5 = harmonic.gridValues[new int[] { accessor[0], accessor[1] + 1, accessor[2] }]; //left
                }
                if (count == 0) throw new System.Exception("Surrounded by boundaries and exteriors!!! Shouldn't happen");

                //float f1 = harmonic.gridValues[new int[] { accessor[0]-1, accessor[1], accessor[2] }]; //left
                //float f2 = harmonic.gridValues[new int[] { accessor[0]+1, accessor[1], accessor[2] }]; //right
                //float f3 = harmonic.gridValues[new int[] { accessor[0], accessor[1], accessor[2]-1 }]; //back
                //float f4 = harmonic.gridValues[new int[] { accessor[0], accessor[1], accessor[2]+1 }]; //forward
                //float f5 = harmonic.gridValues[new int[] { accessor[0], accessor[1] + 1, accessor[2] }]; //up
                float average = (f0 + f1 + f2 + f3 + f4 + f5) / (float)count;
                if (Mathf.Abs(average) > Mathf.Epsilon) Debug.Log("something is more than 0");
                tempGrid[accessor] = average; //six point average
            }
            harmonic.gridValues = new Tensor<float>(tempGrid); //copy updates back to original
        }
    }

    void AddExternalPoints()
    {
        for(int i=0; i<harmonic.grid.dim[0]; i++)
        {
            for (int j = 0; j < harmonic.grid.dim[1]; j++)
            {
                bool foundBoundary = false;
                for (int k = 0; k < harmonic.grid.dim[2]; k++)
                {
                    //we starting setting every gridpoint to external until we find a boundary point
                    if (harmonic.grid[new int[] { i, j, k }].type == CType.Boundary)
                    {
                        foundBoundary = true;
                        break;
                    }
                    harmonic.grid[new int[] { i, j, k }] = new Cell(CType.Exterior);
                }

                if (foundBoundary)//if we found a boundary we now loop backwards to fill the other side
                {
                    foundBoundary = false;
                    for (int k = harmonic.grid.dim[2] - 1; k >= 0; k--)
                    {
                        if (harmonic.grid[new int[] { i, j, k }].type == CType.Boundary)
                        {
                            //break;
                            foundBoundary = true;
                        }
                        if (!foundBoundary)
                        {
                            harmonic.grid[new int[] { i, j, k }] = new Cell(CType.Exterior);
                        }
                        else
                        {
                            if(harmonic.grid[new int[] { i, j, k }].type == CType.Undefined)
                            {
                                harmonic.grid[new int[] { i, j, k }] = new Cell(CType.Interior);
                                Interiors.Add(harmonic.grid.Index(new int[] { i, j, k }));
                            }
                            if(harmonic.grid[new int[] { i, j, k }].type == CType.Exterior)
                            {
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    void DebugVoxelBounds(bool DrawEmpty = false, bool drawExternal = false, bool drawInternal = true)
    {
        for (int i = 0; i < harmonic.grid.dim[0]; i++)
        {
            for (int j = 0; j < harmonic.grid.dim[1]; j++)
            {
                for (int k = 0; k < harmonic.grid.dim[2]; k++)
                {
                    CType type = harmonic.grid[new int[] { i, j, k }].type;

                    if (type == CType.Boundary)
                    {
                        Debug.Log("(i,j,k)=(" + i + "," + j + "," + k + ") => " + harmonic.grid[new int[] { i, j, k }].type.ToString());
                        DebugExtension.DebugWireSphere(harmonic.GridToVector3Center(new int[] { i, j, k }), Color.black, harmonic.d[0] / 2f);
                    }
                    else
                    {
                        if (DrawEmpty)
                        {
                            DebugExtension.DebugWireSphere(harmonic.GridToVector3Center(new int[] { i, j, k }), Color.white, harmonic.d[0] / 2f);
                        }
                        if(drawExternal && type==CType.Exterior)
                        {
                            DebugExtension.DebugWireSphere(harmonic.GridToVector3Center(new int[] { i, j, k }), Color.red, harmonic.d[0] / 2f);
                        }
                        if(drawInternal && type==CType.Interior)
                        {
                            DebugExtension.DebugWireSphere(harmonic.GridToVector3Center(new int[] { i, j, k }), Color.blue, harmonic.d[0] / 2f);
                        }
                    }
                }
            }
        }
    }

    void AddRegularTrianglePointsForMesh(Mesh mesh)
    {
        for (int t = 0; t < mesh.triangles.Length; t+=3)
        {
            //if (t != 2*3) continue;

            Vector3 v1 = mesh.vertices[mesh.triangles[t]];
            Vector3 v2 = mesh.vertices[mesh.triangles[t + 1]];
            Vector3 v3 = mesh.vertices[mesh.triangles[t + 2]];
            Matrix<float> bounds = TriangleBoundingBox(v1, v2, v3);
            Vector<int> num = new Vector<int>(3);
            int max = 0;
            for (int i = 0; i < num.rows; i++)
            {
                num[i] = (int)((bounds[i, 1] - bounds[i, 0]) / harmonic.d[i]);
                Debug.Log("For triangle " + t + ", the num["+i+"] = " + num[i]);
                if (num[i] > max) max = num[i];
            }
            Debug.Log("For triangle " + t + ", the max = " + max);
            AddRegularTrianglePoints(max*4, v1, v2, v3);
        }
    }

    Matrix<float> TriangleBoundingBox(Vector3 v1, Vector3 v2, Vector3 v3)
    {
        Matrix<float> m = new Matrix<float>(3, 2); //[0,0] minX; [0,1] => maxX
        m[0, 0] = v1[0];
        m[0, 1] = v1[0];
        m[1, 0] = v1[1];
        m[1, 1] = v1[1];
        m[2, 0] = v1[2];
        m[2, 1] = v1[2];

        for (int i = 0; i < m.rows; i++)
        {
            if (v2[i] < m[i, 0]) m[i, 0] = v2[i];
            if (v3[i] < m[i, 0]) m[i, 0] = v3[i];
            if (v2[i] > m[i, 1]) m[i, 1] = v2[i];
            if (v3[i] > m[i, 1]) m[i, 1] = v3[i];
        }
        return m;
    }

    void AddTrianglePoints(int count, Vector3 v1, Vector3 v2, Vector3 v3)
    {
        for(int i=0; i<count; i++)
        {
            float r1 = Random.value;
            float r2 = Random.value;
            Vector3 v = MeshUtil.PointOnTriangle(v1, v2, v3, r1, r2);
            int[] p = harmonic.VectorToGrid(v);
            harmonic.grid[p] = new Cell(CType.Boundary);
        }
    }

    void AddRegularTrianglePoints(Vector<int> count, Vector3 v1, Vector3 v2, Vector3 v3)
    {
        float dx = 1f / (count[0]-1);
        float dy = 1f / (count[1]-1);
        //Debug.Log("number of samples = " + (count[0] * count[1]));

        for(int i=0; i<count[0]; i++)
        {
            for(int j=0; j<count[1]; j++)
            {
                Vector3 v = MeshUtil.PointOnTriangle(v1, v2, v3, i * dx, j * dy);
                int[] p = harmonic.VectorToGrid(v);
                harmonic.grid[p] = new Cell(CType.Boundary);

                int p1 = harmonic.grid.Index(p);
                if(!Boundaries.Contains(p1)) //if the boundaries don't already contain this point add this into the list
                {
                    Boundaries.Add(p1);
                }
            }
        }
    }

    void AddRegularTrianglePoints(int count, Vector3 v1, Vector3 v2, Vector3 v3, bool debug=false)
    {
        float dx = 1f / (count - 1);
        //Debug.Log("number of samples = " + (count * count));

        for (int i = 0; i < count; i++)
        {
            for (int j = 0; j < count; j++)
            {
                Vector3 v = MeshUtil.PointOnTriangle(v1, v2, v3, i * dx, j * dx);

                if (debug) DebugExtension.DebugPoint(v, Color.red, 0.2f);

                int[] p = harmonic.VectorToGrid(v);
                harmonic.grid[p] = new Cell(CType.Boundary);
            }
        }
    }

    public void Init()
    {
        mesh = new Mesh();
        mesh.vertices = new Vector3[] { Vector3.zero, Vector3.up, Vector3.right};
        mesh.triangles = new int[] { 0, 1, 2 }; //note unity needs clockwise winding order for visible polygons
        mesh.RecalculateNormals();
        mesh.RecalculateBounds();
        mesh.RecalculateTangents();

        mesh = MeshUtil.GenerateTetrahedron3();

        GetComponent<MeshFilter>().mesh = mesh;
    }

    void Testing2()
    {
        Debug.Log("number of cells = " + harmonic.numCells);
        Debug.Log("difference = " + harmonic.d.ToString());
        Debug.Log("cellsize = " + harmonic.cellSize.ToString());

        Vector3 v1 = mesh.vertices[0];
        Vector3 v2 = mesh.vertices[1];
        Vector3 v3 = mesh.vertices[2];

        Debug.Log("STARTING");
        int[] vp = harmonic.VectorToGrid(v1);
        Debug.Log("Vertex1 position " + v1);
        Debug.Log("vertex1 grid position = (" + vp[0] + "," + vp[1] + "," + vp[2] + ")");
        vp = harmonic.VectorToGrid(v2);
        Debug.Log("Vertex2 position " + v2);
        Debug.Log("vertex2 grid position = (" + vp[0] + "," + vp[1] + "," + vp[2] + ")");
        vp = harmonic.VectorToGrid(v3);
        Debug.Log("Vertex3 position " + v3);
        Debug.Log("vertex3 grid position = (" + vp[0] + "," + vp[1] + "," + vp[2] + ")");
    }

    void MeshTesting()
    {
        for(int i=0; i<mesh.vertices.Length; i++)
        {
            DebugExtension.DebugPoint(mesh.vertices[i], Color.white, 0.5f);
        }
    }

    void Testing()
    {
        Vector3 v1 = mesh.vertices[0];
        Vector3 v2 = mesh.vertices[1];
        Vector3 v3 = mesh.vertices[2];

        Debug.DrawRay(mesh.vertices[0], -Vector3.forward);
        Debug.DrawRay(mesh.vertices[1], -Vector3.forward);
        Debug.DrawRay(mesh.vertices[2], -Vector3.forward);

        Debug.DrawRay(MeshUtil.PointOnTriangle(v1, v2, v3, 0f, 0f), Vector3.forward, Color.red);
        Debug.DrawRay(MeshUtil.PointOnTriangle(v1, v2, v3, 1f, 0f), Vector3.forward, Color.blue);
        Debug.DrawRay(MeshUtil.PointOnTriangle(v1, v2, v3, 0f, 1f), Vector3.forward, Color.green);
        Debug.DrawRay(MeshUtil.PointOnTriangle(v1, v2, v3, 1f, 1f), Vector3.forward, Color.magenta);
        Debug.DrawRay(MeshUtil.PointOnTriangle(v1, v2, v3, 0.5f, 0.5f), Vector3.forward, Color.yellow);
    }
}
