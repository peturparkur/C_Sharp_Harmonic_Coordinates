using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using PUtil.MeshUtil;
using System.IO; //for saving the weights

public class HarmonicScript : MonoBehaviour
{
    public int numCells = 2;
    public VoxelGrid voxel;
    public Mesh boundary;
    //public MeshCollider meshCollider; //just for the lazy boundary detection, this could be done from own code
    public Mesh targetMesh;
    //line intersecting triangle check https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    //or check http://jcgt.org/published/0005/03/03/

    //ATM the interiors are strictly inside the mesh
    public List<int> interiorAccessors; //everything else is boundary or exterior which we don't care about
    public Matrix<float> weights;
    public int maxIteration = 100; //for testing

    public bool clearWeightFile = false;
    public string path = "Assets/Animation/Weight.txt";

    public bool initComplete = false; //for debug
    public bool weightsSet = false; //for debug
    public int t = 0; //for debug
    public float delta = 1f; //for debug
    public float epsilon = 0.0000001f; //for debug


    VertexCalculation gpuCalculator; //testing
    public ComputeShader shader;


    void Init()
    {

        if (boundary == null)
        {
            boundary = MeshUtil.GenerateTetrahedron3(2f);
            GetComponent<MeshFilter>().mesh = boundary;
            targetMesh = MeshUtil.GenerateTetrahedron3(1.5f);
            GetComponent<MeshFilter>().mesh = targetMesh;
        }
        /*
        if (meshCollider == null)
        {
            meshCollider = GetComponent<MeshCollider>();
            meshCollider.sharedMesh = boundary;
            //meshCollider.enabled = true;
            //meshCollider.sharedMesh = boundary;
            //meshCollider.convex = false;
        }*/

        voxel = new VoxelGrid(numCells); //allocate memory for the voxels
        Vector2[] bounds = MeshBoundingBox(boundary.vertices); //boundaries taken from the mesh vertices
        /*
        DebugExtension.DebugPoint(new Vector3(bounds[0].x, 0,0));
        DebugExtension.DebugPoint(new Vector3(bounds[0].y, 0, 0));
        DebugExtension.DebugPoint(new Vector3(0, bounds[1].x, 0));
        DebugExtension.DebugPoint(new Vector3(0, bounds[1].y, 0));
        DebugExtension.DebugPoint(new Vector3(0, 0, bounds[2].x));
        DebugExtension.DebugPoint(new Vector3(0, 0, bounds[2].y));
        */
        //bounding box defined as Vector2[] where Vector2[0] is the bounds for the X coordinates, vector2[0].x => min, vector2[0].y => max
        voxel.SetBoundaries(bounds); //define the space the voxel occupies
        FindBounds(); //Bakes the voxel
        Debug.Log("interval = " + voxel.d.ToString());
/*
        Debug.Log("lower bounds = (" + voxel.BoundingBox[0, 0] + ", " + voxel.BoundingBox[1, 0] + ", " + voxel.BoundingBox[0, 0] + ")");
        Debug.Log("upper bounds = (" + voxel.BoundingBox[0, 1] + ", " + voxel.BoundingBox[1, 1] + ", " + voxel.BoundingBox[0, 1] + ")");
        Debug.Log("voxel d = " + voxel.d.ToString());
        for (int i=0; i<voxel.grid.dim[1]; i++)
        {
            int[] a = new int[] { 0, i, 0 };
            Vector3 v = voxel.GridPointToVector(a);
            Debug.Log("V[" + i + "] = " + v);

            float y = voxel.BoundingBox[0, 0] + (voxel.BoundingBox[0, 1] - voxel.BoundingBox[0, 0]) * ((float)i / (float)(voxel.grid.dim[1] - 1));
            Debug.Log("Y[" + i + "] = " + y);
        }*/
    }

    Vector2[] MeshBoundingBox(Vector3[] verts)
    {
        Vector2[] bounds = new Vector2[3];
        //initial point set
        bounds[0].x = verts[0].x;
        bounds[0].y = verts[0].x;

        bounds[1].x = verts[0].y;
        bounds[1].y = verts[0].y;

        bounds[2].x = verts[0].z;
        bounds[2].y = verts[0].z;

        for (int i=1; i<verts.Length; i++)
        {
            if (verts[i].x < bounds[0].x) bounds[0].x = verts[i].x - 2f * 0.000001f; //less than X min
            if (verts[i].x > bounds[0].y) bounds[0].y = verts[i].x + 2f * 0.000001f; //more than X max

            if (verts[i].y < bounds[1].x) bounds[1].x = verts[i].y - 2f * 0.000001f; //less than Y min
            if (verts[i].y > bounds[1].y) bounds[1].y = verts[i].y + 2f * 0.000001f; //more than Y max

            if (verts[i].z < bounds[2].x) bounds[2].x = verts[i].z - 2f * 0.000001f; //less than Z min
            if (verts[i].z > bounds[2].y) bounds[2].y = verts[i].z + 2f * 0.000001f; //more than Z max
        }
        return bounds; //returns bounding box
    }

    void FindBounds() //finds the bounds in the voxels for the given boundary mesh, //for now it's a lazy raycast
    {
        //we raycast from the lowest coordinates along the X axis looping across rest
        Debug.Log("Find Bounds");
        for(int i=0; i<voxel.grid.dim[1]; i++) //y looping
        {
            for (int j = 0; j < voxel.grid.dim[2]; j++) //z looping
            {
                Vector3 v = voxel.GridPointToVector(new int[] { 0, i, j }); //origin
                //Debug.Log("V = " + v);

                //RaycastHit[] hits = new RaycastHit[2];
                Vector3[] hits;
                hits = MeshUtil.RaycastMesh(v - 10f * Vector3.right * 0.00001f, Vector3.right, voxel.BoundingBox[0, 1] - voxel.BoundingBox[0, 0] + 10f * 0.00001f, boundary);
                MeshUtil.OrderVectorArray(hits, v - 5f * Vector3.right * 0.00001f);
                //hits = Physics.RaycastAll(v-10f*Vector3.right*Mathf.Epsilon, Vector3.right, voxel.BoundingBox[0, 1] - voxel.BoundingBox[0, 0] + 10f*Mathf.Epsilon); //ray from v towards x direction from end-to-end of voxel box
                //bool m = Physics.Raycast(v - 10f * Vector3.right * Mathf.Epsilon, Vector3.right, out hits[0], (voxel.BoundingBox[0, 1] - voxel.BoundingBox[0, 0]) + 10f * Mathf.Epsilon);
                //if (m) {l= Physics.Raycast(hits[0].point, Vector3.right, out hits[1], (voxel.BoundingBox[0, 1] - hits[0].point.x) + 10f * Mathf.Epsilon); }
                //DebugExtension.DebugArrow(v, Vector3.right * (voxel.BoundingBox[0, 1] - voxel.BoundingBox[0, 0]));
                int n = hits.Length; //if n is 0 that line is empty make everything exterior
                //if n is one -> that's odd because it implies we went into the object but never come out
                //if n is even -> as expected we found the entrance and exit of the boundary, maybe even multiple ones

                if (n == 0)
                {
                    //Debug.Log("No collision happened at (i,j) = (" + i + ", " + j + ")");
                    DebugExtension.DebugPoint(v, Color.red, voxel.d[0]);
                    //fill up the entire row with exteriors
                    for (int k = 0; k < voxel.grid.dim[0]; k++)
                    {
                        int[] a = new int[] { k, i, j };
                        voxel.grid[a] = new Cell(CType.Exterior);
                    }
                    continue;
                };
                if(n%2 == 1) //if n is more than 0 as odd
                {
                    DebugExtension.DebugArrow(v -5f * Vector3.right * Mathf.Epsilon, Vector3.right * (voxel.BoundingBox[0, 1] - voxel.BoundingBox[0, 0] + 5f * Mathf.Epsilon), Color.cyan);
                    for (int m = 0; m < hits.Length; m++)
                    {
                        DebugExtension.DebugPoint(hits[m], Color.red, voxel.d[0]);
                        Debug.Log("collisions at " + hits[m]);
                    }
                    throw new Exception("Weird thing happened at (i,j) = (" + i + ", " + j + ") as there was an odd number of collisions " + n);
                } 

                if(n%2 == 0) //we had an even number of collisions => the shape is closed => that's what we want
                {
                    //if(n>2) throw new Exception("Weird thing happened at (i,j) = (" + i + ", " + j + ") as there was an weird number of collisions " + n);
                    //Debug.Log("At point (" + i + ", " + j + "), raycasts had " + n + " hits");
                    //if (i == 0 && j == 2) Debug.Log("number of hits = " + n);
                    for (int k=0; k<n; k+=2)
                    {
                        Vector3 v1 = hits[k];//.point;
                        Vector3 v2 = hits[k + 1];//.point;
                        //if (i == 0 && j == 2) Debug.Log("hits happened at  " + v1 + ", and " + v2);
                        DebugExtension.DebugPoint(v1, Color.green, voxel.d[0] / 2f);
                        DebugExtension.DebugPoint(v2, Color.cyan, voxel.d[0] / 2f);

                        int[] p = voxel.VectorToGridPoint(v1); //the first hit's voxel position rounded down
                        //if (i == 0 && j == 2) { Debug.Log("hits translated to " + new Vector<int>(3, p).ToString()); Debug.Log("hit is meant to be at (" + p[0] + ", 0, 2)"); }
                        //if(i==1 && j == 2) { Debug.Log("Rounded vector = " + voxel.RoundedVectorToGridPoint(v1) + ", and " + voxel.RoundedVectorToGridPoint(v1)); }
                        voxel.grid[new int[] {p[0], i,j}] = new Cell(CType.Boundary); //this is the cells that are just outside the interior

                        p = voxel.VectorToGridPoint(v2); //the second hit's voxel rounded down
                        //if (i == 0 && j == 2) Debug.Log("hits translated to " + new Vector<int>(3, p).ToString());
                        if (p[0] + 1 < voxel.grid.dim[0])
                        {
                            voxel.grid[new int[] { p[0]+1, i, j }] = new Cell(CType.Boundary); //we need the change the other side of the collision so it is rounded up
                        }
                        //Debug.Log("position (" + p[0] + "," + p[1] + "," + p[2] + ") set to boundary");
                    }
                }


                //filling up the voxel with the exteriors and interiors
                bool bound = false; //have we found the boundary yet
                //We could fill up the voxels here with the exteriors
                for(int k=0; k<voxel.grid.dim[0]; k++)
                {
                    int[] a = new int[] { k, i, j };
                    CType type = voxel.grid[a].type;
                    if(type == CType.Boundary) { bound = !bound; continue; } //we flip the bound
                    //we can still have 2 boundaries next to each other which should result in only exteriors being next to them
                    if (type == CType.Undefined)
                    {
                        if(!bound) //we haven't found the bound yet thus we're on the exterior
                        {
                            voxel.grid[a] = new Cell(CType.Exterior);
                        }
                        else
                        {
                            voxel.grid[a] = new Cell(CType.Interior);
                            interiorAccessors.Add(voxel.grid.Index(a)); //converts the accessor to an index
                        }
                    }
                }
            }
        }

        Debug.Log("We have " + interiorAccessors.Count + " interior points");
        
    }

    public int[] FindNearestInteriors(Vector3 v, float delta=0.000001f)
    {
        //we want to find the nearest interior points to v
        //these are just the initial values
        float min = (v - voxel.GridPointToVector(voxel.grid.IndexToAccessor(interiorAccessors[0]))).sqrMagnitude;
        List<int> nearest = new List<int>();
        nearest.Add(interiorAccessors[0]);

        for(int i=1; i<interiorAccessors.Count; i++)
        {
            float d = (v - voxel.GridPointToVector(voxel.grid.IndexToAccessor(interiorAccessors[i]))).sqrMagnitude;
            if(d < min) //if is less than min then it's better
            {
                nearest.Clear();
                nearest.Add(i);
                min = d;
                continue;
            }
            else
            {
                if(d<min + delta) //is the distance almost equal to d
                {
                    //then i is in acceptable range of the current distance
                    nearest.Add(i);
                    continue;
                }
            }
        }
        return nearest.ToArray(); //the equal lengthed nearest points to v
    }

    public void RenderVoxels(float scale = 1f, bool drawInterior=true, bool drawBoundary=true, bool drawExterior=false)
    {
        for(int i=0; i<voxel.grid.dim[0]; i++)
        {
            for (int j = 0; j < voxel.grid.dim[1]; j++)
            {
                for (int k = 0; k < voxel.grid.dim[2]; k++)
                {
                    int[] a = new int[] { i, j, k };
                    Vector3 v = voxel.GridPointToVector(a);
                    Color c = Color.clear;
                    float tmpScale = 1f;
                    if (voxel.grid[a].type == CType.Exterior) { c = Color.white; if (!drawExterior) { tmpScale = 0f; } }
                    if (voxel.grid[a].type == CType.Boundary) { c = Color.black; if (!drawBoundary) { tmpScale = 0f; } }
                    if (voxel.grid[a].type == CType.Interior) { c = Color.blue; if (!drawInterior) { tmpScale = 0f; } }

                    DebugExtension.DebugPoint(v, c, (tmpScale*scale) * (voxel.d[0] / 2f));

                }
            }
        }
    }

    // Start is called before the first frame update
    void Start()
    {

        if (File.Exists(path) && !clearWeightFile) { return; }
        if(File.Exists(path) && clearWeightFile) { File.WriteAllText(path, ""); }
        StreamWriter stream = File.CreateText(path);

        Init();
        initComplete = false;
        delta = 1f;
        CreateWeights(boundary, targetMesh);
        Debug.Log("start baking");
        float time = Time.realtimeSinceStartup;
        weights = BakeWeights(boundary, targetMesh, maxIteration);
        Debug.Log("finished baking after time = " + (Time.realtimeSinceStartup-time));
        Debug.Log("Weights = " + weights);
        stream.WriteLine(weights.rows + ", " + weights.columns);
        stream.WriteLine(weights.ToString());

        /*
        gpuCalculator = new VertexCalculation();
        gpuCalculator.cs = shader;
        gpuCalculator.SetBuffers(weights, boundary.vertices, targetMesh.vertices);
        */
    }

    // Update is called once per frame
    void Update()
    {
        //testing GPU
        Vector3[] vec = boundary.vertices;
        vec[0] += boundary.normals[0];
        boundary.vertices = vec;


        Vector3[] newVerts = new Vector3[targetMesh.vertices.Length];
        //Debug.Log("new verts = " + new Vector<Vector3>(newVerts.Length, newVerts).ToString());
        for (int i=0; i<weights.rows; i++)
        {
            for (int j = 0; j < weights.columns; j++)
            {
                //Debug.Log("at [" + i + "," + j + "] = " + boundary.vertices[j]);
                //Debug.Log("at [" + i + "," + j + "] = " + weights[i, j]);
                //Debug.Log("at [" + i + "," + j + "] = " + (weights[i, j]*boundary.vertices[j]));
                //Debug.Log("at [" + i + "," + j + "] = " + (weights[i, j] * boundary.vertices[j]).x);
                newVerts[i] += weights[i, j] * boundary.vertices[j];
            }
        }

        Debug.Log("boundary vertices = " + new Vector<Vector3>(boundary.vertices.Length, boundary.vertices).ToString());
        Debug.Log("old verts = " + new Vector<Vector3>(newVerts.Length, targetMesh.vertices).ToString());
        Debug.Log("new verts = " + new Vector<Vector3>(newVerts.Length,newVerts).ToString());
        targetMesh.vertices = newVerts;
        
        //MeshFilter filter = GetComponent<MeshFilter>();
        //if (filter.mesh != targetMesh) { filter.mesh = targetMesh; }

        /*
        int[] near = new int[0];
        float[] vals = new float[0];

        //RenderVoxels(0f, false, false, false);
        if (!initComplete)
        {
            near = InitSmoothing(boundary.vertices[0]);
            vals = new float[near.Length];
            for (int i = 0; i < vals.Length; i++)
            {
                vals[i] = 1f;
            }
            initComplete = true;
            Debug.Log("init complete + " + near.Length);
        }

        RenderVoxels(1f, false, false, false);
        DebugExtension.DebugPoint(boundary.vertices[0], 1f);
        RenderInterior(Color.blue, Color.red);
        if (t > maxIteration || delta < epsilon)
        {
            if (weightsSet) return;

            float sum = 0f;
            for(int i=0; i<targetMesh.vertices.Length; i++)
            {
                //this should be trilinearly interpolated
                weights[i, 0] = voxel.gridValues[voxel.VectorToGridPoint(targetMesh.vertices[i])]; //this is just the value of the nearest gridpoint
                sum += weights[i, 0];
            }
            //this is to normalise the weights
            for (int i = 0; i < weights.rows; i++)
            {
                weights[i, 0] /= sum;
            }
            weightsSet = true;
            Debug.Log("Weights for the first boundary vertex = " + weights.ToString());
                return;
        }
        delta = SmoothingStep(near, vals);
        t++;
        //after smoothstep we should have a nicely smoothed representation weights correlated to the boundary vertex
        //we need to set the weights "Matrix" to by the nearest trilinearly interpolated values of the voxels
        Debug.Log("Completed step T = " + t);
        */
    }

    Matrix<float> BakeWeights(Mesh bound, Mesh target, int maxIter)
    {
        Matrix<float> w = new Matrix<float>(target.vertices.Length, bound.vertices.Length);
        for(int i=0; i<bound.vertices.Length; i++)
        {
            //setup initial conditions
            int[] nr = InitSmoothing(bound.vertices[i]);
            for (int j = 0; j < nr.Length; j++)
            {
                voxel.gridValues.values[interiorAccessors[nr[j]]] = 1f; //set the near points to 1
            }

            for (int t=0; t<maxIter; t++)
            {
                SmoothingStep(new int[] { }, new float[] { });
            }
            for(int j=0; j<w.rows; j++)
            {
                Vector<int> vec = new Vector<int>(3, voxel.VectorToGridPoint(target.vertices[j]));
                Debug.Log("vec at j=" + j + " => " + vec.ToString() + " where the vertex is at point " + target.vertices[j]);
                w[j, i] = voxel.gridValues[voxel.VectorToGridPoint(target.vertices[j])]; //the nearest voxel point for now to vertex j
                //sum += w[j, i];
            }
        }
        //Then we normalise rows
        for(int i=0; i<w.rows; i++)
        {
            float sum = 0f;
            for(int j=0; j<w.columns; j++)
            {
                sum += w[i, j];
            }
            for (int j = 0; j < w.columns; j++)
            {
                w[i,j] /= sum;
            }

        }


        return w;
    }

    void RenderInterior(Color low, Color high)
    {
        float maxVal = voxel.gridValues.values[interiorAccessors[0]];
        //Debug.Log("MaxVal = " + maxVal);
        for(int i=1; i<interiorAccessors.Count; i++)
        {
            float val = voxel.gridValues.values[interiorAccessors[i]];
            if (val > maxVal)
            {
                maxVal = val;
            }
        }

        for(int i=0; i<interiorAccessors.Count; i++)
        {
            int[] a = voxel.grid.IndexToAccessor(interiorAccessors[i]);
            float val = voxel.gridValues.values[interiorAccessors[i]];
            float scale = voxel.d[0];
            Color c = Color.Lerp(low, high, val/maxVal);
            if (val / maxVal < 0.1f) continue;
            DebugExtension.DebugPoint(voxel.GridPointToVector(a), c, scale);
        }
    }

    int[] InitSmoothing(Vector3 vertex)
    {
        for(int i=0; i<voxel.gridValues.values.Length; i++) { voxel.gridValues.values[i] = 0f; }//we reset all values to 0
        return FindNearestInteriors(vertex); //this array is in terms of the interior accessors
    }

    float SmoothingStep(int[] index, float[] initValues) //quick CPU implementation of the smoothing for now
    {
        //for the given indexes we keep the initial values
        for(int i=0; i<index.Length; i++)
        {
            voxel.gridValues.values[interiorAccessors[index[i]]] = initValues[i];
        }
        float maxChange = 0f;

        Tensor<float> tmpGrid = new Tensor<float>(voxel.gridValues); //creates a copy of the voxel gridvalues
        for(int i=0; i<interiorAccessors.Count; i++)
        {
            //Debug.Log("completed " + i);
            int[] a = voxel.gridValues.IndexToAccessor(interiorAccessors[i]); //accessor
            //now we want to access the neighbouring values
            int count = 0;
            float average = voxel.gridValues.values[interiorAccessors[i]];

            int[] r1;
            if (a[0] + 1 < voxel.grid.dim[0])
            {
                r1 = new int[] { a[0] + 1, a[1], a[2] };
                if (CType.Interior == voxel.grid[r1].type) { count++; average += voxel.gridValues[r1]; }
            }
            if (a[0] - 1 < voxel.grid.dim[0])
            {
                r1 = new int[] { a[0] - 1, a[1], a[2] };
                if (CType.Interior == voxel.grid[r1].type) { count++; average += voxel.gridValues[r1]; }
            }
            if (a[1] + 1 < voxel.grid.dim[0])
            {
                r1 = new int[] { a[0], a[1] + 1, a[2] };
                if (CType.Interior == voxel.grid[r1].type) { count++; average += voxel.gridValues[r1]; }
            }
            if (a[1] - 1 < voxel.grid.dim[0])
            {
                r1 = new int[] { a[0], a[1] - 1, a[2] };
                if (CType.Interior == voxel.grid[r1].type) { count++; average += voxel.gridValues[r1]; }
            }
            if (a[2] + 1 < voxel.grid.dim[0])
            {
                r1 = new int[] { a[0], a[1], a[2] + 1 };
                if (CType.Interior == voxel.grid[r1].type) { count++; average += voxel.gridValues[r1]; }
            }
            if (a[2] - 1 < voxel.grid.dim[0])
            {
                r1 = new int[] { a[0], a[1], a[2] - 1 };
                if (CType.Interior == voxel.grid[r1].type) { count++; average += voxel.gridValues[r1]; }
            }
            if (count == 0)
            {
                DebugExtension.DebugPoint(voxel.GridPointToVector(a), Color.red, voxel.d[0]);
                DebugExtension.DebugPoint(voxel.GridPointToVector(new int[] { a[0] + 1, a[1], a[2] }), voxel.d[0]);
                DebugExtension.DebugPoint(voxel.GridPointToVector(new int[] { a[0] - 1, a[1], a[2] }), voxel.d[0]);
                DebugExtension.DebugPoint(voxel.GridPointToVector(new int[] { a[0], a[1] + 1, a[2] }), voxel.d[0]);
                DebugExtension.DebugPoint(voxel.GridPointToVector(new int[] { a[0], a[1] - 1, a[2] }), voxel.d[0]);
                DebugExtension.DebugPoint(voxel.GridPointToVector(new int[] { a[0], a[1], a[2] + 1 }), voxel.d[0]);
                DebugExtension.DebugPoint(voxel.GridPointToVector(new int[] { a[0], a[1], a[2] - 1 }), voxel.d[0]);
                throw new Exception("We found no neighbours for " + new Vector<int>(3, a).ToString() + " at index " + i);
            }

            tmpGrid[a] = average / (float)(count+1);
            float change = Mathf.Abs(tmpGrid[a] - voxel.gridValues[a]);
            if (change > maxChange) maxChange = change;
        }
        voxel.gridValues = new Tensor<float>(tmpGrid);
        return maxChange;
    }

    public void CreateWeights(Mesh bound, Mesh targetMesh)
    {
        weights = new Matrix<float>(targetMesh.vertices.Length, bound.vertices.Length); //the weights for vertex i can be accessed as weights[i,x] where x is the boundary vertex
    }
}

public class VoxelGrid
{
    public int numCells; //this is the number of cells that we create within the bounding box - require more than 2
    public Tensor<Cell> grid;
    public Tensor<float> gridValues;
    public Matrix<float> BoundingBox = new Matrix<float>(3, 2); //3 dimension and min-max
    public Vector<float> d = new Vector<float>(3, true); //this defines the space difference between each cell

    public VoxelGrid(int _numCells)
    {
        numCells = _numCells;
        grid = new Tensor<Cell>(new int[] { _numCells, _numCells, _numCells }); //create a NxNxN tensor which represent the voxelisation of the 3D space
        gridValues = new Tensor<float>(new int[] { _numCells, _numCells, _numCells }); //same here
        for(int i=0; i<grid.values.Length; i++) //base values defined
        {
            grid.values[i] = new Cell(CType.Undefined);
            gridValues.values[i] = 0f;
        }
    }
    public void SetBoundaries(Vector2[] box) //x,y,z min-max values for bounding box
    {
        for(int i=0; i<3; i++)
        {
            BoundingBox[i, 0] = box[i][0];
            BoundingBox[i, 1] = box[i][1];
        }

        for(int i=0; i<3; i++)
        {
            d[i] = (BoundingBox[i, 1] - BoundingBox[i, 0]) / (float)(numCells - 1); //the space difference between the voxels
        }
    }

    public Vector3 GridPointToVector(int[] p) //3 length int array
    {
        //we assume that the grids represent the corners of the space
        //0,0,0 represent the minimum value of the bounding box, n-1,n-1,n-1 represents the maximum in the bounding box
        //everything else linearly interpolates
        float x = BoundingBox[0, 0] + (d[0] * (float)p[0]); //technically we wouldn't need d but it makes things easier
        float y = BoundingBox[1, 0] + (d[1] * (float)p[1]);
        float z = BoundingBox[2, 0] + (d[2] * (float)p[2]);
        return new Vector3(x,y,z);

    }

    public int[] VectorToGridPoint(Vector3 v)
    {
        int x = (int)((v.x - BoundingBox[0, 0]) / d[0]); //the distance from the x minimum
        int y = (int)((v.y - BoundingBox[1, 0]) / d[1]); //the distance from the x minimum
        int z = (int)((v.z - BoundingBox[2, 0]) / d[2]); //the distance from the x minimum

        return new int[] { x, y, z };
    }

    public Vector3 RoundedVectorToGridPoint(Vector3 v)
    {
        int x = (int)((v.x - BoundingBox[0, 0]) / d[0]); //the distance from the x minimum
        int y = (int)((v.y - BoundingBox[1, 0]) / d[1]); //the distance from the x minimum
        int z = (int)((v.z - BoundingBox[2, 0]) / d[2]); //the distance from the x minimum

        return new Vector3(x, y, z);
    }

}

public struct Cell //the voxel cells
{
    public CType type; //-1 untyped, 0 Boundary, 1 Exterior, 2 Interior
    //public double value; //this represents the weight

    public Cell(CType _type)//, double _val)
    {
        type = _type;
        //value = _val;
    }
}
