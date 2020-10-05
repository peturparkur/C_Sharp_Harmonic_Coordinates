using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TensorTest : MonoBehaviour
{
    public Tensor<int> tensor;

    // Start is called before the first frame update
    void Start()
    {
        tensor = new Tensor<int>(new int[] { 3, 4, 5 }); //3x3x3 volume?
        for(int i=0; i<tensor.values.Length; i++)
        {
            tensor.values[i] = i;
        }
        Debug.Log("tensor size = " + tensor.values.Length);


        for(int i=0; i<tensor.dim[0]; i++)
        {
            for(int j=0; j< tensor.dim[1]; j++)
            {
                for(int k=0; k< tensor.dim[2]; k++)
                {
                    Debug.Log("Accessing (i,j,k) = (" + i + "," + j + "," + k + ")");
                    Debug.Log("Tensor value = " + tensor[new int[] { i, j, k }] + " this is the " + (k + (j * tensor.dim[2]) +(i* tensor.dim[2] * tensor.dim[1])) + "th element");
                }
            }
        }

        int[] accessor = new int[] { 1, 2, 3 };
        Debug.Log("accessor1 = (" + accessor[0] + ", " + accessor[1] + ", " + accessor[2] + ")");
        int index = tensor.Index(accessor);
        Debug.Log("index = " + index);
        int[] accessor2 = tensor.IndexToAccessor(index);
        Debug.Log("accessor2 = (" + accessor2[0] + ", " + accessor2[1] + ", " + accessor2[2] + ")");
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
