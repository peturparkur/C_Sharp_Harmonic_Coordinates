using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using PUtil.GenericMath;

public class InterpolationTestScript : MonoBehaviour
{
    // Start is called before the first frame update
    void Start()
    {
        Debug.Log("Testing interpolation");
        Debug.Log("value at t=0 => " + GenericMath.Cubic(-1f, 1f, 0f));
        Debug.Log("value at t=1 => " + GenericMath.Cubic(-1f, 1f, 1f));
        Debug.Log("value at t=0.5 => " + GenericMath.Cubic(-1f, 1f, 0.5f));
        Debug.Log("value at t=0.1 => " + GenericMath.Cubic(-1f, 1f, 0.1f));

    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
