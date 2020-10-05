using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ComputeShaderTest : MonoBehaviour
{
    public ComputeShader shader; //compute shader we will use
    public Texture start;


    public RenderTexture texture; //the output of the compute shader
    public RenderTexture inp; //the output of the compute shader

    public int kernel;
    public int kernel2;
    public int kernel3;
    public int[] textureSize = new int[] {256,256};

    public bool ping = true;

    // Start is called before the first frame update
    void Start()
    {
        kernel = shader.FindKernel("CSMain");
        kernel2 = shader.FindKernel("RandomNoise");
        kernel3 = shader.FindKernel("CopyTexture");
        Debug.Log("compute shader kernel = " + kernel);
        texture = new RenderTexture(textureSize[0], textureSize[1], 24); //100x100x1 render texture the last value is the debth buffer
        texture.enableRandomWrite = true;
        //texture.filterMode = FilterMode.Point;
        //texture.wrapMode = TextureWrapMode.Clamp;
        texture.Create();

        shader.SetInts("dim", new int[] { textureSize[0], textureSize[1] }); //dimensions of the texture
        //shader.SetTexture(kernel, "Result", texture); //output texture

        
        //for the input
        inp = new RenderTexture(textureSize[0], textureSize[1], 24); //100x100x1 render texture the last value is the debth buffer
                                                 //inp.enableRandomWrite = true;
                                                 //inp.filterMode = FilterMode.Point;
                                                 //inp.wrapMode = TextureWrapMode.Clamp;
        inp.enableRandomWrite = true;
        inp.Create();
        shader.SetTexture(kernel2, "input", inp);
        shader.SetTexture(kernel2, "Result", texture);

        kernel3 = shader.FindKernel("CopyTexture");
        shader.SetTexture(kernel3, "input", inp);
        shader.SetTexture(kernel3, "Result", texture);
        //Graphics.Blit(start, inp);
        shader.Dispatch(kernel2, textureSize[0] / 8, textureSize[1] / 8, 1); //generating the things
        //shader.SetTexture(kernel, "input", inp);

        int k = shader.FindKernel("Fade");
        shader.SetTexture(k, "Result", texture); //output texture
        shader.SetTexture(k, "input", inp);

    }

    void Swap()
    {
        RenderTexture tmp = inp;
        inp = texture;
        texture = tmp;
    }

    // Update is called once per frame
    void Update()
    {
        if (ping)
        {
            int k = shader.FindKernel("Fade");
            
            shader.SetTexture(k, "Result", texture); //output texture
            shader.SetTexture(k, "input", inp);
            shader.Dispatch(k, textureSize[0] / 8, textureSize[1] / 8, 1);

            //int k2 = shader.FindKernel("CopyTexture");
            //shader.Dispatch(k2, textureSize[0] / 8, textureSize[1] / 8, 1);
            //Swap();
            ping = !ping;
        }
        else
        {
            int k = shader.FindKernel("Fade");
            shader.SetTexture(k, "input", texture); //output texture
            shader.SetTexture(k, "Result", inp);
            shader.Dispatch(k, textureSize[0] / 8, textureSize[1] / 8, 1);

            //int k2 = shader.FindKernel("CopyTexture");
            //shader.Dispatch(k2, textureSize[0] / 8, textureSize[1] / 8, 1);
            //Swap();
            ping = !ping;
        }
        //shader.Dispatch(kernel3, textureSize[0] / 8, textureSize[1] / 8, 1); //copy texture
    }
}
