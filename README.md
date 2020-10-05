# C_Sharp_Harmonic_Coordinates
C# Harmonic Coordinates implementation in Unity Game Engine. Based on the paper https://graphics.pixar.com/library/HarmonicCoordinates/paper.pdf

Current Implementation doesn't use the GPU yet. This could be extended so the diffusion part of the method is done on the GPU as well as the Mesh deformation using the precomputed weights from the Diffusion. This would make it able to precompute more accurate weights in the same time as it currently takes.

The DebugDrawExtension is only used to debug the boundary and weight calculation. Source: https://assetstore.unity.com/packages/tools/debug-drawing-extension-11396
