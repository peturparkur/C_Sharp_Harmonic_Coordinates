# C_Sharp_Harmonic_Coordinates
C# Harmonic Coordinates implementation in Unity Game Engine. Based on the paper https://graphics.pixar.com/library/HarmonicCoordinates/paper.pdf

Current Implementation doesn't use the GPU yet. This could be extended so the diffusion part of the method is done on the GPU as well as the Mesh deformation using the precomputed weights from the Diffusion. This would make it able to precompute more accurate weights.
