# Ray Tracer Project

Welcome to the Ray Tracer project! This ray tracing application is designed to showcase various features and techniques in computer graphics and rendering.

## Features

This ray tracer includes the following features:

- **Acceleration Data-Structure**: Accelerated ray-object intersection using BVH (Bounding Volume Hierarchy) with SAH (Surface Area Heuristic) and binning as the splitting criterion.
- **Shading Models**: Supports various shading models for realistic rendering.
- **Recursive Ray Reflections**: Simulates reflective surfaces with recursive ray tracing.
- **Recursive Ray Transparency**: Achieves transparency effects through recursive ray tracing.
- **Interpolation**: Provides smooth transitions and gradients in rendering.
- **Texture Mapping**: Maps textures onto 3D objects for added realism.
- **Lights and Shadows**: Realistic lighting and shadow effects.
- **Multisampling**: Anti-aliasing technique to reduce aliasing artifacts.
- **Bloom Filter Post-processing**: Adds a bloom effect to the final image.
- **Glossy Reflections**: Simulates glossy reflections and offers visual debugging.
- **Depth of Field**: Realistic depth of field with visual debugging support.

## Getting Started

### Prerequisites

To run the project, you need:

- Visual Studio (in Release mode).
- Graphics card with support for modern OpenGL.
- Input device (keyboard or mouse).

### Building and Running

1. Open the project in Visual Studio.
2. Set the build configuration to Release.
3. Build the project.
4. Run the project.

### User Interface

- The application provides a user interface for enabling/disabling features and adjusting parameters using sliders.
- In rasterization mode, press 'R' to generate a ray (used for debugging transparency, glossy reflections, and reflections).
- In ray tracing mode, enjoy live rendering of the selected scene.
- An option to render the current scene to a file is available.

## Visual Debugging

- Press 'R' to generate a ray for debugging in rasterization mode.
- In-depth debugging support for features like depth of field (press 'D' when enabled).
