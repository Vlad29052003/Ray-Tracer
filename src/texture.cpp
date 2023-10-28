#include "texture.h"
#include "render.h"
#include <framework/image.h>

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    int i = (int)floor(texCoord.x * image.width);
    int j = (int)floor(texCoord.y * image.height);

    return image.pixels[j * image.width + i];
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    //int i = (1.0f - texCoord.x) * image.width;
    //int j = (1.0f - texCoord.y) * image.height;
    auto tempx = (int)floor(texCoord.x * image.width);
    auto tempy = (int)floor(texCoord.y * image.height);

    /* int x0 = static_cast<int>(x);
    int x1 = x0 + 1;
    int y0 = static_cast<int>(y);
    int y1 = y0 + 1;

    // Ensure the coordinates are within the image bounds
    x0 = std::max(0, std::min(x0, image.width - 1));
    x1 = std::max(0, std::min(x1, image.width - 1));
    y0 = std::max(0, std::min(y0, image.height - 1));
    y1 = std::max(0, std::min(y1, image.height - 1));

    // Calculate the interpolation weights
    float sx = x - x0;
    float sy = y - y0;
    float a = (1.0f - sx) * (1.0f - sy);
    float b = sx * (1.0f - sy);
    float c = (1.0f - sx) * sy;
    float d = sx * sy;

    // Sample the texels and interpolate color
    glm::vec3 texel00 = image.pixels[y0 * image.width + x0];
    glm::vec3 texel01 = image.pixels[y0 * image.width + x1];
    glm::vec3 texel10 = image.pixels[y1 * image.width + x0];
    glm::vec3 texel11 = image.pixels[y1 * image.width + x1];

    glm::vec3 interpolatedColor = a * texel00 + b * texel01 + c * texel10 + d * texel11;

    return interpolatedColor;*/
    glm::vec2 v0 = { tempx, tempy };
    glm::vec2 v1 = { tempx+1, tempy };
    glm::vec2 v2 = { tempx+1, tempy+1 };
    glm::vec2 v3 = { tempx, tempy+1 };

    glm::vec2 position = v0 * texCoord.x * texCoord.y + v1 * (1.0f - texCoord.x) * texCoord.y + 
        v2 * (1.0f - texCoord.x) * (1.0f - texCoord.y) + v3 * texCoord.x * (1.0f - texCoord.y);

    auto i = floor(position.x);
    auto j = floor(position.y);
    return image.pixels[j * image.width + i];
}