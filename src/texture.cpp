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
    int i = static_cast<int>(floor(texCoord.x * image.width));
    int j = static_cast<int>(floor(texCoord.y * image.height));

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
    auto tempx = static_cast<int>(floor(texCoord.x * image.width));
    auto tempy = static_cast<int>(floor(texCoord.y * image.height));

    auto alpha = texCoord.x * image.width - tempx;
    auto beta = texCoord.y * image.height - tempy;

    glm::vec2 v0 = { tempx, tempy };
    glm::vec2 v1 = { std::max(0, std::min(tempx+1, image.width)), tempy };
    glm::vec2 v2 = { std::max(0, std::min(tempx+1, image.width)), std::max(0, std::min(tempy+1, image.height)) };
    glm::vec2 v3 = { tempx, std::max(0, std::min(tempy+1, image.height)) };

    return image.pixels[v0.y * image.width + v0.x] * (1.0f - alpha) * (1.0f - beta) + image.pixels[v1.y * image.width + v1.x] * alpha * (1.0f - beta) +
        image.pixels[v2.y * image.width + v2.x] * alpha * beta + image.pixels[v3.y * image.width + v3.x] * (1.0f - alpha) * beta;
}