#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

}

Screen& getBrightAreas(Screen& image) {
    Screen output(image.resolution());
    for (auto y = 0; y < output.resolution().y; y++) {
        for (auto x = 0; x < output.resolution().x; x++) {
            auto color = output.pixels()[output.indexAt(x, y)];
            auto perceived_luminance = 0.2126f * color.r + 0.7152f * color.g + 0.0722f * color.b; //https://en.wikipedia.org/wiki/Relative_luminance

            if (perceived_luminance > 0.9f)
                output.pixels()[output.indexAt(x, y)] = { 1.0f, 1.0f, 1.0f };
            else
                output.pixels()[output.indexAt(x, y)] = { 0.0f, 0.0f, 0.0f };
        }
    }

    return output;
}

float boxFilter(Screen& image, int i, int j, int col, int filterSize)
{
    filterSize = std::max(1, filterSize);
    float sum = 0;
    // Average pixels in the box-filter region
    for (int x = -filterSize; x < filterSize + 1; x++)
        for (int y = -filterSize; y < filterSize + 1; y++)
            sum += image.pixels()[image.indexAt(i + x, j + y)][col];
    sum /= (2 * filterSize + 1) * (2 * filterSize + 1);
    return sum;
}

void filterImage(Screen& image, int filterSize)
{
    // we create a result image
    Screen result(image);
    // and process every channel of every pixel independently
    for (int i = 0; i < image.resolution().x; i++)
        for (int j = 0; j < image.resolution().y; j++) // for each pixel
            for (int col = 0; col < 3; col++) // for each color channel
                result.pixels()[result.indexAt(i, j)][col] = boxFilter(image, i, j, col, filterSize);
}

void combineFilteredWithOriginalImage(const Scene& scene, const Trackball& camera, Screen& originalImage, Screen& filteredImage)
{
    for (int y = 0; y < originalImage.resolution().y; y++) {
        for (int x = 0; x < originalImage.resolution().x; x++) {
            auto filteredPixel = filteredImage.pixels()[filteredImage.indexAt(x, y)];
            auto originalPixel = originalImage.pixels()[originalImage.indexAt(x, y)];

            auto bloomPixel = glm::clamp(originalPixel + filteredPixel, 0.0f, 1.0f);

            originalImage.pixels()[originalImage.indexAt(x, y)] = bloomPixel;
        }
    }
}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (features.extra.enableBloomEffect) {
        Screen& highValues = getBrightAreas(image);

        int filterSize = (image.resolution().x - 1) / 2;
        filterImage(highValues, filterSize);  // Implementation of this helper function got from the class slides
        combineFilteredWithOriginalImage(scene, camera, image, highValues);
    }

    // ...
}

// Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...

    Ray r = generateReflectionRay(ray, hitInfo);

    //generate orthogonal vectors u and v
    //we can generate the first othogonal vector manually and the second using cross product
    glm::vec3 u = glm::vec3(0), v = glm::vec3(0);
    if (r.direction.x != 0) {
        u.z = -r.direction.x;
        u.x = r.direction.z;
    } else if (r.direction.y != 0) {
        u.z = -r.direction.y;
        u.y = r.direction.z;
    } else {
        u.z = r.direction.y;
        u.y = -r.direction.z;
    }

    v = glm::cross(u, ray.direction);

    u = glm::normalize(u);
    v = glm::normalize(v);

    //calculate the regulation factor that will be applied to the initial radius of 1
    float regulationFactor = hitInfo.material.shininess / 64.0f;

    glm::vec3 sumOfInterference = glm::vec3(0);

    for (int i = 0; i < state.features.extra.numGlossySamples; ++i) {
        //map a uniformly distributed 2d sample to the coordinates of a circle, using the regulation factor
        glm::vec2 sample = state.sampler.next_2d();
        float x = regulationFactor * glm::cos(glm::radians(360.f * sample.x));
        float y = regulationFactor * glm::sin(glm::radians(360.f * sample.y));

        glm::vec3 glossyReflectionRayDirection = r.direction + x * u + y * v;
        Ray glossyReflectionRay = Ray(r.origin, glossyReflectionRayDirection, std::numeric_limits<float>::max());

        //sum up the conttributions of each ray
        sumOfInterference += renderRay(state, glossyReflectionRay, rayDepth + 1) * hitInfo.material.ks;
    }

    //hormalize the sum of the rays and add the result to the current color
    sumOfInterference /= state.features.extra.numGlossySamples;
    hitColor += sumOfInterference;
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {
        // Part of your implementation should go here
        return glm::vec3(0.f);
    } else {
        return glm::vec3(0.f);
    }
}


// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// - aabb;       the axis-aligned bounding box around the given triangle set
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires splitting
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    using Primitive = BVH::Primitive;

    return 0; // This is clearly not the solution
}