#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>

// Extra feature
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

    //uniform number generation in the interval [-0.5, 0.5) for generating position on the lens
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<float> distr(-0.5f, 0.5f);


#ifdef NDEBUG // Enable multi threading in Release mode
#pragma omp parallel for schedule(guided)
#endif
    //iterate through all pixels
    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x < screen.resolution().x; x++) {

            glm::vec3 L = glm::vec3(0); //accumulator
            glm::vec2 position = (glm::vec2(x, y) + 0.5f) / glm::vec2(screen.resolution()) * 2.f - 1.f; // position on the screen
            Ray ray = camera.generateRay(position); // the ray from the pixel center
            glm::vec3 focusPoint = ray.origin + features.extra.focusDistance * ray.direction; // calculate the focus point position
            
            for (int samp = 0; samp < features.extra.numDofSamples; samp++) {
                //generate random sample -> represents the position on the square lens
                float xLens = features.extra.lensLength * distr(generator);
                float yLens = features.extra.lensLength * distr(generator);

                //calculates the depth of field ray - which has its origin on the lens and directed towards the focus point
                Ray dofRay;
                dofRay.origin = ray.origin +  xLens * glm::normalize(camera.left()) + yLens * glm::normalize(camera.up()); //shift the origin of the ray on the lens
                dofRay.direction = glm::normalize(focusPoint - dofRay.origin); // calculate the direction (towards the focal point)
                dofRay.t = std::numeric_limits<float>::max();
                RenderState state = {
                    .scene = scene,
                    .features = features,
                    .bvh = bvh,
                    .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
                };
                L += renderRay(state, dofRay, 0); //accumulate colors of the hit point of the rays
            }

            L /= features.extra.numDofSamples; //normalize the color
            screen.setPixel(x, y, L); //set the correcponding color to the pixel
        }
    }
}

std::vector<Ray> generateDofRaysForDebug(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, glm::ivec2 screenResolution, const glm::vec2& pixel, glm::vec3& focusPoint)
{
    //vector of the depth of field rays
    std::vector<Ray> dofRays;

    //random number generators for generating the position on the lens
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<float> distr(-0.5f, 0.5f);

    glm::vec2 position = (glm::vec2(pixel.x, pixel.y) + 0.5f) / glm::vec2(screenResolution) * 2.f - 1.f; // position on the screen
    Ray ray = camera.generateRay(position); // the ray from the pixel center
    focusPoint = ray.origin + features.extra.focusDistance * ray.direction; // calculate the focus point position
    for (int samp = 0; samp < features.extra.numDofSamples; samp++) {
        // generate random sample -> represents the position on the square lens
        float xLens = features.extra.lensLength * distr(generator);
        float yLens = features.extra.lensLength * distr(generator);

        //calculates the depth of field ray - which has its origin on the lens and directed towards the focus point
        Ray dofRay;
        dofRay.origin = ray.origin + xLens * glm::normalize(camera.left()) + yLens * glm::normalize(camera.up()); // shift the origin of the ray on the lens
        dofRay.direction = glm::normalize(focusPoint - dofRay.origin); // calculate the direction (towards the focal point)
        dofRay.t = std::numeric_limits<float>::max();
        dofRays.push_back(dofRay);
    }

    return dofRays;
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

void getBrightAreas(Screen& image, Screen& output) {
    for (auto y = 0; y < output.resolution().y; y++) {
        for (auto x = 0; x < output.resolution().x; x++) {
            auto color = output.pixels()[output.indexAt(x, y)];
            auto perceived_luminance = 0.2126f * color.r + 0.7152f * color.g + 0.0722f * color.b; //https://en.wikipedia.org/wiki/Relative_luminance

            if (perceived_luminance > 0.63) {
                output.pixels()[output.indexAt(x, y)] = { 1.0f, 1.0f, 1.0f };
            } else {
                output.pixels()[output.indexAt(x, y)] = { 0.0f, 0.0f, 0.0f };
            }
        }
    }
}

void applyGaussianFilterHorizontal(Screen& image, const std::vector<double>& filter, int filterSize)
{
    int width = image.resolution().x;
    int height = image.resolution().y;
    int radius = filterSize / 2;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            glm::vec3 color;
            for (int col = 0; col < 3; col++) {
                float sum = 0.0;
                for (int k = -radius; k <= radius; k++) {
                    int idx = glm::clamp(x + k, 0, width - 1);
                    sum += image.pixels()[image.indexAt(idx, y)][col] * filter[k + radius];
                }
                color[col] = sum;
            }

            image.setPixel(x, y, color);
        }
    }
}

void applyGaussianFilterVertical(Screen& image, const std::vector<double>& filter, int filterSize)
{
    int width = image.resolution().x;
    int height = image.resolution().y;
    int radius = filterSize / 2;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            glm::vec3 color(0.0f);
            for (int col = 0; col < 3; col++) {
                float sum = 0.0;
                for (int k = -radius; k <= radius; k++) {
                    int idx = glm::clamp(y + k, 0, height - 1);
                    sum += image.pixels()[image.indexAt(x, idx)][col] * filter[k + radius];
                }
                color[col] = sum;
            }

            image.setPixel(x, y, color);
        }
    }
}

unsigned long long binomialCoefficient(int n, int k)
{
    if (k < 0 || k > n)
        return 0;

    std::vector<unsigned long long> dp(k + 1);
    dp.at(0) = 1; // C(n, 0) = 1

    for (int i = 1; i <= n; i++) {
        for (int j = std::min(i, k); j > 0; j--) {
            dp.at(j) += dp.at(j - 1);
        }
    }

    return dp.at(k);
}

// Function to apply a 2D Gaussian filter to an image
void gaussianFilter(Screen& image, int filterSize, float sigma)
{
    std::vector<double> filter(filterSize);
    float sum = 0.0f;

    for (int k = 0; k < filterSize; k++) {
        int n = filterSize;        

        filter[k] = static_cast<double>(binomialCoefficient(n, k));
        sum += filter[k];
    }

    for (int k = 0; k < filterSize; k++) {
        filter[k] /= sum;
    }

    applyGaussianFilterHorizontal(image, filter, filterSize);
    applyGaussianFilterVertical(image, filter, filterSize);
}

void combineFilteredWithOriginalImage(const Scene& scene, const Trackball& camera, Screen& originalImage, Screen& filteredImage)
{
    for (int y = 0; y < originalImage.resolution().y; y++) {
        for (int x = 0; x < originalImage.resolution().x; x++) {
            auto filteredPixel = filteredImage.pixels()[filteredImage.indexAt(x, y)];
            auto originalPixel = originalImage.pixels()[originalImage.indexAt(x, y)];

            auto bloomPixel = glm::clamp(originalPixel + filteredPixel, 0.0f, 1.0f);

            originalImage.setPixel(x, y, bloomPixel);
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
        Screen filteredImage(image);
        getBrightAreas(image, filteredImage);
        gaussianFilter(filteredImage, 100, 4.0f);  // Implementation of this helper function got from the class slides
        combineFilteredWithOriginalImage(scene, camera, image, filteredImage);
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
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<float> distrRadius(0.0f, 0.5f);
    std::uniform_real_distribution<float> distrAngle(0.0f, 360.0f);

    Ray reflectedRay = generateReflectionRay(ray, hitInfo);

    // Generate orthogonal vectors u and v
    glm::vec3 u, v;
    u = glm::normalize(glm::cross(reflectedRay.direction, reflectedRay.direction + glm::vec3(0, 1, 0)));
    v = glm::normalize(glm::cross(u, reflectedRay.direction)); 

    // Calculate the regulation factor that will be applied to the initial radius of 1
    float regulationFactor = hitInfo.material.shininess / 64.0f;

    glm::vec3 sumOfInterference = glm::vec3(0);

    for (int i = 0; i < state.features.extra.numGlossySamples; ++i) {
        // Map a uniformly distributed angle to coordinates in disk with radius 0.5 * regulationFactor
        // The points are on a circle whose radius is smaller than 0.5 * regulation factor, which is randomly generated.
        float angle = distrAngle(generator);
        float radius = distrRadius(generator);
        float x = glm::cos(angle) * radius * regulationFactor;
        float y = glm::sin(angle) * radius * regulationFactor;

        // Shift the direction of the perfect reflection
        glm::vec3 glossyReflectionRayDirection = reflectedRay.direction + x * u + y * v;
        Ray glossyReflectionRay = Ray(reflectedRay.origin, glossyReflectionRayDirection, std::numeric_limits<float>::max());

        // Sum up the conttributions of each ray
        sumOfInterference += renderRay(state, glossyReflectionRay, rayDepth + 1) * hitInfo.material.ks;
    }
    // Normalize the sum of the rays and add the result to the current color
    sumOfInterference /= state.features.extra.numGlossySamples;
    hitColor += sumOfInterference;
}

std::vector<glm::vec3> glossyDebug(Ray ray, const HitInfo& hitInfo) {
    std::vector<glm::vec3> vectors = std::vector<glm::vec3>();
    // Calculates the perfect reflected Ray
    Ray reflectedRay = generateReflectionRay(ray, hitInfo);

    // Calculates the othonormal basis for the disk
    glm::vec3 u, v;
    u = glm::normalize(glm::cross(reflectedRay.direction, reflectedRay.direction + glm::vec3(0, 1, 0)));
    v = glm::normalize(glm::cross(u, reflectedRay.direction));

    // Calculates the radius which is 0.5 * regulationFactor so shininess / 128
    float radius = hitInfo.material.shininess / 128.0f;
    glm::vec3 rad = glm::vec3(radius);

    vectors.push_back(reflectedRay.origin);
    vectors.push_back(u);
    vectors.push_back(v);
    vectors.push_back(rad);
    return vectors;
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


float areaOfAABB(AxisAlignedBox aabb)
{
    float l1 = abs(aabb.lower[0] - aabb.upper[0]);
    float l2 = abs(aabb.lower[1] - aabb.upper[1]);
    float l3 = abs(aabb.lower[2] - aabb.upper[2]);

    return 2 * (l1 * l2 + l2 * l3 + l3 * l1);
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

    float global_cost = std::numeric_limits<float>::infinity();
    float left_size = 0;

    std::sort(primitives.begin(), primitives.end(),
        [axis](const Primitive& a, const Primitive& b) {
            glm::vec3 centroidA = computePrimitiveCentroid(a);
            glm::vec3 centroidB = computePrimitiveCentroid(b);
            return centroidA[axis] < centroidB[axis];});

    uint32_t num_bins = 10;

    float size_of_bin = (aabb.upper[axis] - aabb.lower[axis]) / num_bins;

    for(int i = 1; i < num_bins; ++i) {
        float dist = aabb.lower[axis] + i * size_of_bin;

        std::vector<BVH::Primitive> left_primitives;
        std::vector<BVH::Primitive> right_primitives;

        for(auto& primitive : primitives) {
            if(computePrimitiveCentroid(primitive)[axis] < dist) {
                left_primitives.push_back(primitive);
            }
            else {
                right_primitives.push_back(primitive);
            }
        }

        if(left_primitives.size() == 0 || right_primitives.size() == 0) continue;

        std::span<BVH::Primitive> left_span = left_primitives;
        AxisAlignedBox left_aabb = computeSpanAABB(left_span);
        std::span<BVH::Primitive> right_span = right_primitives;
        AxisAlignedBox right_aabb = computeSpanAABB(right_span);

        float cost = left_span.size() * areaOfAABB(left_aabb) + right_span.size() * areaOfAABB(right_aabb);

        if(cost < global_cost) {
            global_cost = cost;
            left_size = left_span.size();
        }

    }

    if(left_size == 0) {
        return splitPrimitivesByMedian(aabb, axis, primitives);
    }

    return left_size;
}