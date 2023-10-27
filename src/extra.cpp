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
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(-0.5f, 0.5f);


#ifdef NDEBUG // Enable multi threading in Release mode
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x < screen.resolution().x; x++) {
            glm::vec3 L = glm::vec3(0);
            for (int samp = 0; samp < features.extra.numDofSamples; samp++) {
                //generate random sample -> represents the position on the square lens
                float xLens = features.extra.aperture * distribution(generator);
                float yLens = features.extra.aperture * distribution(generator);

                glm::vec2 position = (glm::vec2(x, y) + 0.5f) / glm::vec2(screen.resolution()) * 2.f - 1.f; //position on the screen
                Ray ray = camera.generateRay(position); //the ray from the pixel center
                glm::vec3 focalPoint = ray.origin + features.extra.focusDistance * ray.direction; //calculate the focal point
                ray.origin += xLens * glm::normalize(camera.left()) + yLens * glm::normalize(camera.up()); //shift the origin of the ray on the lens aperture
                ray.direction = glm::normalize(focalPoint - ray.origin); // calculate the direction (towards the focal point)
                RenderState state = {
                    .scene = scene,
                    .features = features,
                    .bvh = bvh,
                    .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
                };
                L += renderRay(state, ray, 0);
            }

            L /= features.extra.numDofSamples;
            screen.setPixel(x, y, L);
        }
    }
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

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
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

    Ray reflectedRay = generateReflectionRay(ray, hitInfo);

    //generate orthogonal vectors u and v
    //we can generate the first othogonal vector manually and the second using cross product
    glm::vec3 u = glm::vec3(0), v = glm::vec3(0);
    if (reflectedRay.direction.x != 0) {
        u.z = -reflectedRay.direction.x;
        u.x = reflectedRay.direction.z;
    } else if (reflectedRay.direction.y != 0) {
        u.z = -reflectedRay.direction.y;
        u.y = reflectedRay.direction.z;
    } else {
        u.z = reflectedRay.direction.y;
        u.y = -reflectedRay.direction.z;
    }
    u = glm::normalize(u);

    v = glm::cross(u, ray.direction);
    v = glm::normalize(v);

    //calculate the regulation factor that will be applied to the initial radius of 1
    float regulationFactor = hitInfo.material.shininess / 64.0f;

    glm::vec3 sumOfInterference = glm::vec3(0);

    for (int i = 0; i < state.features.extra.numGlossySamples; ++i) {
        //map a uniformly distributed 2d sample (essentially a square) into coordinates of a disc
        //using FG-Squircular mapping (https://arxiv.org/ftp/arxiv/papers/1709/1709.07875.pdf)
        glm::vec2 sample = state.sampler.next_2d() - glm::vec2(0.5f, 0.5f); // [0,1) to [-0.5, 0.5)
        float sqRoot = sqrtf(sample.x * sample.x + sample.y * sample.y - sample.x * sample.y);
        float length = sqrtf(sample.x * sample.x + sample.y * sample.y);
        float ratio = (length != 0) ? (sqRoot / length) : sqRoot;
        float x = regulationFactor * sample.x * ratio;
        float y = regulationFactor * sample.y * ratio;

        //shift the direction of the perfect reflection
        glm::vec3 glossyReflectionRayDirection = reflectedRay.direction + x * u + y * v;
        Ray glossyReflectionRay = Ray(reflectedRay.origin, glossyReflectionRayDirection, std::numeric_limits<float>::max());

        //sum up the conttributions of each ray
        sumOfInterference += renderRay(state, glossyReflectionRay, rayDepth + 1) * hitInfo.material.ks;
    }
    //normalize the sum of the rays and add the result to the current color
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