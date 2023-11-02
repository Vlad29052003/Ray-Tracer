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