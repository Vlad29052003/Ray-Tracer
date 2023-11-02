#pragma once

#include "fwd.h"
#include "bvh_interface.h"
#include "render.h"
#include "scene.h"
#include "screen.h"

// Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of depth of field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen);

// Generates a set of rays from the given pixel position directed towards the focus points.
std::vector<Ray> generateDofRaysForDebug(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, glm::ivec2 screenResolution, const glm::vec2& pixel, glm::vec3& focusPoint);

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// allowing objects to move during a render, and visualize the appearance of movement.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen);

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& screen);

// Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// For a description of the method's arguments, refer to 'extra.cpp'
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth);

// This method calculates the coordinates of the disk on which the rays are sampled.
// It returns an std::vector containing 4 glm::vectors, in this order:
//   the origin position
//   two vectors representing an orthonormal basis for the plane that contains the disk
//   a vector containing the radius as all 3 coordinates (just to make it simpler to return all 4)
std::vector<glm::vec3> glossyDebug(Ray ray, const HitInfo& hitInfo);

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You may have to add support for environment maps
// to the Scene object, and provide a scene with the right data for this.
// For a description of the method's arguments, refer to 'extra.cpp'
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray);

// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// For a description of the method's arguments, refer to 'bounding_volume_hierarchy.cpp'
// NOTE: this method is unit-tested, so do not change the function signature.
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVHInterface::Primitive> primitives);