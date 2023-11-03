#include "intersect.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <limits>


bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p) {
    glm::vec3 v0_to_p = p - v0;
    glm::vec3 v0_to_v1 = v1 - v0;
    glm::vec3 v0_to_v2 = v2 - v0;

    float d00 = glm::dot(v0_to_v1, v0_to_v1);
    float d01 = glm::dot(v0_to_v1, v0_to_v2);
    float d11 = glm::dot(v0_to_v2, v0_to_v2);
    float d20 = glm::dot(v0_to_p, v0_to_v1);
    float d21 = glm::dot(v0_to_p, v0_to_v2);

    float denom = d00 * d11 - d01 * d01;

    float beta = (d11 * d20 - d01 * d21) / denom;
    float gamma = (d00 * d21 - d01 * d20) / denom;
    float alpha = 1.0f - beta - gamma;

    if (alpha >= 0 && beta >= 0 && gamma >= 0 && alpha < 1 && beta < 1 && gamma < 1)
    {
        return true;
    }

    return false;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    float dot = glm::dot(plane.normal, ray.direction);
    if (std::abs(dot) < 1e-6)
        return false;
    float t = ((-glm::dot(ray.origin, plane.normal) + plane.D) / glm::dot(plane.normal, ray.direction));
    if(t <= 0) {
        return false;
    }
    if(t <= ray.t) {
        ray.t = t;
        return true;
    }
    return false;
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;
    glm::vec3 edge1 = v1 - v0;
    glm::vec3 edge2 = v2 - v0;
    plane.normal = glm::normalize(glm::cross(edge1, edge2));
    plane.D = glm::dot(plane.normal, v0);
    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{
    Plane plane = trianglePlane(v0, v1, v2);
    // Store previous value
    float t = ray.t;

    bool intersects = intersectRayWithPlane(plane, ray);

    if(!intersects) {
        // If ray doesn't intersect with plane, then rollback t
        ray.t = t;
        return false;
    }
    // Ray intersects with plane
    glm::vec3 point = ray.origin + ray.t * ray.direction;
    if(pointInTriangle(v0, v1, v2, plane.normal, point)) {
        return true;
    }
    // If ray doesn't intersect with triangle, then rollback t
    ray.t = t;
    return false;
}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    glm::vec3 oc = ray.origin - sphere.center;
    float A = glm::dot(ray.direction, ray.direction);
    float B = 2.0f * glm::dot(oc, ray.direction);
    float C = glm::dot(oc, oc) - sphere.radius * sphere.radius;

    float discriminant = B * B - 4 * A * C;
    if (discriminant < 0) {
        return false;
    } else {
        float t1 = (-B - std::sqrt(discriminant)) / (2.0f * A);
        float t2 = (-B + std::sqrt(discriminant)) / (2.0f * A);
        if (t1 > 0 && t1 < ray.t) {
            ray.t = t1;
            return true;
        }
        if (t2 > 0 && t2 < ray.t) {
            ray.t = t2;
            return true;
        }
        return false;
    }
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    float tmin = -std::numeric_limits<float>::infinity();
    float tmax = std::numeric_limits<float>::infinity();

    glm::vec3 bounds[2] = {box.lower, box.upper};

    for (int i = 0; i < 3; ++i) {
        float t1 = (bounds[0][i] - ray.origin[i]) / ray.direction[i];
        float t2 = (bounds[1][i] - ray.origin[i]) / ray.direction[i];
        if (t1 > t2) std::swap(t1, t2);
        tmin = std::max(tmin, t1);
        tmax = std::min(tmax, t2);
        if (tmin > tmax) return false;
    }
    if (ray.t < 0) {
        return false;
    }
    if (tmin >= 0 && tmin < ray.t) {
        ray.t = tmin;
        return true;
    }
    if (tmax < ray.t) {
        ray.t = tmax;
        return true;
    }

    return false;
}
