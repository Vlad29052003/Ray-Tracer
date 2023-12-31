#include "bvh.h"
#include "draw.h"
#include "extra.h"
#include "interpolate.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "texture.h"
#include <algorithm>
#include <bit>
#include <chrono>
#include <framework/opengl_includes.h>
#include <iostream>
#include <queue>
#include <queue>

// Helper method to fill in hitInfo object. This can be safely ignored (or extended).
// Note: many of the functions in this helper tie in to standard/extra features you will have
// to implement separately, see interpolate.h/.cpp for these parts of the project
void updateHitInfo(RenderState& state, const BVHInterface::Primitive& primitive, const Ray& ray, HitInfo& hitInfo)
{
    const auto& [v0, v1, v2] = std::tie(primitive.v0, primitive.v1, primitive.v2);
    const auto& mesh = state.scene.meshes[primitive.meshID];
    const auto n = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
    const auto p = ray.origin + ray.t * ray.direction;

    // First, fill in default data, unrelated to separate features
    hitInfo.material = mesh.material;
    hitInfo.normal = n;
    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, p);

    // Next, if `features.enableNormalMapping` is true, generate smoothly interpolated vertex normals
    if (state.features.enableNormalInterp) {
        hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
    }

    // Next, if `features.enableTextureMapping` is true, generate smoothly interpolated vertex uvs
    if (state.features.enableTextureMapping) {
        hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
    }

    // Finally, catch flipped normals
    if (glm::dot(ray.direction, n) > 0.0f) {
        hitInfo.normal = -hitInfo.normal;
    }
}

// BVH constructor; can be safely ignored. You should not have to touch this
// NOTE: this constructor is tested, so do not change the function signature.
BVH::BVH(const Scene& scene, const Features& features)
{
#ifndef NDEBUG
    // Store start of bvh build for timing
    using clock = std::chrono::high_resolution_clock;
    const auto start = clock::now();
#endif

    // Count the total nr. of triangles in the scene
    size_t numTriangles = 0;
    for (const auto& mesh : scene.meshes)
        numTriangles += mesh.triangles.size();

    // Given the input scene, gather all triangles over which to build the BVH as a list of Primitives
    std::vector<Primitive> primitives;
    primitives.reserve(numTriangles);
    for (uint32_t meshID = 0; meshID < scene.meshes.size(); meshID++) {
        const auto& mesh = scene.meshes[meshID];
        for (const auto& triangle : mesh.triangles) {
            primitives.push_back(Primitive {
                .meshID = meshID,
                .v0 = mesh.vertices[triangle.x],
                .v1 = mesh.vertices[triangle.y],
                .v2 = mesh.vertices[triangle.z] });
        }
    }

    // Tell underlying vectors how large they should approximately be
    m_primitives.reserve(numTriangles);
    m_nodes.reserve(numTriangles + 1);

    // Recursively build BVH structure; this is where your implementation comes in
    m_nodes.emplace_back(); // Create root node
    m_nodes.emplace_back(); // Create dummy node s.t. children are allocated on the same cache line
    buildRecursive(scene, features, primitives, RootIndex);

    // Fill in boilerplate data
    buildNumLevels();
    buildNumLeaves();

#ifndef NDEBUG
    // Output end of bvh build for timing
    const auto end = clock::now();
    std::cout << "BVH construction time: " << std::chrono::duration<double, std::milli>(end - start).count() << "ms" << std::endl;
#endif
}

// BVH helper method; allocates a new node and returns its index
// You should not have to touch this
uint32_t BVH::nextNodeIdx()
{
    const auto idx = static_cast<uint32_t>(m_nodes.size());
    m_nodes.emplace_back();
    return idx;
}

// TODO: Standard feature
// Given a BVH triangle, compute an axis-aligned bounding box around the primitive
// - primitive; a single triangle to be stored in the BVH
// - return;    an axis-aligned bounding box around the triangle
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computePrimitiveAABB(const BVHInterface::Primitive primitive)
{
    glm::vec3 min = glm::min(primitive.v0.position, glm::min(primitive.v1.position, primitive.v2.position));
    glm::vec3 max = glm::max(primitive.v0.position, glm::max(primitive.v1.position, primitive.v2.position));

    return { .lower = min, .upper = max };
}

// TODO: Standard feature
// Given a range of BVH triangles, compute an axis-aligned bounding box around the range.
// - primitive; a contiguous range of triangles to be stored in the BVH
// - return;    a single axis-aligned bounding box around the entire set of triangles
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computeSpanAABB(std::span<const BVHInterface::Primitive> primitives)
{
    glm::vec3 min(FLT_MAX, FLT_MAX, FLT_MAX);
    glm::vec3 max(-FLT_MAX, -FLT_MAX, -FLT_MAX);

    for (const auto &bvhTriangle : primitives){
        AxisAlignedBox box = computePrimitiveAABB(bvhTriangle);
        min.x = std::min(min.x, box.lower.x);
        min.y = std::min(min.y, box.lower.y);
        min.z = std::min(min.z, box.lower.z);

        max.x = std::max(max.x, box.upper.x);
        max.y = std::max(max.y, box.upper.y);
        max.z = std::max(max.z, box.upper.z);
    }

    return { .lower = min, .upper = max };
}

// TODO: Standard feature
// Given a BVH triangle, compute the geometric centroid of the triangle
// - primitive; a single triangle to be stored in the BVH
// - return;    the geometric centroid of the triangle's vertices
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePrimitiveCentroid(const BVHInterface::Primitive primitive)
{
    return (primitive.v0.position + primitive.v1.position + primitive.v2.position) / 3.0f;
}

// TODO: Standard feature
// Given an axis-aligned bounding box, compute the longest axis; x = 0, y = 1, z = 2.
// - aabb;   the input axis-aligned bounding box
// - return; 0 for the x-axis, 1 for the y-axis, 2 for the z-axis
//           if several axes are equal in length, simply return the first of these
// This method is unit-tested, so do not change the function signature.
uint32_t computeAABBLongestAxis(const AxisAlignedBox& aabb)
{
    glm::vec3 dimensions = aabb.upper - aabb.lower;

    if (dimensions.x >= dimensions.y && dimensions.x >= dimensions.z) {
        return 0;
    } else if (dimensions.y >= dimensions.x && dimensions.y >= dimensions.z) {
        return 1;
    } else {
        return 2;
    }
}

// TODO: Standard feature
// Given a range of BVH triangles, sort these along a specified axis based on their geometric centroid.
// Then, find and return the split index in the range, such that the subrange containing the first element
// of the list is at least as big as the other, and both differ at most by one element in size.
// Hint: you should probably reuse `computePrimitiveCentroid()`
// - aabb;       the axis-aligned bounding box around the given triangle range
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires sorting/splitting along an axis
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesByMedian(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVHInterface::Primitive> primitives)
{
    using Primitive = BVHInterface::Primitive;

    std::sort(primitives.begin(), primitives.end(),
        [axis](const Primitive& a, const Primitive& b) {
            glm::vec3 centroidA = computePrimitiveCentroid(a);
            glm::vec3 centroidB = computePrimitiveCentroid(b);
            return centroidA[axis] < centroidB[axis];
        });
    size_t medianIndex = primitives.size() / 2;

    return primitives.size() % 2 == 1 ? medianIndex + 1 : medianIndex;
}


// TODO: Standard feature
// Hierarchy traversal routine; called by the BVH's intersect(),
// you must implement this method and implement it carefully!
//
// If `features.enableAccelStructure` is not enabled, the method should just iterate the BVH's
// underlying primitives (or the scene's geometry). The default imlpementation already does this.
// You will have to implement the part which actually traverses the BVH for a faster intersect,
// given that `features.enableAccelStructure` is enabled.
//
// This method returns `true` if geometry was hit, and `false` otherwise. On first/closest hit, the
// distance `t` in the `ray` object is updated, and information is updated in the `hitInfo` object.
//
// - state;    the active scene, and a user-specified feature config object, encapsulated
// - bvh;      the actual bvh which should be traversed for faster intersection
// - ray;      the ray intersecting the scene's geometry
// - hitInfo;  the return object, with info regarding the hit geometry
// - return;   boolean, if geometry was hit or not
//
// This method is unit-tested, so do not change the function signature.
bool intersectRayWithBVH(RenderState& state, const BVHInterface& bvh, Ray& ray, HitInfo& hitInfo)
{

    using Node = BVHInterface::Node;
    using Primitive = BVHInterface::Primitive;

    // Relevant data in the constructed BVH
    std::span<const Node> nodes = bvh.nodes();
    std::span<const Primitive> primitives = bvh.primitives();

    // Return value
    bool is_hit = false;

    if (state.features.enableAccelStructure) {
        float prev = ray.t;
        std::queue<uint32_t> s;
        s.push(0);

        while (!s.empty()) {
            Node node = nodes[s.front()];
            s.pop();

            if (node.isLeaf()) {
                for (int i = node.primitiveOffset(); i < node.primitiveOffset() + node.primitiveCount(); ++i) {
                    if (intersectRayWithTriangle(primitives[i].v0.position, primitives[i].v1.position, primitives[i].v2.position, ray, hitInfo)) {
                        updateHitInfo(state, primitives[i], ray, hitInfo);
                    }
                }
            } else {
                float previous = ray.t;
                    if(intersectRayWithShape(nodes[node.rightChild()].aabb, ray)) {
                        ray.t = previous;
                        s.push(node.rightChild());
                    }
                    if(intersectRayWithShape(nodes[node.leftChild()].aabb, ray)) {
                        ray.t = previous;
                        s.push(node.leftChild());
                    }
            }
        }
        if(prev > ray.t)
            is_hit = true;
    } else {
        // Naive implementation; simply iterates over all primitives
        for (const auto& prim : primitives) {
            const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                updateHitInfo(state, prim, ray, hitInfo);
                is_hit = true;
            }
        }
    }

    // Intersect with spheres.
    for (const auto& sphere : state.scene.spheres)
        is_hit |= intersectRayWithShape(sphere, ray, hitInfo);
    return is_hit;
}

// TODO: Standard feature
// Leaf construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and a range of triangles, generate a valid leaf object
// and store the triangles in the `m_primitives` vector.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;      the active scene
// - features;   the user-specified features object
// - aabb;       the axis-aligned bounding box around the primitives beneath this leaf
// - primitives; the range of triangles to be stored for this leaf
BVH::Node BVH::buildLeafData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, std::span<Primitive> primitives)
{
    Node node;
    node.aabb = aabb;
    uint32_t primitiveOffset = static_cast<uint32_t>(m_primitives.size());
    uint32_t primitiveCount = static_cast<uint32_t>(primitives.size());
    node.data[0] = primitiveOffset | Node::LeafBit;
    node.data[1] = primitiveCount;
    std::copy(primitives.begin(), primitives.end(), std::back_inserter(m_primitives));
    return node;
}

// TODO: Standard feature
// Node construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and left/right child indices, generate a valid node object.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;           the active scene
// - features;        the user-specified features object
// - aabb;            the axis-aligned bounding box around the primitives beneath this node
// - leftChildIndex;  the index of the node's left child in `m_nodes`
// - rightChildIndex; the index of the node's right child in `m_nodes`
BVH::Node BVH::buildNodeData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, uint32_t leftChildIndex, uint32_t rightChildIndex)
{
    Node node;
    node.aabb = aabb;
    node.data[0] = leftChildIndex;
    node.data[1] = rightChildIndex;
    return node;
}

// TODO: Standard feature
// Hierarchy construction routine; called by the BVH's constructor,
// you must implement this method and implement it carefully!
//
// You should implement the other BVH standard features first, and this feature last, as you can reuse
// most of the other methods to assemble this part. There are detailed instructions inside the
// method which we recommend you follow.
//
// Arguments:
// - scene;      the active scene
// - features;   the user-specified features object
// - primitives; a range of triangles to be stored in the BVH
// - nodeIndex;  index of the node you are currently working on, this is already allocated
//
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildRecursive(const Scene& scene, const Features& features, std::span<Primitive> primitives, uint32_t nodeIndex)
{

    AxisAlignedBox aabb = computeSpanAABB(primitives);

    if(primitives.size() <= LeafSize) {
        m_nodes[nodeIndex] = buildLeafData(scene, features, aabb, primitives);
    }
    else {
        size_t idx;
        if(features.extra.enableBvhSahBinning) {
            idx = splitPrimitivesBySAHBin(aabb, computeAABBLongestAxis(aabb), primitives);
        } else {
            idx = splitPrimitivesByMedian(aabb, computeAABBLongestAxis(aabb), primitives);
        }
        uint32_t left = nextNodeIdx();
        buildRecursive(scene, features, primitives.subspan(0, idx), left);
        uint32_t right = nextNodeIdx();
        buildRecursive(scene, features, primitives.subspan(idx), right);
        m_nodes[nodeIndex] = buildNodeData(scene, features, aabb, left, right);
    }
}

int maxDepth(BVHInterface::Node node, std::vector<BVHInterface::Node> nodes)
{
    if (node.isLeaf()){
        return 1;
    }
    else {
        int lDepth = maxDepth(nodes[node.leftChild()], nodes);
        int rDepth = maxDepth(nodes[node.rightChild()], nodes);
        if (lDepth > rDepth)
            return (lDepth + 1);
        else
            return (rDepth + 1);
    }
}

void BVH::buildNumLevels()
{
    m_numLevels = maxDepth(m_nodes[0], m_nodes);
}


// Compute the nr. of leaves in your hierarchy after construction; useful for `debugDrawLeaf()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLeaves() {
    uint32_t cnt = 0;
    for(auto node: m_nodes) {
        if(node.isLeaf()) cnt++;
    }
    m_numLeaves = cnt;
}


// Draw the bounding boxes of the nodes at the selected level. Use this function to visualize nodes
// for debugging. You may wish to implement `buildNumLevels()` first. We suggest drawing the AABB
// of all nodes on the selected level.
// You are free to modify this function's signature.
void BVH::debugDrawLevel(int l)
{
    // Example showing how to draw an AABB as a (white) wireframe box.
    // Hint: use draw functions (see `draw.h`) to draw the contained boxes with different
    // colors, transparencies, etc.
    std::vector<Node> levelNodes;
    std::queue<std::pair<Node, int>> q;
    q.push({m_nodes[0], 0});

    while (!q.empty()) {
        Node current_node = q.front().first;
        int current_level = q.front().second;
        q.pop();

        if (current_level == l) {
            drawAABB(current_node.aabb, DrawMode::Wireframe, glm::vec3(0.5f, 1.0f, 0.5f), 1.0f);
        }

        if (current_level > l) {
            break;
        }

        if (!current_node.isLeaf()) {
            if (current_node.leftChild() != -1) {
                q.push({m_nodes[current_node.leftChild()], current_level + 1});
            }
            if (current_node.rightChild() != -1) {
                q.push({m_nodes[current_node.rightChild()], current_level + 1});
            }
        }
    }
}


// Draw data of the leaf at the selected index. Use this function to visualize leaf nodes
// for debugging. You may wish to implement `buildNumLeaves()` first. We suggest drawing the AABB
// of the selected leaf, and then its underlying primitives with different colors.
// - leafIndex; index of the selected leaf.
//              (Hint: not the index of the i-th node, but of the i-th leaf!)
// You are free to modify this function's signature.
void BVH::debugDrawLeaf(int leafIndex)
{
    int i = 0;
    while(leafIndex > 0 && i < m_nodes.size()) {
        if(m_nodes[i].isLeaf()) leafIndex--;
        i++;
    }
    drawAABB(m_nodes[i-1].aabb, DrawMode::Wireframe, glm::vec3(1.0f, 0, 0), 0.9f);
}