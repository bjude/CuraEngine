//Copyright (c) 2017 Ultimaker B.V.
//CuraEngine is released under the terms of the AGPLv3 or higher.

#ifndef TREESUPPORT_H
#define TREESUPPORT_H

#include <forward_list>
#include <unordered_set>

#include <memory>
#include <vector>
#include <utility>
#include <unordered_map>

#include "settings/types/AngleRadians.h" //Creating the correct branch angles.

#include "sliceDataStorage.h"

namespace cura
{

/*!
 * \brief Generates a tree structure to support your models.
 */
class TreeSupport
{
public:
    /*!
     * \brief Creates an instance of the tree support generator.
     *
     * \param storage The data storage to get global settings from.
     */
    TreeSupport(const SliceDataStorage& storage);

    /*!
     * \brief Create the areas that need support.
     *
     * These areas are stored inside the given SliceDataStorage object.
     * \param storage The data storage where the mesh data is gotten from and
     * where the resulting support areas are stored.
     */
    void generateSupportAreas(SliceDataStorage& storage);

    /*!
     * \brief Represents the metadata of a node in the tree.
     */
    struct Node
    {
        static constexpr Node* NO_PARENT = nullptr;

        Node()
         : distance_to_top(0)
         , position(Point(0, 0))
         , skin_direction(false)
         , support_roof_layers_below(0)
         , to_buildplate(true)
         , parent(nullptr)
        {}

        Node(const Point position, const size_t distance_to_top, const bool skin_direction, const int support_roof_layers_below, const bool to_buildplate, Node* const parent)
         : distance_to_top(distance_to_top)
         , position(position)
         , skin_direction(skin_direction)
         , support_roof_layers_below(support_roof_layers_below)
         , to_buildplate(to_buildplate)
         , parent(parent)
        {}

#ifdef DEBUG // Clear the delete node's data so if there's invalid access after, we may get a clue by inspecting that node.
        ~Node()
        {
            parent = nullptr;
            merged_neighbours.clear();
        }
#endif // DEBUG

        /*!
         * \brief The number of layers to go to the top of this branch.
         */
        mutable size_t distance_to_top;

        /*!
         * \brief The position of this node on the layer.
         */
        Point position;

        /*!
         * \brief The direction of the skin lines above the tip of the branch.
         *
         * This determines in which direction we should reduce the width of the
         * branch.
         */
        mutable bool skin_direction;

        /*!
         * \brief The number of support roof layers below this one.
         *
         * When a contact point is created, it is determined whether the mesh
         * needs to be supported with support roof or not, since that is a
         * per-mesh setting. This is stored in this variable in order to track
         * how far we need to extend that support roof downwards.
         */
        mutable int support_roof_layers_below;

        /*!
         * \brief Whether to try to go towards the build plate.
         *
         * If the node is inside the collision areas, it has no choice but to go
         * towards the model. If it is not inside the collision areas, it must
         * go towards the build plate to prevent a scar on the surface.
         */
        mutable bool to_buildplate;

        /*!
         * \brief The originating node for this one, one layer higher.
         *
         * In order to prune branches that can't have any support (because they
         * can't be on the model and the path to the buildplate isn't clear),
         * the entire branch needs to be known.
         */
        Node *parent;

        /*!
        * \brief All neighbours (on the same layer) that where merged into this node.
        *
        * In order to prune branches that can't have any support (because they
        * can't be on the model and the path to the buildplate isn't clear),
        * the entire branch needs to be known.
        */
        mutable std::forward_list<Node*> merged_neighbours;

        bool operator==(const Node& other) const
        {
            return position == other.position;
        }
    };

private:
    /*!
     * \brief The border of the printer where we may not put tree branches.
     *
     * Lest they produce g-code that goes outside the build volume.
     */
    Polygons machine_volume_border;

    /*!
     * \brief Creates the areas that have to be avoided by the tree's branches.
     *
     * The result is a vector of 3D volumes that have to be avoided, where each
     * volume consists of a number of layers where the branch would collide with
     * the model.
     * There will be a volume for each sample of branch radius. The radii of the
     * branches are unknown at this point (there will be several radii at any
     * given layer too), so a collision area is generated for every possible
     * radius.
     *
     * \param storage The settings storage to get settings from.
     * \param model_collision[out] A vector to fill with the output collision
     * areas.
     */
    void collisionAreas(const SliceDataStorage& storage, std::vector<std::vector<Polygons>>& model_collision);

    /*!
     * \brief Draws circles around each node of the tree into the final support.
     *
     * This also handles the areas that have to become support roof, support
     * bottom, the Z distances, etc.
     *
     * \param storage[in, out] The settings storage to get settings from and to
     * save the resulting support polygons to.
     * \param contact_nodes The nodes to draw as support.
     * \param model_collision The model infill with the X/Y distance already
     * subtracted.
     */
    void drawCircles(SliceDataStorage& storage, const std::vector<std::unordered_set<Node*>>& contact_nodes, const std::vector<std::vector<Polygons>>& model_collision);

    /*!
     * \brief Drops down the nodes of the tree support towards the build plate.
     *
     * This is where the cleverness of tree support comes in: The nodes stay on
     * their 2D layers but on the next layer they are slightly shifted. This
     * causes them to move towards each other as they are copied to lower layers
     * which ultimately results in a 3D tree.
     *
     * \param contact_nodes[in, out] The nodes in the space that need to be
     * dropped down. The nodes are dropped to lower layers inside the same
     * vector of layers.
     * \param model_collision For each sample of radius, a list of layers with
     * the polygons of the collision areas of the model. Any node in there will
     * collide with the model.
     * \param model_avoidance For each sample of radius, a list of layers with
     * the polygons that must be avoided if the branches wish to go towards the
     * build plate.
     * \param model_internal_guide For each sample of radius, a list of layers
     * with the polygons that must be avoided if the branches wish to go towards
     * the model.
     */
    void dropNodes(std::vector<std::unordered_set<Node*>>& contact_nodes, const std::vector<std::vector<Polygons>>& model_collision, const std::vector<std::vector<Polygons>>& model_avoidance, const std::vector<std::vector<Polygons>>& model_internal_guide);

    /*!
     * \brief Creates points where support contacts the model.
     *
     * A set of points is created for each layer.
     * \param mesh The mesh to get the overhang areas to support of.
     * \param contact_nodes[out] A vector of mappings from contact points to
     * their tree nodes.
     * \param collision_areas For every layer, the areas where a generated
     * contact point would immediately collide with the model due to the X/Y
     * distance.
     * \return For each layer, a list of points where the tree should connect
     * with the model.
     */
    void generateContactPoints(const SliceMeshStorage& mesh, std::vector<std::unordered_set<Node*>>& contact_nodes, const std::vector<Polygons>& collision_areas);

    /*!
     * \brief Add a node to the next layer.
     *
     * If a node is already at that position in the layer, the nodes are merged.
     */
    void insertDroppedNode(std::unordered_set<Node*>& nodes_layer, Node* node);

    /*!
     * \brief Creates the areas that have to be avoided by the tree's branches
     * in order to reach the build plate.
     *
     * The result is a vector of 3D volumes that have to be avoided, where each
     * volume consists of a number of layers where the branch would collide with
     * the model.
     * There will be a volume for each sample of branch radius. The radii of the
     * branches are unknown at this point (there will be several radii at any
     * given layer too), so a collision area is generated for every possible
     * radius.
     *
     * The input collision areas are inset by the maximum move distance and
     * propagated upwards. This generates volumes so that the branches can
     * predict in time when they need to be moving away in order to avoid
     * hitting the model.
     * \param storage The settings storage to get settings from.
     * \param model_collision The collision areas that may not be hit by the
     * model.
     * \param model_avoidance[out] A vector to fill with the output avoidance
     * areas.
     */
    void propagateCollisionAreas(const SliceDataStorage& storage, const std::vector<std::vector<Polygons>>& model_collision, std::vector<std::vector<Polygons>>& model_avoidance);
};

}

namespace std
{
    template<> struct hash<cura::TreeSupport::Node>
    {
        size_t operator()(const cura::TreeSupport::Node& node) const
        {
            return hash<cura::Point>()(node.position);
        }
    };
}

namespace std {
template <typename T1, typename T2>
struct hash<std::pair<T1, T2>> {
    size_t operator()(const std::pair<T1, T2>& val) const {
        const auto first = std::hash<T1>{}(val.first);
        const auto second = std::hash<T2>{}(val.second);
        return first ^ (second + 0x9e3779b9 + (first << 6) + (first >> 2));
    }
};
}

namespace cura
{
namespace Tree
{
class Node;

struct TreeParams
{
    coord_t branch_radius;
    coord_t radius_sample_resolution;
    coord_t layer_height;
    coord_t xy_distance;
    coord_t max_move;
    coord_t radius_increment;
    coord_t point_spread;
    coord_t z_gap;
    size_t support_roof_layers;
    AngleRadians support_angle;
    coord_t initial_radius;
    bool can_support_on_model;
    BuildPlateShape buildplate_shape;
    EPlatformAdhesion adhesion_type;
    coord_t brim_size;
    coord_t raft_margin;
    coord_t skirt_size;
    coord_t line_width;
    int wall_count;
};

Polygons calculate_machine_border(const SliceDataStorage& storage, const TreeParams& params);
Point moveTowards(const Point& point, const Point& target, const Polygons& invalid, coord_t move_limit);
std::vector<Polygons> circlePolygons(const std::vector<std::unique_ptr<Node>>& nodes);

class ModelVolumes
{
public:
    ModelVolumes(const TreeParams& params, const SliceDataStorage& storage);

    ModelVolumes(const ModelVolumes&) = delete;
    ModelVolumes& operator=(const ModelVolumes&) = delete;

    const Polygons& collision(coord_t radius, int layer) const;
    const Polygons& avoidance(coord_t radius, int layer) const;
    const Polygons& internal_model(coord_t radius, int layer) const;

private:
    using RadiusLayerPair = std::pair<coord_t, int>;

    const TreeParams params_;
    Polygons machine_border_;
    std::vector<Polygons> layer_outlines_;
    mutable std::unordered_map<RadiusLayerPair, Polygons> collision_cache_;
    mutable std::unordered_map<RadiusLayerPair, Polygons> avoidance_cache_;
    mutable std::unordered_map<RadiusLayerPair, Polygons> internal_model_cache_;
};

class Node
{
public:
    Node() = default;
    Node(const Point& pos, coord_t radius, int layer, std::vector<std::unique_ptr<Node>> children = {},
         Node* parent = nullptr);

    Node& operator=(Node&&) = default;
    ~Node() = default;

    void merge(std::unique_ptr<Node> other);
    void merge(std::vector<std::unique_ptr<Node>> others);
    template <typename It>
    void merge(It first, It last)
    {
        while (first != last)
        {
            merge(std::move(*first));
            ++first;
        }
    }

    Node(const Node&) = delete;
    Node& operator=(const Node&) = delete;

    const Point& position() const { return position_; }
    const coord_t& radius() const { return radius_; }
    const int& layer() const { return layer_; }
    const std::vector<std::unique_ptr<Node>>& children() const { return children_; }
    const Node* parent() { return parent_; }

    void position(const Point& pos) { position_ = pos; }
    void parent(Node* parent) { parent_ = parent; }

private:
    Point position_{0, 0};
    coord_t radius_{0};
    int layer_{0};
    std::vector<std::unique_ptr<Node>> children_{};
    Node* parent_{nullptr};
};

class TreeSupport
{
    using NodePtr = std::unique_ptr<Node>;
    using NodePtrVec = std::vector<NodePtr>;

public:
    TreeSupport(const TreeParams& params, const SliceDataStorage& storage);
    void generateSupportAreas(SliceDataStorage& storage);

private:
    void processLayer();
    void dropNodes();
    void removeUnsupportableByBuildPlate();
    NodePtrVec generateContactPoints(const SliceDataStorage& data) const;
    NodePtrVec generateContactPoints(const SliceMeshStorage& mesh) const;
    std::vector<Point> generateContactSamplePoints(const SliceMeshStorage& mesh) const;
    void drawCircles(SliceDataStorage& storage) const;

    std::vector<Node*> gatherNodes(int layer) const;

    std::vector<NodePtrVec::iterator> groupNodes();
    int currentLayer() const;

    TreeParams params_;
    ModelVolumes volumes_;
    NodePtrVec trees_;
};
}
}

#endif /* TREESUPPORT_H */
