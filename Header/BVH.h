#pragma once 
#include "../System_Files/support_files/parser.h"
#include "../Header/Prototypes.h"
#include "../Header/Ray.h"

#include <vector>
namespace BVH
{
    class BoundingBox
    {
        public:
        BoundingBox( float x_left , float x_right , float y_down , float y_up , float z_near , float z_far );
    
        float x_left ; float x_right ; float y_up ; float y_down ; float z_near ; float z_far ;
        int object_id;
        float t; 
    };
    struct BoundingBoxTree
    {
        std::vector< BoundingBox > box_vec; 
        BoundingBoxTree *left; 
        BoundingBoxTree *right; 
        BoundingBox global_bounding_box; 

    };
    //check if bounding box got intersected 
    bool is_bounding_box_intersected( Ray & ray , BoundingBox & box );
    //generate every bounding box
    std::vector< BoundingBox > generate_bounding_boxes(parser::Scene &scene  );
    //setup every node 
    void setup_bvh_tree( BoundingBoxTree * node   );
    // first thing to do 
    void setup_bvh_head_node( BoundingBoxTree * head_node ,   std::vector< BoundingBox >&  global_bounding_box  );
    // traversing tree
    //returns bo containing object id 
    BoundingBox traverse_tree(BoundingBoxTree * head_node , Ray &ray );

    // calll only this 
    static void generate_BVH_tree(parser::Scene &scene  , BoundingBoxTree * head_node );

    //test
    static void test_inorder_traversal(BoundingBoxTree * head_node);
}