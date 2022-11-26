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
        BoundingBox();
        float x_left ; float x_right ; float y_up ; float y_down ; float z_near ; float z_far ;
        int object_id;
        float t; 
    };
    struct BoundingBoxTree
    {
        std::vector< BoundingBox > box_vec; 
        BoundingBoxTree *left; 
        BoundingBoxTree *right; 
        BoundingBoxTree *mid; 

        BoundingBox global_bounding_box; 
    };
    //check if bounding box got intersected 
    bool is_bounding_box_intersected( Ray & ray , BoundingBox & box );
    //generate every bounding box
    std::vector< BoundingBox > generate_bounding_boxes(parser::Scene &scene  );
    //setup every node 
    void setup_bvh_tree( BoundingBoxTree * &node , int splitting_axis    );
    // first thing to do 
    void setup_bvh_head_node( BoundingBoxTree * &head_node ,   std::vector< BoundingBox >&  global_bounding_box  );
    // traversing tree
    //returns bo containing object id 
    BoundingBox traverse_tree(BoundingBoxTree * head_node , Ray &ray , std::vector<std::pair<int, float> > &  object_no_t_pair  );
    // generate global bounding box
    void generate_global_bounding_box(BoundingBoxTree * &  head_node);


    // !!!!!!!! 
    // calll only this 
     void generate_BVH_tree(parser::Scene &scene  , BoundingBoxTree * & head_node );
    // delete bo 
    void delete_boxes( BoundingBoxTree * & node );
    //test
     void test_inorder_traversal(BoundingBoxTree * head_node);
}