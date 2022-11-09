#include "../Header/BVH.h"
#include <iostream>
BVH::BoundingBox::BoundingBox( float x_left , float x_right , float y_up , float y_down , float z_near , float z_far)
{
    //init 
    this->x_left = x_left;
    this->x_right = x_right;
    this->y_down = y_down;
    this->y_up = y_up;
    this->z_far = z_far;
    this->z_near = z_near;

}

bool BVH::is_bounding_box_intersected( Ray & ray , BVH::BoundingBox & box )
{
    //init t as inf
    box.t = INFINITY;
    //for x 
    float tmin_x =  (box.x_left - ray.origin.x ) / ray.direction.x;
    float tmax_x =  (box.x_right - ray.origin.x ) / ray.direction.x;

    if( tmin_x > tmax_x)
    {
        float temp = tmax_x;
        tmax_x = tmin_x;
        tmin_x = temp; 
    }

    //for y 
    float tmin_y =  (box.y_down - ray.origin.y ) / ray.direction.y;
    float tmax_y =  (box.y_up - ray.origin.y ) / ray.direction.y;

    if( tmin_y > tmax_y)
    {
        float temp = tmax_y;
        tmax_y = tmin_y;
        tmin_y = temp; 
    }


   
    //for z 
    float tmin_z =  (box.z_far - ray.origin.z ) / ray.direction.z;
    float tmax_z =  (box.z_near - ray.origin.z ) / ray.direction.z;
    if( tmin_z > tmax_z)
    {
        float temp = tmax_z;
        tmax_z = tmin_z;
        tmin_z = temp; 
     }

    
    // find maximu of minimums
    float maximum_of_minimum; 
    if( tmin_x > tmin_y && tmin_x > tmin_z  )
    {
        maximum_of_minimum = tmin_x;
    }
    else if(tmin_z > tmin_x && tmin_z > tmin_y  )
    {
        maximum_of_minimum = tmin_z;
    }
    else
    {
        maximum_of_minimum = tmin_y;
    }

     // find minimum  of maxmimum
    float minimum_of_maximum; 
    if( tmax_x < tmax_y && tmax_x < tmax_z  )
    {
        minimum_of_maximum = tmax_x;
    }
    else if(tmax_z <  tmax_x && tmax_z < tmax_y  )
    {
        minimum_of_maximum = tmax_z;
    }
    else
    {
        minimum_of_maximum = tmax_y;
    }
    std::cout << minimum_of_maximum <<  "   " << maximum_of_minimum << std::endl; 
    if( maximum_of_minimum > minimum_of_maximum )
    {
        return false; 
    }
    box.t = maximum_of_minimum; //init t 
    return true; 
}

std::vector< BVH::BoundingBox > BVH::generate_bounding_boxes(parser::Scene &scene  )
{
    std::vector< BVH::BoundingBox > bounding_box_vec; 
    std::cout << "scene mesh size " << scene.meshes.size() << std::endl;
    //for each mesh 
    for (size_t i = 0; i < scene.meshes.size(); i++)
    {
        float x_min = INFINITY;
        float x_max = -INFINITY;

        float y_min = INFINITY;
        float y_max = -INFINITY;

        float z_min = INFINITY;
        float z_max= -INFINITY;

        for (size_t j = 0; j < scene.meshes[i].faces.size(); j++)
        {
            //get points
            parser::Vec3f p1 = parser::Vec3f( scene.vertex_data[ (scene.meshes[i].faces[j].v0_id-1)    ].x , scene.vertex_data[ (scene.meshes[i].faces[j].v0_id-1)    ].y , scene.vertex_data[ (scene.meshes[i].faces[j].v0_id-1)  ].z); 
            parser::Vec3f p2 = parser::Vec3f( scene.vertex_data[ (scene.meshes[i].faces[j].v1_id-1)   ].x , scene.vertex_data[ (scene.meshes[i].faces[j].v1_id-1)  ].y , scene.vertex_data[ (scene.meshes[i].faces[j].v1_id-1)   ].z);
            parser::Vec3f p3 = parser::Vec3f( scene.vertex_data[ (scene.meshes[i].faces[j].v2_id-1)   ].x , scene.vertex_data[ (scene.meshes[i].faces[j].v2_id-1)    ].y , scene.vertex_data[ (scene.meshes[i].faces[j].v2_id-1)   ].z);
            // if transformation change the points
            if( scene.meshes[i].transformations.size() > 0 )
            {
                parser::Vec4f temp1; temp1.x = p1.x; temp1.y = p1.y; temp1.z = p1.z;temp1.w = 1; 
                parser::Vec4f temp2; temp2.x = p2.x; temp2.y = p2.y; temp2.z = p2.z;temp2.w = 1;
                parser::Vec4f temp3; temp3.x = p3.x; temp3.y = p3.y; temp3.z = p3.z;temp3.w = 1;
                for (size_t k = 0; i < scene.meshes[i].transformations.size(); i++)
                {
                    temp1 = scene.meshes[i].transformations[k] * temp1;
                    temp2 = scene.meshes[i].transformations[k] * temp2;
                    temp3 = scene.meshes[i].transformations[k] * temp3;
                }
                p1.x = temp1.x / temp1.w;
                p1.y = temp1.y / temp1.w;
                p1.z = temp1.z / temp1.w;

                p2.x = temp2.x / temp2.w;
                p2.y = temp2.y / temp2.w;
                p2.z = temp2.z / temp2.w;

                p3.x = temp3.x / temp3.w;
                p3.y = temp3.y / temp3.w;
                p3.z = temp3.z / temp3.w;
            
            }
            
            // update box
            //p1        
            if( p1.x >  x_max )
            {
                x_max = p1.x; 
            }
            if( p1.x <  x_min )
            {
                x_min = p1.x; 
            }
            if( p1.y >  y_max )
            {
                y_max = p1.y; 
            }
            if( p1.y <  y_min )
            {
                y_min = p1.y; 
            }
            if( p1.z >  z_max )
            {
                z_max = p1.z; 
            }
            if( p1.z <  z_min )
            {
                z_min = p1.z; 
            }
            // p2 
             if( p2.x >  x_max )
            {
                x_max = p2.x; 
            }
            if( p2.x <  x_min )
            {
                x_min = p2.x; 
            }
            if( p2.y >  y_max )
            {
                y_max = p2.y; 
            }
            if( p2.y <  y_min )
            {
                y_min = p2.y; 
            }
            if( p2.z >  z_max )
            {
                z_max = p2.z; 
            }
            if( p2.z <  z_min )
            {
                z_min = p2.z; 
            }
            //p3 
             if( p3.x >  x_max )
            {
                x_max = p3.x; 
            }
            if( p3.x <  x_min )
            {
                x_min = p3.x; 
            }
            if( p3.y >  y_max )
            {
                y_max = p3.y; 
            }
            if( p3.y <  y_min )
            {
                y_min = p3.y; 
            }
            if( p3.z >  z_max )
            {
                z_max = p3.z; 
            }
            if( p3.z <  z_min )
            {
                z_min = p3.z; 
            }
        }
        //create a box
        BVH::BoundingBox box(x_min , x_max , y_min , y_max , z_min , z_max);
        int object_id = i;
        //push the object
        box.object_id = object_id; 
        bounding_box_vec.push_back(box);        
    }
    
    // sphere 
    for (size_t i = 0; i < scene.spheres.size(); i++)
    {
        parser::Vec3f center = scene.vertex_data[scene.spheres[i].center_vertex_id];
        //transformation 
        if( scene.spheres[i].transformations.size() >  0 )
        {
            parser::Vec4f temp; 
            temp.x= center.x;
            temp.y= center.y;
            temp.z= center.z;
            temp.w= 1;
            parser::Sphere temp_sphere = (parser::Sphere)scene.spheres[i]; 
            parser::Matrix inverse_scale(4,4);
            for (size_t i = 0; i < temp_sphere.transformations.size(); i++)
            {
                temp =   temp_sphere.transformations[i] * temp ; 
            }
            center.x = temp.x/ temp.w; 
            center.y = temp.y/ temp.w; 
            center.z = temp.z/ temp.w; 
        }
        
        
        float x_min = center.x - scene.spheres[i].radius;
        float x_max = center.x + scene.spheres[i].radius;

        float y_min = center.y - scene.spheres[i].radius;
        float y_max = center.y + scene.spheres[i].radius;

        float z_min =center.z - scene.spheres[i].radius;
        float z_max= center.z + scene.spheres[i].radius;     

         //create a box
        BVH::BoundingBox box(x_min , x_max , y_min , y_max , z_min , z_max);
        int object_id = i + scene.meshes.size();
        //push the object

        box.object_id = object_id; 
        bounding_box_vec.push_back(box);
    }

    //for each triangle
    for (size_t i = 0; i < scene.triangles.size(); i++)
    {
        //get points
            parser::Vec3f p1 = parser::Vec3f( scene.vertex_data[ (scene.triangles[i].indices.v0_id-1)    ].x , scene.vertex_data[(scene.triangles[i].indices.v0_id-1)    ].y , scene.vertex_data[ (scene.triangles[i].indices.v0_id-1)  ].z); 
            parser::Vec3f p2 = parser::Vec3f( scene.vertex_data[ (scene.triangles[i].indices.v1_id-1)   ].x , scene.vertex_data[ (scene.triangles[i].indices.v1_id-1)  ].y , scene.vertex_data[ (scene.triangles[i].indices.v1_id-1)   ].z);
            parser::Vec3f p3 = parser::Vec3f( scene.vertex_data[ (scene.triangles[i].indices.v2_id-1)   ].x , scene.vertex_data[ (scene.triangles[i].indices.v2_id-1)    ].y , scene.vertex_data[ (scene.triangles[i].indices.v2_id-1)   ].z);
            // update box
            //transformations
            if( scene.triangles[i].transformations.size() > 0 )
            {
                parser::Vec4f temp1; temp1.x = p1.x; temp1.y = p1.y; temp1.z = p1.z;temp1.w = 1; 
                parser::Vec4f temp2; temp2.x = p2.x; temp2.y = p2.y; temp2.z = p2.z;temp2.w = 1;
                parser::Vec4f temp3; temp3.x = p3.x; temp3.y = p3.y; temp3.z = p3.z;temp3.w = 1;
                for (size_t j = 0; j < scene.triangles[i].transformations.size(); i++)
                {
                    temp1 = scene.triangles[i].transformations[j] * temp1;
                    temp2 = scene.triangles[i].transformations[j] * temp2;
                    temp3 = scene.triangles[i].transformations[j] * temp3;

                }
                p1.x = temp1.x / temp1.w;
                p1.y = temp1.y / temp1.w;
                p1.z = temp1.z / temp1.w;

                p2.x = temp2.x / temp2.w;
                p2.y = temp2.y / temp2.w;
                p2.z = temp2.z / temp2.w;

                p3.x = temp3.x / temp3.w;
                p3.y = temp3.y / temp3.w;
                p3.z = temp3.z / temp3.w;
            }
            //p1        
            float x_min = INFINITY;
            float x_max = -INFINITY;

            float y_min = INFINITY;
            float y_max = -INFINITY;

            float z_min = INFINITY;
            float z_max= -INFINITY;
            if( p1.x >  x_max )
            {
                x_max = p1.x; 
            }
            if( p1.x <  x_min )
            {
                x_min = p1.x; 
            }
            if( p1.y >  y_max )
            {
                y_max = p1.y; 
            }
            if( p1.y <  y_min )
            {
                y_min = p1.y; 
            }
            if( p1.z >  z_max )
            {
                z_max = p1.z; 
            }
            if( p1.z <  z_min )
            {
                z_min = p1.z; 
            }
            // p2 
             if( p2.x >  x_max )
            {
                x_max = p2.x; 
            }
            if( p2.x <  x_min )
            {
                x_min = p2.x; 
            }
            if( p2.y >  y_max )
            {
                y_max = p2.y; 
            }
            if( p2.y <  y_min )
            {
                y_min = p2.y; 
            }
            if( p2.z >  z_max )
            {
                z_max = p2.z; 
            }
            if( p2.z <  z_min )
            {
                z_min = p2.z; 
            }
            //p3 
             if( p3.x >  x_max )
            {
                x_max = p3.x; 
            }
            if( p3.x <  x_min )
            {
                x_min = p3.x; 
            }
            if( p3.y >  y_max )
            {
                y_max = p3.y; 
            }
            if( p3.y <  y_min )
            {
                y_min = p3.y; 
            }
            if( p3.z >  z_max )
            {
                z_max = p3.z; 
            }
            if( p3.z <  z_min )
            {
                z_min = p3.z; 
            }    

         //create a box
        BVH::BoundingBox box(x_min , x_max , y_min , y_max , z_min , z_max);
        int object_id = i + scene.meshes.size() + scene.spheres.size();

        box.object_id = object_id; 
        //push the object
       bounding_box_vec.push_back(box);
    }
    return bounding_box_vec; 
}


void BVH::setup_bvh_tree( BVH::BoundingBoxTree * node   )
{
   
    //step1 compute global bounding box of all objects
    float min_x = INFINITY;
    float max_x = -INFINITY; 
    float min_y = INFINITY;
    float max_y = -INFINITY;
    float min_z = INFINITY;
    float max_z = -INFINITY;
    for (size_t i = 0; i < node->box_vec.size(); i++)
    {
        BVH::BoundingBox b = node->box_vec[i];
        if( min_x > b.x_left )
        {
            min_x = b.x_left; 
        }
        if( max_x > b.x_right )
        {
            max_x = b.x_right; 
        }
        if( min_y > b.y_down )
        {
            min_y = b.y_down; 
        }
        if( max_y > b.y_up )
        {
            max_y = b.y_up; 
        }
        if( min_z > b.z_near )
        {
            min_z = b.z_near; 
        }
        if( max_z > b.z_far )
        {
            max_z = b.z_far; 
        }
    }
    BVH::BoundingBox global(min_x , max_x , min_y , max_y , min_z , max_z );
    node->global_bounding_box = global; 
    if( node->box_vec.size() ==  1 ) // if leaf node  
    {
        node->left = NULL;
        node->right = NULL;
        return; 
    }
    else if( node->box_vec.size()  == 0 ) // if nothing 
    {
        node = NULL; 
        return;
    }

    // step 2 decide on a splitting axis
    int splitting_axis = -1; // 0 - x 1 - y 2 - z 
    float x_dif = std::abs(global.x_right  - global.x_left);
    float y_dif = std::abs(global.y_up  - global.y_down);
    float z_dif = std::abs(global.z_far  - global.z_near);
    if( x_dif  > y_dif &&  x_dif  > z_dif )
    {
        splitting_axis = 0; 
    }
    else if(  y_dif  > x_dif &&  y_dif  > z_dif  )
    {
        splitting_axis = 1; 
    }
    else
    {
        splitting_axis = 2; 
    }


    // simple method choose middle point
    if( splitting_axis == 0 )
    {
        float mid_point_x = (global.x_right  + global.x_left) / 2 ;
        for (size_t i = 0; i < node->box_vec.size(); i++)
        {
            BVH::BoundingBox b = node->box_vec[i];
            //left
            if(  b.x_left <= mid_point_x && b.x_right <= mid_point_x )
            {
                node->left->box_vec.push_back(node->box_vec[i]);
            }
            // right
            else
            {
                node->right->box_vec.push_back(node->box_vec[i]);
            }
        }
        
    }
    else if( splitting_axis == 1 )
    {
        float mid_point_y = (global.y_up  + global.y_down) / 2 ;
        for (size_t i = 0; i < node->box_vec.size(); i++)
        {
            BVH::BoundingBox b = node->box_vec[i];
            //down
            if(  b.y_down <= mid_point_y && b.y_up <= mid_point_y )
            {
                node->left->box_vec.push_back(node->box_vec[i]);
            }
            // up
            else
            {
                node->right->box_vec.push_back(node->box_vec[i]);
            }
        }
    }
    else // splitting axis = 2 
    {
        float mid_point_z = (global.z_far  + global.z_near) / 2 ;
         for (size_t i = 0; i < node->box_vec.size(); i++)
        {
            BVH::BoundingBox b = node->box_vec[i];
            //left
            if(  b.z_near <= mid_point_z && b.z_far <= mid_point_z )
            {
                node->left->box_vec.push_back(node->box_vec[i]);
            }
            // right
            else
            {
                node->left->box_vec.push_back(node->box_vec[i]);
            }
        }
    }
    // handled leaf and empty node above
    setup_bvh_tree( node->left );
    setup_bvh_tree( node->right );

}
void BVH::setup_bvh_head_node( BVH::BoundingBoxTree * head_node ,   std::vector< BVH::BoundingBox > & global_bounding_box  )
{
    // all of the elements
    head_node->box_vec = global_bounding_box; 
}

//return an object id
//if no intersection return -1 
BVH::BoundingBox BVH::traverse_tree(BVH::BoundingBoxTree * head_node , Ray &ray )
{
    BVH::BoundingBox box(INFINITY , INFINITY , INFINITY ,INFINITY  , INFINITY ,INFINITY);
    box.object_id = -1;
    box.t = INFINITY;

    if( !is_bounding_box_intersected(ray , head_node->global_bounding_box ))
    {
        return box; 
    }
    BVH::BoundingBox left_box(INFINITY , INFINITY , INFINITY ,INFINITY  , INFINITY ,INFINITY);
    left_box.object_id = -1;
    left_box.t = INFINITY;

    BVH::BoundingBox right_box(INFINITY , INFINITY , INFINITY ,INFINITY  , INFINITY ,INFINITY);
    right_box.object_id = -1;
    right_box.t = INFINITY;
    if( head_node->left != NULL  )
    {
        left_box = traverse_tree(head_node->left , ray );
    }
    if( head_node->right != NULL  )
    {
        right_box = traverse_tree(head_node->right , ray );
    }
    
    if( right_box.t == INFINITY && left_box.t == INFINITY)
    {
        return  right_box;
    }
    if( left_box.t == INFINITY)
    {
        return right_box;
    }
    if ( right_box.t == INFINITY )
    {
        return left_box;
    }
    // now the hard case
    //return smallest t
    if( right_box.t < left_box.t )
    {
        return right_box;
    }
    else
    {
        return left_box; 
    }
    

}

static void BVH::generate_BVH_tree(parser::Scene &scene  , BVH::BoundingBoxTree * head_node )
{
    std::vector<BVH::BoundingBox> bounding_boxes = BVH::generate_bounding_boxes(scene );
    BVH::setup_bvh_head_node(head_node , bounding_boxes);
    BVH::setup_bvh_tree(head_node);

}
static void BVH::test_inorder_traversal(BoundingBoxTree * head_node)
{
    if( head_node != NULL )
    {
    std::cout << head_node->global_bounding_box.x_left << " " << head_node->global_bounding_box.x_right << " " << head_node->global_bounding_box.y_down << " " << head_node->global_bounding_box.y_up << " " << head_node->global_bounding_box.z_near << " " << head_node->global_bounding_box.z_far << " " << std::endl;
    }
    if( head_node->left != NULL )
    {
        BVH::test_inorder_traversal(head_node->left);
    }
    if( head_node->right != NULL )
    {
        BVH::test_inorder_traversal(head_node->right);
    }
}