#include "../Header/BVH.h"
#include <iostream>
#include <algorithm>
#include <functional>
BVH::BoundingBox::BoundingBox( float x_left , float x_right , float y_down , float y_up , float z_near , float z_far)
{
    //init 
    this->x_left = x_left;
    this->x_right = x_right;
    this->y_down = y_down;
    this->y_up = y_up;
    this->z_far = z_far;
    this->z_near = z_near;

}
BVH::BoundingBox::BoundingBox()
{

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
    float tmin_y =  (box.y_up - ray.origin.y ) / ray.direction.y;
    float tmax_y =  (box.y_down - ray.origin.y ) / ray.direction.y;

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
    if( maximum_of_minimum > minimum_of_maximum )
    {
        return false; 
    }
    if( maximum_of_minimum < 0 && minimum_of_maximum < 0 )
    {
        return false;
    }
    
    box.t = maximum_of_minimum; //init t 
    if( maximum_of_minimum < 0 )
    {
        box.t = minimum_of_maximum; 
    }
    parser::Vec3f intersection_p;
    intersection_p.x = maximum_of_minimum * ray.direction.x + ray.origin.x;
    intersection_p.y = maximum_of_minimum * ray.direction.y + ray.origin.y;
    intersection_p.z = maximum_of_minimum * ray.direction.z + ray.origin.z;


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
            if(  !( scene.meshes[i].faces[j].v0_id-1 >= 0  && scene.meshes[i].faces[j].v1_id-1 >= 0  && scene.meshes[i].faces[j].v2_id-1 >= 0 ) )
            {
                continue; 
            }
            parser::Vec3f p1 = parser::Vec3f( scene.vertex_data[ (scene.meshes[i].faces[j].v0_id-1)    ].x , scene.vertex_data[ (scene.meshes[i].faces[j].v0_id-1)    ].y , scene.vertex_data[ (scene.meshes[i].faces[j].v0_id-1)  ].z); 
            parser::Vec3f p2 = parser::Vec3f( scene.vertex_data[ (scene.meshes[i].faces[j].v1_id-1)   ].x , scene.vertex_data[ (scene.meshes[i].faces[j].v1_id-1)  ].y , scene.vertex_data[ (scene.meshes[i].faces[j].v1_id-1)   ].z);
            parser::Vec3f p3 = parser::Vec3f( scene.vertex_data[ (scene.meshes[i].faces[j].v2_id-1)   ].x , scene.vertex_data[ (scene.meshes[i].faces[j].v2_id-1)    ].y , scene.vertex_data[ (scene.meshes[i].faces[j].v2_id-1)   ].z);

            //transformation 
            parser::Vec4f temp1; temp1.x = p1.x; temp1.y = p1.y; temp1.z = p1.z;temp1.w = 1; 
            parser::Vec4f temp2; temp2.x = p2.x; temp2.y = p2.y; temp2.z = p2.z;temp2.w = 1;
            parser::Vec4f temp3; temp3.x = p3.x; temp3.y = p3.y; temp3.z = p3.z;temp3.w = 1;

            temp1 = scene.meshes[i].transformation * temp1 ; 
            temp2 = scene.meshes[i].transformation * temp2 ; 
            temp3 = scene.meshes[i].transformation * temp3 ; 

            p1.x = temp1.x / temp1.w;
            p1.y = temp1.y / temp1.w;
            p1.z = temp1.z / temp1.w;

            p2.x = temp2.x / temp2.w;
            p2.y = temp2.y / temp2.w;
            p2.z = temp2.z / temp2.w;

            p3.x = temp3.x / temp3.w;
            p3.y = temp3.y / temp3.w;
            p3.z = temp3.z / temp3.w;
            
            

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

        parser::Vec3f center = scene.vertex_data[scene.spheres[i].center_vertex_id - 1 ];

        parser::Vec4f temp_radius;
        temp_radius.x = scene.spheres[i].radius;
        temp_radius.y = scene.spheres[i].radius; 
        temp_radius.z = scene.spheres[i].radius; 
        temp_radius.w = 1.0f; 


        //transformation 
        float maximum_radius_multiplier = scene.spheres[i].radius;

        parser::Vec4f temp; 
        temp.x= center.x;
        temp.y= center.y;
        temp.z= center.z;
        temp.w= 1;


        for (size_t j = 0; j < scene.spheres[i].transformations.size(); j++)
        {
            bool is_scale = false; 
            for (size_t k = 0; k < scene.spheres[i].scale_indices.size(); k++)
            {
                if(scene.spheres[i].scale_indices[k] == j )
                {
                    is_scale = true; 
                    break; 
                }
            }
            if( is_scale )
            {
                
                float max = 1 ;
                if( max < scene.spheres[i].transformation.elements[0])
                {
                    max = scene.spheres[i].transformation.elements[0];
                } 
                if( max < scene.spheres[i].transformation.elements[5])
                {
                    max = scene.spheres[i].transformation.elements[5];

                }
                if( max < scene.spheres[i].transformation.elements[10])
                {
                    max = scene.spheres[i].transformation.elements[10];
                }
                maximum_radius_multiplier *=  max;

            }
            else
            {
                temp =  scene.spheres[i].transformations[j] * temp ; 
            }
        }

        center.x = temp.x/ temp.w; 
        center.y = temp.y/ temp.w; 
        center.z = temp.z/ temp.w; 
        
        
        
        float x_min = center.x - maximum_radius_multiplier;
        float x_max = center.x + maximum_radius_multiplier;

        float y_min = center.y - maximum_radius_multiplier;
        float y_max = center.y + maximum_radius_multiplier;

        float z_min = center.z - maximum_radius_multiplier;
        float z_max = center.z + maximum_radius_multiplier;     

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
            parser::Vec4f temp1; temp1.x = p1.x; temp1.y = p1.y; temp1.z = p1.z;temp1.w = 1; 
            parser::Vec4f temp2; temp2.x = p2.x; temp2.y = p2.y; temp2.z = p2.z;temp2.w = 1;
            parser::Vec4f temp3; temp3.x = p3.x; temp3.y = p3.y; temp3.z = p3.z;temp3.w = 1;


            temp1 = scene.triangles[i].transformation * temp1;
            temp2 = scene.triangles[i].transformation * temp2;
            temp3 = scene.triangles[i].transformation * temp3;

            p1.x = temp1.x / temp1.w;
            p1.y = temp1.y / temp1.w;
            p1.z = temp1.z / temp1.w;

            p2.x = temp2.x / temp2.w;
            p2.y = temp2.y / temp2.w;
            p2.z = temp2.z / temp2.w;

            p3.x = temp3.x / temp3.w;
            p3.y = temp3.y / temp3.w;
            p3.z = temp3.z / temp3.w;
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
    //for each mesh instance
    for (size_t i = 0; i < scene.mesh_instances.size(); i++)
    {
        float x_min = INFINITY;
        float x_max = -INFINITY;

        float y_min = INFINITY;
        float y_max = -INFINITY;

        float z_min = INFINITY;
        float z_max= -INFINITY;

        
        for (size_t j = 0; j < scene.mesh_instances[i].mesh_ptr->faces.size(); j++)
        {
            //get points
            parser::Vec3f p1 = parser::Vec3f( scene.vertex_data[ (scene.mesh_instances[i].mesh_ptr->faces[j].v0_id-1)    ].x , scene.vertex_data[ (scene.mesh_instances[i].mesh_ptr->faces[j].v0_id-1)    ].y , scene.vertex_data[ (scene.mesh_instances[i].mesh_ptr->faces[j].v0_id-1)  ].z); 
            parser::Vec3f p2 = parser::Vec3f( scene.vertex_data[ (scene.mesh_instances[i].mesh_ptr->faces[j].v1_id-1)   ].x , scene.vertex_data[ (scene.mesh_instances[i].mesh_ptr->faces[j].v1_id-1)  ].y , scene.vertex_data[ (scene.mesh_instances[i].mesh_ptr->faces[j].v1_id-1)   ].z);
            parser::Vec3f p3 = parser::Vec3f( scene.vertex_data[ (scene.mesh_instances[i].mesh_ptr->faces[j].v2_id-1)   ].x , scene.vertex_data[ (scene.mesh_instances[i].mesh_ptr->faces[j].v2_id-1)    ].y , scene.vertex_data[ (scene.mesh_instances[i].mesh_ptr->faces[j].v2_id-1)   ].z);

            //transformation 
            parser::Vec4f temp1; temp1.x = p1.x; temp1.y = p1.y; temp1.z = p1.z;temp1.w = 1; 
            parser::Vec4f temp2; temp2.x = p2.x; temp2.y = p2.y; temp2.z = p2.z;temp2.w = 1;
            parser::Vec4f temp3; temp3.x = p3.x; temp3.y = p3.y; temp3.z = p3.z;temp3.w = 1;

            temp1 = scene.mesh_instances[i].transformation * temp1 ; 
            temp2 = scene.mesh_instances[i].transformation * temp2 ; 
            temp3 = scene.mesh_instances[i].transformation * temp3 ; 

            p1.x = temp1.x / temp1.w;
            p1.y = temp1.y / temp1.w;
            p1.z = temp1.z / temp1.w;

            p2.x = temp2.x / temp2.w;
            p2.y = temp2.y / temp2.w;
            p2.z = temp2.z / temp2.w;

            p3.x = temp3.x / temp3.w;
            p3.y = temp3.y / temp3.w;
            p3.z = temp3.z / temp3.w;
            
            

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
        int object_id = i + scene.meshes.size() + scene.spheres.size() + scene.triangles.size();
        //push the object
        box.object_id = object_id; 
        bounding_box_vec.push_back(box);        
    }
    return bounding_box_vec; 
}


void BVH::setup_bvh_tree( BVH::BoundingBoxTree * & node , int splitting_axis    )
{
   
    if( node->box_vec.size() ==  1 ) // if leaf node  
    {
        node->global_bounding_box = node->box_vec[0];
        node->left = NULL;
        node->right = NULL;
    return; 
    }
    else if( node->box_vec.size()  == 0 ) // if nothing 
    {
        node = NULL; 
        return;
    }
    //initialzie  left and right
    node->left = new BVH::BoundingBoxTree; 
    node->right = new BVH::BoundingBoxTree; 
    //sort all the objects in that axis
    if( splitting_axis == 0 )
    {


    
        float mid_point_x = (node->global_bounding_box.x_left + node->global_bounding_box.x_right) / 2;
          
        for (size_t i = 0; i < node->box_vec.size(); i++)
        {
            BVH::BoundingBox b = node->box_vec[i];
            //down
            if(  b.x_left <= mid_point_x && b.x_right <= mid_point_x )
            {
                node->left->box_vec.push_back(node->box_vec[i]);
            }
            // up
            else
            {
                node->right->box_vec.push_back(node->box_vec[i]);
            }
            
        }
        generate_global_bounding_box(node->left);
        generate_global_bounding_box(node->right);

        
        splitting_axis = 1; 
    }
    else if( splitting_axis == 1 )
    {

        float mid_point_y = (node->global_bounding_box.y_down + node->global_bounding_box.y_up) / 2;
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
        generate_global_bounding_box(node->left);
        generate_global_bounding_box(node->right);

        splitting_axis = 2; 

        
    }
    else if( splitting_axis == 2 ) // splitting axis = 2 
    {

    
        float mid_point_z = (node->global_bounding_box.z_far + node->global_bounding_box.z_near) / 2;

        for (size_t i = 0; i < node->box_vec.size(); i++)
        {
            BVH::BoundingBox b = node->box_vec[i];
           //down
            if(  b.z_near <= mid_point_z && b.z_far <= mid_point_z )
            {
                node->left->box_vec.push_back(node->box_vec[i]);
            }
            // up
            else
            {
                node->right->box_vec.push_back(node->box_vec[i]);
            }
            

        }
        generate_global_bounding_box(node->left);
        generate_global_bounding_box(node->right);
        splitting_axis = 0; 

    }
    std::cout << " 1 ";
    // handled leaf and empty node above
    setup_bvh_tree( node->left , splitting_axis  );
    setup_bvh_tree( node->right , splitting_axis);

}
/*void BVH::setup_bvh_tree( BVH::BoundingBoxTree * & node , int splitting_axis    )
{
    if( node->box_vec.size() ==  1 ) // if leaf node  
    {
        node->global_bounding_box = node->box_vec[0];
        node->global_bounding_box.object_id = node->box_vec[0].object_id;

        node->left = NULL;
        node->right = NULL;
        node->mid = NULL;

    return; 
    }
    else if( node->box_vec.size()  == 0 ) // if nothing 
    {
        node = NULL; 
        return;
    }
     //initialzie  left and right
    node->left = new BVH::BoundingBoxTree; 
    node->right = new BVH::BoundingBoxTree; 
    node->mid = new BVH::BoundingBoxTree; 
    
     //if divided at least two
    int sep_no[3];
    sep_no[0] = 0; 
    sep_no[1] = 0; 
    sep_no[2] = 0; 


    if(splitting_axis == 0 ) //x
    {
        float mid_point = (node->global_bounding_box.x_right + node->global_bounding_box.x_left) / 2;


        for (size_t i = 0; i < node->box_vec.size(); i++)
        {
            if( mid_point >= node->box_vec[i].x_left && mid_point >= node->box_vec[i].x_right )
            {
                node->left->box_vec.push_back( node->box_vec[i]);
                sep_no[0] = 1;
            }
            else if(mid_point <= node->box_vec[i].x_left && mid_point <= node->box_vec[i].x_right )
            {
                node->right->box_vec.push_back( node->box_vec[i]);
                sep_no[1] = 1;

            }
            else if(mid_point >= node->box_vec[i].x_left && mid_point <= node->box_vec[i].x_right )
            {
                node->mid->box_vec.push_back( node->box_vec[i]);
                sep_no[2] = 1;

            }
        }
         int total_sep_no = sep_no[0] + sep_no[1] +  sep_no[2];
        if( total_sep_no < 2 ) // in only 1 axis 
        {
            //go numeric
            for (float j = node->global_bounding_box.x_left; j < node->global_bounding_box.x_right; j+= 0.1f )
            {
                sep_no[0] = 0;
                sep_no[1] = 0;
                sep_no[2] = 0;

                mid_point = j;
                node->left->box_vec.clear();
                node->right->box_vec.clear();
                node->mid->box_vec.clear();
                for (size_t i = 0; i < node->box_vec.size(); i++)
                {
                    if( mid_point >= node->box_vec[i].x_left && mid_point >= node->box_vec[i].x_right )
                    {
                        node->left->box_vec.push_back( node->box_vec[i]);
                        sep_no[0] = 1;
                    }
                    else if(mid_point <= node->box_vec[i].x_left && mid_point <= node->box_vec[i].x_right )
                    {
                        node->right->box_vec.push_back( node->box_vec[i]);
                        sep_no[1] = 1;

                    }
                    else if(mid_point >= node->box_vec[i].x_left && mid_point <= node->box_vec[i].x_right )
                    {
                        node->mid->box_vec.push_back( node->box_vec[i]);
                        sep_no[2] = 1;

                    }
                }
                total_sep_no = sep_no[0] + sep_no[1] +  sep_no[2];
                if( total_sep_no > 1 )
                {
                    break; 
                }

            }
            //if no break yapacak birsey yok 
            
        }
        generate_global_bounding_box(node->left);
        generate_global_bounding_box(node->right);
        generate_global_bounding_box(node->mid);

        node->left->global_bounding_box.x_right = mid_point;
        node->right->global_bounding_box.x_left = mid_point;

        splitting_axis = 1; 
            
    }
    else if( splitting_axis == 1 )//y 
    {
        float mid_point = (node->global_bounding_box.y_up + node->global_bounding_box.y_down) / 2;
        // check if intersection
     
        for (size_t i = 0; i < node->box_vec.size(); i++)
        {
            if( mid_point >= node->box_vec[i].y_down && mid_point >= node->box_vec[i].y_up )
            {
                node->left->box_vec.push_back( node->box_vec[i]);
                sep_no[0] = 1;

            }
            else if(mid_point <= node->box_vec[i].y_down && mid_point <= node->box_vec[i].y_up)
            {
                node->right->box_vec.push_back( node->box_vec[i]);
                sep_no[1] = 1;

            }
            else if(mid_point >= node->box_vec[i].y_down && mid_point <= node->box_vec[i].y_up)
            {
                node->mid->box_vec.push_back( node->box_vec[i]);
                sep_no[2] = 1;

            }
        }
        int total_sep_no = sep_no[0] + sep_no[1] +  sep_no[2];
        if( total_sep_no < 2 ) // in only 1 axis 
        {
            //go numeric
            for (float j = node->global_bounding_box.y_down; j < node->global_bounding_box.y_up; j+= 0.1f )
            {
                sep_no[0] = 0;
                sep_no[1] = 0;
                sep_no[2] = 0;

                mid_point = j;
                node->left->box_vec.clear();
                node->right->box_vec.clear();
                node->mid->box_vec.clear();
                for (size_t i = 0; i < node->box_vec.size(); i++)
                {
                    if( mid_point >= node->box_vec[i].y_down && mid_point >= node->box_vec[i].y_up )
                    {
                        node->left->box_vec.push_back( node->box_vec[i]);
                        sep_no[0] = 1;

                    }
                    else if(mid_point <= node->box_vec[i].y_down && mid_point <= node->box_vec[i].y_up)
                    {
                        node->right->box_vec.push_back( node->box_vec[i]);
                        sep_no[1] = 1;

                    }
                    else if(mid_point >= node->box_vec[i].y_down && mid_point <= node->box_vec[i].y_up)
                    {
                        node->mid->box_vec.push_back( node->box_vec[i]);
                        sep_no[2] = 1;

                    }
                }
                total_sep_no = sep_no[0] + sep_no[1] +  sep_no[2];
                if( total_sep_no > 1 )
                {
                    break; 
                }

            }
            //if no break yapacak birsey yok 
            
        }
        generate_global_bounding_box(node->left);
        generate_global_bounding_box(node->right);
        generate_global_bounding_box(node->mid);

        node->left->global_bounding_box.y_up = mid_point;
        node->right->global_bounding_box.y_down = mid_point;
        
        
        splitting_axis = 2; 

    }
    else if( splitting_axis == 2 ) //z 
    {
        float mid_point = (node->global_bounding_box.z_far + node->global_bounding_box.z_near) / 2;
        
        for (size_t i = 0; i < node->box_vec.size(); i++)
        {
            if( mid_point >= node->box_vec[i].z_near && mid_point >= node->box_vec[i].z_far )
            {
                node->left->box_vec.push_back( node->box_vec[i]);
                sep_no[0] = 1;

            }
            else if(mid_point <= node->box_vec[i].z_near && mid_point <= node->box_vec[i].z_far   )
            {
                node->right->box_vec.push_back( node->box_vec[i]);
                sep_no[1] = 1;

            }
            else if(mid_point >= node->box_vec[i].z_near && mid_point <= node->box_vec[i].z_far   )
            {
                node->mid->box_vec.push_back( node->box_vec[i]);
                sep_no[2] = 1;

            }
        }
        
        int total_sep_no = sep_no[0] + sep_no[1] +  sep_no[2];
        if( total_sep_no < 2 ) // in only 1 axis 
        {
            //go numeric
            for (float j = node->global_bounding_box.z_near; j < node->global_bounding_box.z_far; j+= 0.1f )
            {
                sep_no[0] = 0;
                sep_no[1] = 0;
                sep_no[2] = 0;

                mid_point = j;
                node->left->box_vec.clear();
                node->right->box_vec.clear();
                node->mid->box_vec.clear();
                for (size_t i = 0; i < node->box_vec.size(); i++)
                {
                    if( mid_point >= node->box_vec[i].z_near && mid_point >= node->box_vec[i].z_far )
                    {
                        node->left->box_vec.push_back( node->box_vec[i]);
                        sep_no[0] = 1;

                    }
                    else if(mid_point <= node->box_vec[i].z_near && mid_point <= node->box_vec[i].z_far   )
                    {
                        node->right->box_vec.push_back( node->box_vec[i]);
                        sep_no[1] = 1;

                    }
                    else if(mid_point >= node->box_vec[i].z_near && mid_point <= node->box_vec[i].z_far   )
                    {
                        node->mid->box_vec.push_back( node->box_vec[i]);
                        sep_no[3] = 1;

                    }
                }
                total_sep_no = sep_no[0] + sep_no[1] +  sep_no[2];
                if( total_sep_no > 1 )
                {
                    break; 
                }
               

            }
            //if no break yapacak birsey yok 
            
        }
        
        generate_global_bounding_box(node->left);
        generate_global_bounding_box(node->right);
        generate_global_bounding_box(node->mid);

        node->left->global_bounding_box.z_far = mid_point;
        node->right->global_bounding_box.z_near = mid_point;
        splitting_axis = 0; 

    }

    
    // handled leaf and empty node above
    setup_bvh_tree( node->left , splitting_axis  );
    setup_bvh_tree( node->right , splitting_axis);
    setup_bvh_tree( node->mid , splitting_axis);

}*/
void BVH::setup_bvh_head_node( BVH::BoundingBoxTree * & head_node ,   std::vector< BVH::BoundingBox > & global_bounding_box  )
{
    // all of the elements
    std::cout << global_bounding_box.size() << std::endl; 

    head_node = new BoundingBoxTree; 
    head_node->box_vec = global_bounding_box; 
    generate_global_bounding_box(head_node);
    //get the  global out of it
    

}

//return an object id
//if no intersection return -1 
BVH::BoundingBox BVH::traverse_tree(BVH::BoundingBoxTree * head_node , Ray &ray , std::vector<std::pair<int, float> > &  object_no_t_pair  )
{
    BVH::BoundingBox box(INFINITY , INFINITY , INFINITY ,INFINITY  , INFINITY ,INFINITY);
    box.object_id = -1;
    box.t = INFINITY;

    if( !is_bounding_box_intersected(ray , head_node->global_bounding_box ))
    {
        return box; 
    }
    // box is global box
    box = head_node->global_bounding_box;
   
    
    // get obj id
    if(head_node->box_vec.size() == 1 )
    {
        box.object_id = head_node->box_vec[0].object_id;

        object_no_t_pair.push_back(std::make_pair( box.object_id , box.t));

        return box; 
    }


    BVH::BoundingBox left_box(INFINITY , INFINITY , INFINITY ,INFINITY  , INFINITY ,INFINITY);
    left_box.object_id = -1;
    left_box.t = INFINITY;

    BVH::BoundingBox right_box(INFINITY , INFINITY , INFINITY ,INFINITY  , INFINITY ,INFINITY);
    right_box.object_id = -1;
    right_box.t = INFINITY;

    BVH::BoundingBox mid_box(INFINITY , INFINITY , INFINITY ,INFINITY  , INFINITY ,INFINITY);
    mid_box.object_id = -1;
    mid_box.t = INFINITY;

    if( head_node->left == NULL && head_node->right == NULL && head_node->mid == NULL  )
    {
        return box; 
    }
    
    if( head_node->left != NULL  )
    {
        left_box = traverse_tree(head_node->left , ray , object_no_t_pair );
    }
    if( head_node->right != NULL  )
    {
        right_box = traverse_tree(head_node->right , ray , object_no_t_pair );
    }
    if( head_node->mid != NULL )
    {
        mid_box = traverse_tree(head_node->mid , ray , object_no_t_pair );
    }
  
    if( right_box.t <= mid_box.t && right_box.t <=  left_box.t )
    {
        return  right_box;
    }
    else if( left_box.t <= mid_box.t && left_box.t <=  right_box.t )
    {
        return  left_box;
    }
    else if( mid_box.t <= left_box.t && mid_box.t <= right_box.t )
    {
        return mid_box; 
    }
    else
    {
        BVH::BoundingBox box(INFINITY , INFINITY , INFINITY ,INFINITY  , INFINITY ,INFINITY);
        box.object_id = -1;
        box.t = INFINITY;
        return box; 
    }
    

}

 void BVH::generate_BVH_tree(parser::Scene &scene  , BVH::BoundingBoxTree * &head_node )
{
    std::vector<BVH::BoundingBox> bounding_boxes = BVH::generate_bounding_boxes(scene );
    std::cout << " end of genenrating boundng boxes " << std::endl; 
    BVH::setup_bvh_head_node(head_node , bounding_boxes);
    std::cout << " end of setup head node " << std::endl; 

    BVH::setup_bvh_tree(head_node , 0 );
    std::cout << " end of setup bvh tree  " << std::endl; 

}
void BVH::generate_global_bounding_box(BVH::BoundingBoxTree * &  node)
{
    //step1 compute global bounding box of all objects
    float min_x = INFINITY;
    float max_x = -INFINITY; 
    float min_y = INFINITY;
    float max_y = -INFINITY;
    float min_z = INFINITY;
    float max_z = -INFINITY;
    //init node
    for (size_t i = 0; i < node->box_vec.size(); i++)
    {
        BVH::BoundingBox b = node->box_vec[i];
        if( min_x > b.x_left )
        {
            min_x = b.x_left; 
        }
        if( max_x < b.x_right )
        {
            max_x = b.x_right; 
        }
        if( min_y > b.y_down )
        {
            min_y = b.y_down; 
        }
        if( max_y < b.y_up )
        {
            max_y = b.y_up; 
        }
        if( min_z > b.z_near )
        {
            min_z = b.z_near; 
        }
        if( max_z < b.z_far )
        {
            max_z = b.z_far; 
        }
    }
    BVH::BoundingBox global(min_x , max_x , min_y , max_y , min_z , max_z );
    node->global_bounding_box = global;
}

 void BVH::test_inorder_traversal(BVH::BoundingBoxTree * head_node)
{
    if( head_node != NULL )
    {
       std::cout << head_node->global_bounding_box.x_left << " " << head_node->global_bounding_box.x_right << " " << head_node->global_bounding_box.y_down << " " << head_node->global_bounding_box.y_up << " " << head_node->global_bounding_box.z_near << " " <<head_node->global_bounding_box.z_far << " size == "<<  head_node->box_vec.size() <<  std::endl ;
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

void BVH::delete_boxes( BVH::BoundingBoxTree * & node )
{
    if( node->left == NULL  && node->right == NULL && node->mid == NULL )
    {
        delete node;
    }
    if( node->left != NULL )
    {
        delete_boxes(node->left );
    }
    if( node->right != NULL )
    {
        delete_boxes(node->right );
    }
    if( node->mid != NULL )
    {
        delete_boxes(node->mid );
    }
}
