#pragma once
//#define DEBUG
#include "../Header/Prototypes.h"
#include "../Dependencies/glm_headers/glm.hpp"
#include "../Header/Ray.h"
// the basic ray tracing block 
static void generate_image(parser::Scene & scene)
{
    for (size_t camera_no = 0; camera_no < scene.cameras.size() ; camera_no++)
    {
        current_camera = scene.cameras[camera_no];  // get the current_camera

        int width = current_camera.image_width, height = current_camera.image_height;
        unsigned char* image = new unsigned char [width * height * 3];

        //calculate the starting point of  (0 , 0 ) in world space
        parser::Vec3f starting_point; // up left  (0 , 0 )
        parser::Vec3f interval_row; // real space interval in two consec pixels in same row 
        parser::Vec3f interval_col;// real space interval in two consec pixels in same  col

        calculate_image_plane(current_camera , starting_point , interval_row , interval_col );
        //we now have a startin point
        glm::vec3 starting_point_glm = glm::vec3( starting_point.x , starting_point.y , starting_point.z );  
        for (size_t y = 0; y < height; y++)
        {
            for (size_t x = 0; x < width; x++)
            {
                // initialize the ray
                glm::vec3 current_pixel_world_space = starting_point_glm + glm::vec3( interval_row.x ,  interval_row.y ,  interval_row.z) * (float) x + (float ) y * glm::vec3( interval_col.x ,  interval_col.y ,  interval_col.z) ; 
                Ray ray(glm::vec3( current_camera.position.x , current_camera.position.y , current_camera.position.z  )  , current_pixel_world_space );
                float color = color_pixel(scene  , ray);
            }
        }
        //save the written image
        write_ppm( current_camera.image_name.c_str(), image, width, height);
        //in the end just reallocate the image
        delete[] image; 
    }
}

static void calculate_image_plane(const parser::Camera &current_camera , parser::Vec3f & starting_point,  parser::Vec3f & interval_row, parser::Vec3f & interval_col)
{
    //calculate right vector
    glm::vec3 cam_up = glm::normalize(glm::vec3(current_camera.up.x,current_camera.up.y,current_camera.up.z));
    glm::vec3 cam_gaze = glm::normalize(glm::vec3(current_camera.gaze.x,current_camera.gaze.y,current_camera.gaze.z));

    glm::vec3 right_vec = glm::normalize(glm::cross(cam_gaze , cam_up ));
    std::cout << right_vec.x << " " << right_vec.y <<  "  " << right_vec.z << std::endl; 
    //find the intersection with the image
    glm::vec3 intersection =  glm::vec3(current_camera.position.x , current_camera.position.y , current_camera.position.z   )  +  current_camera.near_distance *  cam_gaze;

    //now get to the (0 ,0 ) 
    
    // but we cannot use right we need left
    glm::vec3 left_vec = glm::vec3( -1 * right_vec.x , -1 * right_vec.y , -1 * right_vec.z  );
    #ifdef DEBUG 
    std::cout << " left vec" << left_vec.x << " " << left_vec.y << " " << left_vec.z << std::endl;
    std::cout << " right vec" << right_vec.x << " " << right_vec.y << " " << right_vec.z << std::endl;
    #endif
    // check how much we need to go left
    float left = current_camera.near_plane.x;
    float right = current_camera.near_plane.y;
    float bottom = current_camera.near_plane.z;
    float top = current_camera.near_plane.w;
    
    glm::vec3 top_left_corner = intersection + left_vec * left + cam_up * top;
    starting_point.x = top_left_corner.x;
    starting_point.y = top_left_corner.y;
    starting_point.z = top_left_corner.z;

    //lastly calculate the intervals

    glm::vec3 most_left = intersection + right_vec * left; 
    glm::vec3 most_right = intersection + right_vec * right;

    

    glm::vec3 most_top = intersection + top * cam_up; 
    glm::vec3 most_bottom = intersection + bottom * cam_up;

    #ifdef DEBUG
    std::cout << "most right " <<  most_right.x  << " " << most_right.y <<" " << most_right.z << std::endl;
    std::cout << "most left " <<  most_left.x  << " " << most_left.y <<" " << most_left.z << std::endl;
    std::cout << "most top " <<  most_top.x  << " " << most_top.y <<" " << most_top.z << std::endl;
    std::cout << "most bottom " <<  most_bottom.x  << " " << most_bottom.y <<" " << most_bottom.z << std::endl;
    #endif DEBUG
    glm::vec3 interval_row_glm  = (most_right - most_left  );
    interval_row_glm = interval_row_glm / (float)current_camera.image_width;
    glm::vec3 interval_col_glm = ( most_bottom - most_top );
    interval_col_glm = interval_col_glm / (float)current_camera.image_height;

    #ifdef DEBUG
    std::cout << "interval row " << interval_row_glm.x <<  " " << interval_row_glm.y << " " << interval_row_glm.z << std::endl;
    std::cout << "interval col " << interval_col_glm.x <<  " " << interval_col_glm.y << " " << interval_col_glm.z << std::endl;
    #endif DEBUG

    interval_row.x = interval_row_glm.x;
    interval_row.y = interval_row_glm.y;
    interval_row.z = interval_row_glm.z;

    interval_col.x = interval_col_glm.x;
    interval_col.y = interval_col_glm.y;
    interval_col.z = interval_col_glm.z;




    //glm::vec3 interval_col_glm  = (most_ - most_right) / current_camera.image_width;  







}

static float color_pixel(parser::Scene& scene , Ray & ray )
{
    // 1 - check intersections with messhes if no intersection return 0
    bool is_intersected = false;
    for (size_t i = 0; i < scene.meshes.size(); i++) // travere each object
    {
        glm::vec3 intersection_point(0.0f , 0.0f ,0.0f);  
        bool is_intersected_with_this_object = false; 
        is_intersected_with_this_object = calculate_intersection(scene , scene.meshes[i] , intersection_point);
    }
    // 2 - check intersections with spheres
    for (size_t i = 0; i < scene.spheres.size(); i++) // travere each object
    {
        glm::vec3 intersection_point(0.0f , 0.0f ,0.0f);  
        bool is_intersected_with_this_sphere = false; 
        is_intersected_with_this_sphere = calculate_intersection(scene , scene.spheres[i] , intersection_point);
    }    
    return 1.0f;
}


// intersections and stuff 
static bool ray_triangle_intersection( const Ray &ray , const parser::Scene& scene ,  glm::vec3 intersection_point  )
{
     
}
static bool ray_plane_intersection(const Ray& ray , const Plane& plane , glm::vec3 & hitpoint  )
{

}
static bool ray_sphere_intersection(const Ray& ray , const Sphere& sphere , glm::vec3 & hitpoint  )
{
    // at^2 + bT + c = 0
    glm::vec3 O = ray.origin; 
    glm::vec3 k = ray.direction; 

    float a  = 3.0f;
    float b = 2 * ( O.x * k.x + O.y * k.y + O.z * k.z ) -2*(k.x * sphere.center.x + k.y * sphere.center.y + k.z * sphere.center.z ); 
    float c = glm::dot(O , O ) + glm::dot(sphere.center,sphere.center) + -2*(O.x*sphere.center.x + O.y*sphere.center.y + O.z*sphere.center.z);
    float discriminant =  b*b - 4 * a * c; 
    if( discriminant < 0 )
    {
        return false; 
    }
    if( discriminant == 0 )
    {
        float root = -b-glm::sqrt(discriminant) / (2*a);

        hitpoint = O + root * k;
        return true;  
    }
    // else there are two roots

    float root_1 = -b-glm::sqrt(discriminant) / (2*a);
    float root_2 = -b-glm::sqrt(discriminant) / (2*a);

    //two hitpoints
    glm::vec3 h1 = O + root_1 * k;
    glm::vec3 h2 = O + root_2 * k;
    //get the smaller length root * k;
    float d1 = glm::distance( O ,  h1 );
    float d2 = glm::distance( O ,  h2 );

    if( d1 > d2 )
    {
        hitpoint = h2; 
    }
    else
    {
        hitpoint = h1;
    }
    return true; 

}


static bool calculate_intersection(parser::Scene& scene ,parser::Mesh& object ,  glm::vec3 & intersection_point )
{
    return true; 
}
static bool calculate_intersection(parser::Scene& scene ,parser::Sphere& object ,  glm::vec3 & intersection_point )
{
    return true; 
}