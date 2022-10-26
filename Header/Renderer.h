#pragma once
#include "../Header/Prototypes.h"
#include "../Dependencies/glm_headers/glm.hpp"
#include "../Header/Ray.h"
#include "../Header/Helper.h"
// the basic ray tracing block 
static void generate_image(parser::Scene & scene)
{
    int i = 0; 
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
        std::cout << "starting point " <<  starting_point_glm.x << " " << starting_point_glm.y  << starting_point_glm.z << std::endl; 
        for (size_t y = 0; y < height; y++)
        {
            for (size_t x = 0; x < width; x++)
            {
                //std::cout << " pixel no " << i /3 << std::endl; 
                // initialize the ray
                glm::vec3 current_pixel_world_space = starting_point_glm + glm::vec3( interval_row.x ,  interval_row.y ,  interval_row.z) * (float) x + (float ) y * glm::vec3( interval_col.x ,  interval_col.y ,  interval_col.z) ; 
                Ray ray(glm::vec3( current_camera.position.x , current_camera.position.y , current_camera.position.z  )  , current_pixel_world_space );
                glm::vec3  color = color_pixel(scene  , ray);
                
                // color cast 
                image[i++] = (unsigned char) (color.x);
                image[i++] = (unsigned char) (color.y);
                image[i++] = (unsigned char) (color.z);
                //std::cout << (unsigned int)color.x << " " <<  (unsigned int)color.y << " " <<  (unsigned int)color.z << std::endl; 
                //std::cout << "image " << (unsigned int )image[i-3 ] << " " << (unsigned int )image[i -2 ]  << " " << (unsigned int )image[i-1] << std::endl; 
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
    glm::vec3 left_vec = glm::normalize(glm::cross(cam_up , cam_gaze));
    #ifdef DEBUG 
    std::cout << " left vec" << left_vec.x << " " << left_vec.y << " " << left_vec.z << std::endl;
    std::cout << " right vec" << right_vec.x << " " << right_vec.y << " " << right_vec.z << std::endl;
    std::cout << " intersection " << intersection.x << " " << intersection.y << " " << intersection.z << std::endl; 
    std::cout << " cam_up " << cam_up.x << " " << cam_up.y << " " << cam_up.z << std::endl; 
    #endif
    // check how much we need to go left
    float left = current_camera.near_plane.x;
    float right = current_camera.near_plane.y;
    float bottom = current_camera.near_plane.z;
    float top = current_camera.near_plane.w;
    
    glm::vec3 top_left_corner = intersection + right_vec * left + cam_up * top;
    

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

static glm::vec3 color_pixel(parser::Scene& scene , Ray & ray )
{

    
    glm::vec3 color(0.0f , 0.0f , 0.0f );
    
    glm::vec3 hit_point; 
    glm::vec3 normal;
    parser::Material material; 
    
    //find the hitpoint 
    bool is_shadow_rays_active = false;
    glm::vec3 no_meaning_hit_point; 
    int object_id = -1; 
    bool  is_object_hit = ray_object_intersection( ray , scene ,  hit_point , normal  , material   , object_id ,  is_shadow_rays_active);
    
    if( !is_object_hit)
    {
        return glm::vec3(scene.background_color.x, scene.background_color.y , scene.background_color.z );
    }
    //normalize normal 
    //now we got the  nearest hitpoint and normal of that hitpoint. we can calculate color and cast shadow rays
    
    //we can add ambient shading
    color = glm::vec3(scene.ambient_light.x , scene.ambient_light.y , scene.ambient_light.z);
    //  for each light
    is_shadow_rays_active = true; 

    float least_cosine  = 9999; 
    for (size_t i = 0; i < scene.point_lights.size(); i++) 
    {
        glm::vec3 light_pos = glm::vec3(scene.point_lights[i].position.x , scene.point_lights[i].position.y ,scene.point_lights[i].position.z);
        
        //add shadow ray epsilon in direction of normal 
        hit_point += normal * scene.shadow_ray_epsilon;
        // cast  shadow ray
        Ray shadow_ray(hit_point , light_pos );
        
        glm::vec3 hit_point_temp; 
        glm::vec3 normal_temp;
        parser::Material material_temp;  
        
        glm::vec3 diffuse_1 = glm::vec3( material.diffuse.x , material.diffuse.y, material.diffuse.z);
        bool  is_object_hit = ray_object_intersection( shadow_ray , scene ,  hit_point_temp , normal_temp , material_temp , object_id , is_shadow_rays_active);
        
        if( is_object_hit ) // it is in shadow no contribution from light
        {

            continue; 
        }
        // else you directly went to a light
        if(material.is_mirror ) // do some weird stuff 
        {

        }
        
        // diffuse 
        glm::vec3 diffuse = glm::vec3( material.diffuse.x , material.diffuse.y, material.diffuse.z);
        float cosine_alpha = glm::dot(shadow_ray.direction  ,  normal ) / (glm::length(shadow_ray.direction) * glm::length(normal)  ) ; //normal vectors
        //clamp cosine alpha
        cosine_alpha = glm::max(0.0f , cosine_alpha);
        diffuse = diffuse * cosine_alpha; 
        diffuse = glm::vec3( diffuse.x * scene.point_lights[i].intensity.x , diffuse.y *  scene.point_lights[i].intensity.y ,diffuse.z *   scene.point_lights[i].intensity.z) / (glm::distance(light_pos , hit_point )*glm::distance(light_pos, hit_point)  );
        //clamp diffuse
        diffuse.x = glm::min(255.0f , diffuse.x );
        diffuse.y = glm::min(255.0f , diffuse.y );
        diffuse.z = glm::min(255.0f , diffuse.z );


        // specular 
        glm::vec3 specular = glm::vec3( material.specular.x , material.specular.y, material.specular.z);
        float phong_exponent = material.phong_exponent; 

        glm::vec3 h =  ( shadow_ray.direction + -1.0f*(ray.direction) )  / glm::length(shadow_ray.direction + -1.0f*(ray.direction));  // this might be flawed 
        std::cout << " h " << glm::length(h) << std::endl; 
        print_vec3(h);
        //std::cout << " normal " << std::endl; 
       // print_vec3(normal);
        float cosine_h_n = glm::dot(h , normal) /  ( glm::length(h) * glm::length(normal) ); // they re both normalised but just in case
        //clamp cosine
        cosine_h_n = glm::max(0.0f , cosine_h_n );
        if( least_cosine > cosine_h_n )
        {
            least_cosine = cosine_h_n;
        } 
        specular =  glm::pow( cosine_h_n , phong_exponent   ) *  glm::vec3( specular.x * scene.point_lights[i].intensity.x , specular.y * scene.point_lights[i].intensity.y  , specular.z * scene.point_lights[i].intensity.z );  
        //clamp diffuse
        specular.x = glm::min(255.0f , specular.x );
        specular.y = glm::min(255.0f , specular.y );
        specular.z = glm::min(255.0f , specular.z );
        
        color += diffuse  + specular  ;   


    }
    //std::cout << " leats cosine  " << least_cosine << std::endl; 
    //clamp color
    color.x = glm::min( color.x , 255.0f );
    color.y = glm::min( color.y , 255.0f );
    color.z = glm::min( color.z , 255.0f );

    //clamp to nearest integer
    color.x = (float) ( (int ) color.x );
    color.y = (float) ( (int ) color.y );
    color.z = (float) ( (int ) color.z );

    return color; 

    

}


