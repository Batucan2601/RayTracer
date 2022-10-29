#pragma once
#include "../Header/Prototypes.h"
#include "../Dependencies/glm_headers/glm.hpp"
#include "../Header/Ray.h"
#include "../Header/Helper.h"
// the basic ray tracing block 
int max_recursion_depth; 
static void generate_image(parser::Scene & scene)
{
    int i = 0; 
    
    std::cout <<  "camera size  " << scene.cameras.size()  << std::endl; 
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
                //initilaize recursion depth
                max_recursion_depth = scene.max_recursion_depth;
                // initialize the ray
                glm::vec3 current_pixel_world_space = starting_point_glm + glm::vec3( interval_row.x ,  interval_row.y ,  interval_row.z) * (float) x + (float ) y * glm::vec3( interval_col.x ,  interval_col.y ,  interval_col.z) ; 
                Ray ray(glm::vec3( current_camera.position.x , current_camera.position.y , current_camera.position.z  )  , current_pixel_world_space );
                glm::vec3  color = color_pixel(scene  , ray);
                
                // color cast 
                image[i++] = (unsigned char) (color.x);
                image[i++] = (unsigned char) (color.y);
                image[i++] = (unsigned char) (color.z);
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
    
    glm::vec3 top_left_corner = intersection + left_vec * right + cam_up * top;
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
    material.is_conductor = false; 
    material.is_dielectric = false; 
    material.is_mirror = false; 

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
    for (size_t i = 0; i < scene.point_lights.size(); i++) 
    {
        //each light pos 
        glm::vec3 light_pos = glm::vec3(scene.point_lights[i].position.x , scene.point_lights[i].position.y ,scene.point_lights[i].position.z);
        
        //add shadow ray epsilon in direction of normal 
        hit_point += normal * scene.shadow_ray_epsilon;
        // cast  shadow ray
        Ray shadow_ray(hit_point , light_pos );

        
        // else you directly went to a light
        
        glm::vec3 hit_point_temp; 
        glm::vec3 normal_temp;
        parser::Material material_temp;  
        bool  is_object_hit = ray_object_intersection( shadow_ray , scene ,  hit_point_temp , normal_temp , material_temp , object_id , is_shadow_rays_active);
        
        if( is_object_hit ) // it is in shadow no contribution from light
        {
            if( glm::distance(hit_point_temp , hit_point) < glm::distance(hit_point, light_pos) ) //before light source 
            {
                continue; 
            }
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

        glm::vec3 viewpos = glm::normalize(ray.origin - hit_point); 
        glm::vec3 h =  glm::normalize( shadow_ray.direction + viewpos ); // this might be flawed 
        float cosine_h_n = glm::dot(h , normal) /  ( glm::length(h) * glm::length(normal) ); // they re both normalised but just in case
        //clamp cosine
        cosine_h_n = glm::max(0.0f , cosine_h_n );

        specular =  glm::pow( cosine_h_n , phong_exponent   ) *  glm::vec3( specular.x * scene.point_lights[i].intensity.x , specular.y * scene.point_lights[i].intensity.y  , specular.z * scene.point_lights[i].intensity.z ) / (glm::distance(light_pos , hit_point )*glm::distance(light_pos, hit_point)  );  
        //clamp diffuse
        specular.x = glm::min(255.0f , specular.x );
        specular.y = glm::min(255.0f , specular.y );
        specular.z = glm::min(255.0f , specular.z );
        
        color += diffuse  + specular  ;   


    }
    if(material.is_mirror ) // do recursion
    {
        if( max_recursion_depth > 0 )
        {
            max_recursion_depth -= 1; // max recursion depth is reduced 
            glm::vec3 view_pos = glm::vec3( ray.direction.x  * -1 , ray.direction.y  * -1 , ray.direction.z  * -1  ); 
            float cosine_omega  = glm::dot( view_pos , normal ) / ( glm::length( normal ) * glm::length(view_pos));
            glm::vec3 vr = - view_pos +  normal * cosine_omega * 2.0f ;  
            Ray recursion_ray( hit_point , vr );
            glm::vec3 recursion_color = color_pixel(scene , recursion_ray );
            color =  color + glm::vec3( recursion_color.x  * material.mirror.x ,  recursion_color.y  * material.mirror.y  ,  recursion_color.z  * material.mirror.z ); 
        }
        
    }
    else if( material.is_dielectric)
    {
        //for wt however...
        // assume refractive index of air is 1
        float n1 = 1.0f; 
        float n2 = material.refraction_index;
        glm::vec3 view_pos = glm::vec3( ray.direction.x  * -1 , ray.direction.y  * -1 , ray.direction.z  * -1  ); 
        float cosine_omega  = glm::dot( view_pos , normal ) / ( glm::length( normal ) * glm::length(view_pos));
        glm::vec3 vr = - view_pos +  normal * cosine_omega * 2.0f ;  
        float sine_omega = glm::sqrt(1 -  ( cosine_omega*cosine_omega) );
        float sine_phi = n1/n2 *  sine_omega;
        float cosine_phi = glm::sqrt(1 -  ( sine_phi*sine_phi) );
        glm::vec3 wt = (ray.direction + ( cosine_omega * normal )   ) * ( n1 / n2) -  ( normal * cosine_phi ) ;
        // compute reflection ratio
        float r_parallel  = (n2 * cosine_omega - n1 * cosine_phi) / ( n2* cosine_omega + n1 * cosine_phi); 
        float r_ortho =  (n1 * cosine_omega - n2 * cosine_phi) / (n1 * cosine_omega - n2 * cosine_phi);
        float reflection_ratio = 0.5f * ( r_parallel * r_parallel + r_ortho * r_ortho);

        // compute transmission ratio 
        float refraction_ratio = 1 - reflection_ratio; 
        
        Ray reflection_ray( hit_point , hit_point + vr  );
        glm::vec3 iterated_hit_point =  hit_point +  (ray.direction * 0.0001f); // iterate a little bit
        
        Ray refraction_ray( iterated_hit_point , iterated_hit_point + wt );
        glm::vec3 reflection_color(0.0f , 0.0f , 0.0f );
        glm::vec3 refraction_color(0.0f , 0.0f , 0.0f );
        int max_rec_for_ref_rec = max_recursion_depth;
        if( max_recursion_depth > 0 )
        {
            max_recursion_depth -= 1; 
            reflection_color =  reflection_ratio * color_pixel(scene ,reflection_ray);
            //refraction_color =  refraction_ratio * color_pixel(scene ,refraction_ray);

        }
        glm::vec3 second_hit_point(0.0f , 0.0f , 0.0f );
        glm::vec3 second_normal(0.0f ,0.0f , 0.0f );
        // calculate attenuation
        bool is_hit = calculate_second_hitpoint_in_same_object( scene , refraction_ray , hit_point , normal , object_id  , second_hit_point , second_normal);
        if ( is_hit )
        {
            float distance = sqrt( ( hit_point.x - second_hit_point.x) * ( hit_point.x - second_hit_point.x)  + ( hit_point.y - second_hit_point.y) * ( hit_point.y - second_hit_point.y) + ( hit_point.z - second_hit_point.z) * ( hit_point.z - second_hit_point.z)    ); 
            // now exit the medium 
            // now n1 and n2 are swaped
            float temp_n = n2;
            n2 = n1;
            n1 = temp_n;

            glm::vec3 incoming_ray =  refraction_ray.direction;
            float cosine_omega  = glm::dot( incoming_ray , normal ) / ( glm::length( normal ) * glm::length(incoming_ray));
            float sine_omega = glm::sqrt(1 -  ( cosine_omega*cosine_omega) );
            float sine_phi = n1/n2 *  sine_omega;
            float cosine_phi = glm::sqrt(1 -  ( sine_phi*sine_phi) );
            glm::vec3 wt = (incoming_ray + (cosine_omega * normal)  ) * ( n1 / n2) -  ( normal * cosine_phi ) ;
            Ray out_ray( second_hit_point + second_normal * ( 0.01f) ,  second_hit_point + wt );

            if( max_recursion_depth > 0 )
            {
                max_recursion_depth -= 1; 
                refraction_color = refraction_ratio *  color_pixel(scene , out_ray  ); // L(0)
                glm::vec3 c  = glm::vec3( material.absorptionCoefficient.x ,  material.absorptionCoefficient.y  ,  material.absorptionCoefficient.z )  ;
                float e = 2.7182818f;
                refraction_color.x = refraction_color.x * powf32(e , -1 * c.x * distance  );
                refraction_color.y = refraction_color.y * powf32(e , -1 * c.y * distance  );
                refraction_color.z = refraction_color.z * powf32(e , -1 * c.z * distance  );
            }
        }
        color +=  reflection_color + refraction_color;

    }
    else if( material.is_conductor)
    {
        
        //for wt however...
        // assume rfeeractive index of air is 1
        glm::vec3  k = glm::vec3( material.absorptionCoefficient.x , material.absorptionCoefficient.y , material.absorptionCoefficient.z  ); ; 
        float n1 = 1.0f; 
        float n2 = material.refraction_index;
        glm::vec3 view_pos = glm::vec3( ray.direction.x  * -1 , ray.direction.y  * -1 , ray.direction.z  * -1  ); 
        float cosine_omega  = glm::dot( view_pos , normal ) / ( glm::length( normal ) * glm::length(view_pos));
        glm::vec3 vr = - view_pos +  normal * cosine_omega * 2.0f ;  
        float sine_omega = glm::sqrt(1 -  ( cosine_omega*cosine_omega) );
        float sine_phi = n1/n2 *  sine_omega;
        float cosine_phi = glm::sqrt(1 -  ( sine_phi*sine_phi) );
        glm::vec3 wt = (ray.direction + cosine_omega * normal  ) * ( n1 / n2) - normal * cosine_phi;
        // compute reflection ratio
        float r_s  =  ( ( n2 * n2 + glm::length(k) * glm::length(k) ) - 2 * n2 * cosine_omega + (cosine_omega * cosine_omega ) )  / ( ( n2 * n2 + glm::length(k) * glm::length(k) ) + 2 * n2 * cosine_omega + (cosine_omega * cosine_omega ));  
        float r_p =   ( ( n2 * n2 + glm::length(k) * glm::length(k) ) * (cosine_omega * cosine_omega )  - 2 * n2 * cosine_omega + 1  ) / (( n2 * n2 + glm::length(k) * glm::length(k) ) * (cosine_omega * cosine_omega )  + 2 * n2 * cosine_omega + 1 );
        float reflection_ratio = 0.5f * ( r_s + r_p);

        // compute transmission ratio 
        if( max_recursion_depth > 0 )
        {
            max_recursion_depth -= 1; 
            Ray reflection_ray( hit_point , hit_point + vr  );
            glm::vec3 iterated_hit_point =  hit_point + (normal * 0.01f); // iterate a little bit
            glm::vec3 reflection_color =  reflection_ratio * color_pixel(scene ,reflection_ray);
            reflection_color = glm::vec3( reflection_color.x * material.mirror.x , reflection_color.y * material.mirror.y , reflection_color.z * material.mirror.z);
            color += reflection_color; 
            
        }
       

    }
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


