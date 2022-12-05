#pragma once
#define STB_IMAGE_WRITE_IMPLEMENTATION // stb 
#include "../System_Files/support_files/stb_image_write.h"
#include "../System_Files/support_files/parser.h"
#include "../Header/Prototypes.h"
#include "../Header/Ray.h"
#include "../Header/BVH.h"
#include <random>
static parser::Vec3f color_pixel_BVH(parser::Scene& scene , Ray & ray, BVH::BoundingBoxTree *&  node  );

// the basic ray tracing block 
int max_recursion_depth; 
// uniform dist for area light
std::random_device rd_area_light;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_area(rd_area_light()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_area_light(-0.5 , 0.5);
static void generate_image(parser::Scene & scene)
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0 , 1.0);
    
    std::cout <<  "camera size  " << scene.cameras.size()  << std::endl; 
    for (size_t camera_no = 0; camera_no < scene.cameras.size() ; camera_no++)
    {
        int i = 0; 

        current_camera = scene.cameras[camera_no];  // get the current_camera

        int width = current_camera.image_width, height = current_camera.image_height;

        //calculate the starting point of  (0 , 0 ) in world space
        parser::Vec3f starting_point; // up left  (0 , 0 )
        parser::Vec3f interval_row; // real space interval in two consec pixels in same row 
        parser::Vec3f interval_col;// real space interval in two consec pixels in same  col
        calculate_image_plane(current_camera , starting_point , interval_row , interval_col );
        //we now have a startin point
        parser::Vec3f starting_point_parser = parser::Vec3f( starting_point.x , starting_point.y , starting_point.z );  
        //std::cout << "starting point " <<  starting_point_parser.x << " " << starting_point_parser.y  << starting_point_parser.z << std::endl; 
        //std::cout << "starting point " <<  interval_row.x << " " << interval_row.y  << interval_row.z << std::endl; 
        //std::cout << "starting point " <<  interval_col.x << " " << interval_col.y  << interval_col.z << std::endl; 

        // HW3
        int num_samples = current_camera.number_of_samples; 
        int num_samples_sqrt = std::sqrt(num_samples); // we will use this
        unsigned char* image = new unsigned char [width * height * 3  ];
        //HW2
        //create bvh tree
        //BVH::BoundingBoxTree boxtree;
       //BVH::BoundingBoxTree* head_node = &boxtree; 
        //BVH::generate_BVH_tree(scene ,head_node );
        parser::Vec3f  color;
        color.x = 0.0f;
        color.y = 0.0f;
        color.z = 0.0f;
        for (size_t y = 0; y < height; y++)
        {
            for (size_t x = 0; x < width; x++)
            {
                
                //initilaize recursion depth
                
                // initialize the ray
                parser::Vec3f current_pixel_world_space =  (starting_point_parser +  ( parser::Vec3f( interval_row.x ,  interval_row.y ,  interval_row.z) * (float)x )  )  +   ( parser::Vec3f( interval_col.x ,  interval_col.y ,  interval_col.z) * (float )y  )  ; 
                
                if( num_samples_sqrt > 1 )
                {
                    for (size_t j = 0; j < num_samples_sqrt ; j++)
                    {

                        for (size_t k = 0; k < num_samples_sqrt ; k++)
                        {
                            // new boundary 
                            parser::Vec3f sub_pixel = current_pixel_world_space +  ( (interval_row/num_samples_sqrt ) * (float)(k+1) ) + ( (interval_col/num_samples_sqrt ) * (float) ( j +1)  );
                            parser::Vec3f prev_pixel = current_pixel_world_space +  ( (interval_row/num_samples_sqrt ) * (float)k ) + ( (interval_col/num_samples_sqrt ) * (float)j );

                            // choose a point between them
                            float alpha = dis(gen); 
                            // pixel = alpha * ( prev_pixel (our pixel) ) + (1-alpha ) * sub_pixel
                            parser::Vec3f stratified_pixel =   ( ( prev_pixel * alpha ) +  (sub_pixel * (1-alpha) )  ) ; 

                            //refresh recursion depth
                            max_recursion_depth = scene.max_recursion_depth;

                            Ray ray(parser::Vec3f( current_camera.position.x , current_camera.position.y , current_camera.position.z  )  , stratified_pixel );
                            ray.sample_no = j * k + k; 
                            ray.generate_time();
                            color = color +  color_pixel(scene  , ray);

                        }
                        
                    }
                }

                else
                {
                    //refresh recursion depth
                    max_recursion_depth = scene.max_recursion_depth;
                    Ray ray(parser::Vec3f( current_camera.position.x , current_camera.position.y , current_camera.position.z  )  , current_pixel_world_space );
                    color = color +  color_pixel(scene  , ray);
                }
                
                color = color / num_samples ;
                color.x = std::min( 249.0f , color.x );
                color.y = std::min( 249.0f , color.y );
                color.z = std::min( 249.0f , color.z );

                image[i ++ ] = (unsigned char) (color.x);
                image[i ++ ] = (unsigned char) (color.y);
                image[i ++ ] = (unsigned char) (color.z);
               /*  //initilaize recursion depth
                max_recursion_depth = scene.max_recursion_depth;
                // initialize the ray
                parser::Vec3f current_pixel_world_space =  (starting_point_parser +  ( parser::Vec3f( interval_row.x ,  interval_row.y ,  interval_row.z) * (float)x )  )  +   ( parser::Vec3f( interval_col.x ,  interval_col.y ,  interval_col.z) * (float )y  )  ; 
                Ray ray(parser::Vec3f( current_camera.position.x , current_camera.position.y , current_camera.position.z  )  , current_pixel_world_space );
                parser::Vec3f  color = color_pixel(scene  , ray);
                // color cast 
                image[i++] = (unsigned char) (color.x);
                image[i++] = (unsigned char) (color.y);
                image[i++] = (unsigned char) (color.z);
               */
                // color cast 
            
            }
        }
        //BVH::delete_boxes(head_node);
        //save the written image
        //stbi_write_png(current_camera.image_name.c_str() , width * num_samples_sqrt, height * num_samples_sqrt,3 , image , width * num_samples_sqrt *  3  );
        stbi_write_png(current_camera.image_name.c_str() , width , height,3 , image , width  *  3  );
        
        
        //in the end just reallocate the image
        delete[] image; 
    }
}

static void calculate_image_plane(const parser::Camera &current_camera , parser::Vec3f & starting_point,  parser::Vec3f & interval_row, parser::Vec3f & interval_col)
{
    //calculate right vector
    parser::Vec3f cam_up = parser::normalize(parser::Vec3f(current_camera.up.x,current_camera.up.y,current_camera.up.z));
    parser::Vec3f cam_gaze = parser::normalize(parser::Vec3f(current_camera.gaze.x,current_camera.gaze.y,current_camera.gaze.z));

    
    parser::Vec3f right_vec = parser::normalize(parser::cross(cam_gaze , cam_up ));
    //find the intersection with the image
    parser::Vec3f intersection =  parser::Vec3f(current_camera.position.x , current_camera.position.y , current_camera.position.z   )  +  cam_gaze *  current_camera.near_distance ;

    //now get to the (0 ,0 ) 
    
    // but we cannot use right we need left
    parser::Vec3f left_vec = parser::Vec3f( -1 * right_vec.x , -1 * right_vec.y , -1 * right_vec.z  );

    // check how much we need to go left
    float left = current_camera.near_plane.x;
    float right = current_camera.near_plane.y;
    float bottom = current_camera.near_plane.z;
    float top = current_camera.near_plane.w;
    
    parser::Vec3f top_left_corner = intersection + left_vec * right + cam_up * top;
    starting_point.x = top_left_corner.x;
    starting_point.y = top_left_corner.y;
    starting_point.z = top_left_corner.z;
    

    //lastly calculate the intervals

    parser::Vec3f most_left = intersection + right_vec * left; 
    parser::Vec3f most_right = intersection + right_vec * right;

    

    parser::Vec3f most_top = intersection +  cam_up * top ; 
    parser::Vec3f most_bottom = intersection + cam_up * bottom ;

    #ifdef DEBUG
    std::cout << "most right " <<  most_right.x  << " " << most_right.y <<" " << most_right.z << std::endl;
    std::cout << "most left " <<  most_left.x  << " " << most_left.y <<" " << most_left.z << std::endl;
    std::cout << "most top " <<  most_top.x  << " " << most_top.y <<" " << most_top.z << std::endl;
    std::cout << "most bottom " <<  most_bottom.x  << " " << most_bottom.y <<" " << most_bottom.z << std::endl;
    #endif DEBUG
    parser::Vec3f interval_row_parser  = (most_right - most_left  );
    std::cout << "interval row " << interval_row_parser.x <<  " " << interval_row_parser.y << " " << interval_row_parser.z << std::endl;
    
    interval_row_parser = interval_row_parser / (float)current_camera.image_width;
    std::cout << "interval row " << interval_row_parser.x <<  " " << interval_row_parser.y << " " << interval_row_parser.z << std::endl;
    
    parser::Vec3f interval_col_parser = ( most_bottom - most_top );
    interval_col_parser = interval_col_parser / (float)current_camera.image_height;

    #ifdef DEBUG
    std::cout << "interval row " << interval_row_parser.x <<  " " << interval_row_parser.y << " " << interval_row_parser.z << std::endl;
    std::cout << "interval col " << interval_col_parser.x <<  " " << interval_col_parser.y << " " << interval_col_parser.z << std::endl;
    #endif DEBUG

    interval_row.x = interval_row_parser.x;
    interval_row.y = interval_row_parser.y;
    interval_row.z = interval_row_parser.z;

    interval_col.x = interval_col_parser.x;
    interval_col.y = interval_col_parser.y;
    interval_col.z = interval_col_parser.z;
}

static parser::Vec3f color_pixel(parser::Scene& scene , Ray & ray )
{

    
    parser::Vec3f color(0.0f , 0.0f , 0.0f );
    
    parser::Vec3f hit_point; 
    parser::Face hit_face; 

    parser::Vec3f normal;
    parser::Material material; 
    material.is_conductor = false; 
    material.is_dielectric = false; 
    material.is_mirror = false; 
    
    if( current_camera.aperture_size != -1 )
    {
        //find  focal point p
        parser::Vec3f focal_point_p = ray.origin  + ray.direction * current_camera.focus_distance; 
        // assume your ray origin is middle of aperture
        //generate two random numbers
        float random_num_u = current_camera.aperture_size * dis_area_light(gen_area); //right 
        float random_num_v = current_camera.aperture_size *  dis_area_light(gen_area); //up 
        //recalculate ray origin
        ray.origin = ray.origin +  parser::cross(current_camera.gaze , current_camera.up) * random_num_u + current_camera.up * random_num_v ;
        ray.direction = parser::normalize( focal_point_p - ray.origin);
    }
    
    //find the hitpoint 
    bool is_shadow_rays_active = false;
    int object_id = -1; 
    bool  is_object_hit = ray_object_intersection( ray , scene ,  hit_point , normal  , material   , object_id , hit_face, is_shadow_rays_active);

    if( !is_object_hit)
    {
        return parser::Vec3f(scene.background_color.x, scene.background_color.y , scene.background_color.z );
    }
    //normalize normal 
    //now we got the  nearest hitpoint and normal of that hitpoint. we can calculate color and cast shadow rays
    
    //we can add ambient shading
    color = parser::Vec3f(scene.ambient_light.x , scene.ambient_light.y , scene.ambient_light.z);
    //  for each point light
    is_shadow_rays_active = true; 
    for (size_t i = 0; i < scene.point_lights.size(); i++) 
    {
        //each light pos 
        parser::Vec3f light_pos = parser::Vec3f(scene.point_lights[i].position.x , scene.point_lights[i].position.y ,scene.point_lights[i].position.z);
        
        //add shadow ray epsilon in direction of normal 
        hit_point = hit_point +  normal * scene.shadow_ray_epsilon;
        // cast  shadow ray
        Ray shadow_ray(hit_point , light_pos );
        shadow_ray.time_parameter = ray.time_parameter;
        
        // else you directly went to a light
        
        parser::Vec3f hit_point_temp; 
        parser::Face hit_face_temp; 
        parser::Vec3f normal_temp;
        parser::Material material_temp;  
        bool  is_object_hit = ray_object_intersection( shadow_ray , scene ,  hit_point_temp , normal_temp , material_temp , object_id , hit_face_temp ,  is_shadow_rays_active);
        
        if( is_object_hit ) // it is in shadow no contribution from light
        {
            if( parser::distance(hit_point_temp , hit_point) < parser::distance(hit_point, light_pos) ) //before light source 
            {
                continue; 
            }
        }
        
        // diffuse 
        parser::Vec3f diffuse = parser::Vec3f( material.diffuse.x , material.diffuse.y, material.diffuse.z);
        float cosine_alpha = parser::dot(shadow_ray.direction  ,  normal ) / (parser::length(shadow_ray.direction) * parser::length(normal)  ) ; //normal vectors
        //clamp cosine alpha
        cosine_alpha = std::max(0.0f , cosine_alpha);
        diffuse = diffuse * cosine_alpha; 
        diffuse = parser::Vec3f( diffuse.x * scene.point_lights[i].intensity.x , diffuse.y *  scene.point_lights[i].intensity.y ,diffuse.z *   scene.point_lights[i].intensity.z) / (parser::distance(light_pos , hit_point )*parser::distance(light_pos, hit_point)  );
        //clamp diffuse
        diffuse.x = std::min(255.0f , diffuse.x );
        diffuse.y = std::min(255.0f , diffuse.y );
        diffuse.z = std::min(255.0f , diffuse.z );


        // specular 
        // specular 
        parser::Vec3f specular = parser::Vec3f( material.specular.x , material.specular.y, material.specular.z);
        float phong_exponent = material.phong_exponent; 

        parser::Vec3f viewpos = parser::normalize(ray.origin - hit_point); 
        parser::Vec3f h =  parser::normalize( shadow_ray.direction + viewpos ); // this might be flawed 
        float cosine_h_n = parser::dot(h , normal) /  ( parser::length(h) * parser::length(normal) ); // they re both normalised but just in case
        //clamp cosine
        cosine_h_n = std::max(0.0f , cosine_h_n );

        specular =    parser::Vec3f( specular.x * scene.point_lights[i].intensity.x , specular.y * scene.point_lights[i].intensity.y  , specular.z * scene.point_lights[i].intensity.z ) / (parser::distance(light_pos , hit_point )*parser::distance(light_pos, hit_point)  ) * std::pow( cosine_h_n , phong_exponent   );  
        //clamp diffuse
        specular.x = std::min(255.0f , specular.x );
        specular.y = std::min(255.0f , specular.y );
        specular.z = std::min(255.0f , specular.z );
        
         
        color =  color + diffuse  + specular  ;   


    }
    // for each area light
    is_shadow_rays_active = true; 
    for (size_t i = 0; i < scene.area_lights.size(); i++)
    {
        
        float random_num_u = scene.area_lights[i].size * dis_area_light(gen_area);
        float random_num_v = scene.area_lights[i].size *  dis_area_light(gen_area);
        parser::Vec3f p = scene.area_lights[i].position + scene.area_lights[i].ortho_basis.u * random_num_u  +  scene.area_lights[i].ortho_basis.v * random_num_v;  
    
        parser::Vec3f l = parser::normalize( hit_point - p );
        
        float cosine_theta = parser::dot( l , scene.area_lights[i].normal) /   (parser::length(scene.area_lights[i].normal) * parser::length(l) );
        if( cosine_theta < 0 )
        {
            cosine_theta = std::abs( cosine_theta); //sin beta
            cosine_theta =  std::sqrt(1 - (cosine_theta * cosine_theta) );
            //cosine_theta = std::sqrt(1 - cosine_theta * cosine_theta);
        }
        //clamp cosine alpha
        //cosine_theta = std::abs( cosine_theta);
        parser::Vec3f L =  scene.area_lights[i].radiance *  (cosine_theta  * scene.area_lights[i].size * scene.area_lights[i].size   )/ (parser::distance(p , hit_point )*parser::distance(p, hit_point));        
        

        //add shadow ray epsilon in direction of normal 
        hit_point = hit_point +  normal * scene.shadow_ray_epsilon;
        // cast  shadow ray
        Ray shadow_ray(hit_point , p );
        shadow_ray.time_parameter = ray.time_parameter;

        
        // else you directly went to a light
        
        parser::Vec3f hit_point_temp; 
        parser::Vec3f normal_temp;
        parser::Material material_temp;  
        parser::Face hit_face_temp;
        

        bool  is_object_hit = ray_object_intersection( shadow_ray , scene ,  hit_point_temp , normal_temp , material_temp , object_id , hit_face_temp , is_shadow_rays_active);
        
        //before that calculate t of light pos //also pass  
        float light_pos_t = (p.x - shadow_ray.origin.x )  / shadow_ray.direction.x;
        //calculate hitpoint t 
        float hit_point_t = (hit_point_temp.x - shadow_ray.origin.x )  / shadow_ray.direction.x;

        if( is_object_hit && ( light_pos_t > hit_point_t) ) // it is in shadow no contribution from light
        {
            if( parser::distance(hit_point_temp , hit_point) < parser::distance(hit_point, p) ) //before light source 
            {
                continue; 
            }
        }
        // diffuse 
        parser::Vec3f diffuse = parser::Vec3f( material.diffuse.x , material.diffuse.y, material.diffuse.z);
        float cosine_alpha = parser::dot(shadow_ray.direction  ,  normal ) / (parser::length(shadow_ray.direction) * parser::length(normal)  ) ; //normal vectors

        //clamp cosine alpha
        cosine_alpha = std::max(0.0f , cosine_alpha);
        diffuse = diffuse * cosine_alpha; 
        diffuse = parser::Vec3f( diffuse.x * L.x , diffuse.y *  L.y ,diffuse.z *   L.z);
        //clamp diffuse
        diffuse.x = std::min(255.0f , diffuse.x );
        diffuse.y = std::min(255.0f , diffuse.y );
        diffuse.z = std::min(255.0f , diffuse.z );


        // specular 
        // specular 
        parser::Vec3f specular = parser::Vec3f( material.specular.x , material.specular.y, material.specular.z);
        float phong_exponent = material.phong_exponent; 

        parser::Vec3f viewpos = parser::normalize(ray.origin - hit_point); 
        parser::Vec3f h =  parser::normalize( shadow_ray.direction + viewpos ); // this might be flawed 
        float cosine_h_n = parser::dot(h , normal) /  ( parser::length(h) * parser::length(normal) ); // they re both normalised but just in case
        //clamp cosine
        cosine_h_n = std::max(0.0f , cosine_h_n );

        specular =    parser::Vec3f( specular.x * L.x, specular.y * L.y  , specular.z * L.z )  * std::pow( cosine_h_n , phong_exponent   );  
        //clamp diffuse
        specular.x = std::min(255.0f , specular.x );
        specular.y = std::min(255.0f , specular.y );
        specular.z = std::min(255.0f , specular.z );
        
         
        color =  color + diffuse  + specular  ; 

    }
    
    if(material.is_mirror ) // do recursion
    {
        if( max_recursion_depth > 0 )
        {
            max_recursion_depth -= 1; // max recursion depth is reduced 
            parser::Vec3f W0 = parser::Vec3f( ray.direction.x  * -1 , ray.direction.y  * -1 , ray.direction.z  * -1  ); 
            float cosine_omega  = parser::dot( W0 , normal ) / ( parser::length( normal ) * parser::length(W0));
            parser::Vec3f Wr = W0 * -1.0f +  normal * cosine_omega * 2.0f ;  
            if( material.roughness != 0 )
            {
                parser::Vec3f Wr_ = Wr; 
                //calculate orthonormal bases
                // n'
                if( std::abs(Wr.x) < std::abs(Wr.y) && std::abs(Wr.x) < std::abs(Wr.z))
                {
                    // x is smallest
                    Wr_.x =  1.0f; 
                }
                else if( std::abs(Wr.y) < std::abs(Wr.x) && std::abs(Wr.y) < std::abs(Wr.z))
                {
                    // y is smallest
                    Wr_.y =  1.0f; 
                }
                else
                {
                    // z is smallest
                    Wr_.z =  1.0f; 
                }
                //compute u
                parser::Vec3f u = parser::normalize(parser::cross(Wr_ , Wr ));
                parser::Vec3f v = parser::normalize(parser::cross(Wr , u ));
                float random_num_u = material.roughness * dis_area_light(gen_area);
                float random_num_v = material.roughness * dis_area_light(gen_area);
                Wr_ = Wr + u * random_num_u + v * random_num_v; 
                Wr = Wr_; 
            }
            Ray recursion_ray( hit_point + Wr * 0.001f , hit_point + Wr );
            recursion_ray.generate_time();
            parser::Vec3f recursion_color = color_pixel(scene , recursion_ray );
            color =  color + parser::Vec3f( recursion_color.x  * material.mirror.x ,  recursion_color.y  * material.mirror.y  ,  recursion_color.z  * material.mirror.z ); 
        }
        
    }
    else if( material.is_dielectric)
    {
        //for wt however...
        // assume refractive index of air is 1
        float n1 = 1.0f; 
        float n2 = material.refraction_index;
        parser::Vec3f W0 = parser::Vec3f( ray.direction.x  * -1 , ray.direction.y  * -1 , ray.direction.z  * -1  ); // wr we want
        float cosine_omega  = parser::dot( W0 , normal ) / ( parser::length( normal ) * parser::length(W0));
        parser::Vec3f Wr = W0 *-1  +  normal * cosine_omega * 2.0f ;  
        float sine_omega = std::sqrt(1 -  ( cosine_omega*cosine_omega) );
        float sine_phi = n1/n2 *  sine_omega;
        float cosine_phi = (1 -   ( ( (n1/n2) * (n1/n2) )  *(1 -  (cosine_omega*cosine_omega )) )  );
        parser::Vec3f wt =  ((W0 * -1.0f  + ( normal * cosine_omega )   ) * ( n1 / n2) )  -  ( normal * cosine_phi ) ;
        // compute reflection ratio
        float r_parallel  = (n2 * cosine_omega - n1 * cosine_phi) / ( n2* cosine_omega + n1 * cosine_phi); 
        float r_ortho =  (n1 * cosine_omega - n2 * cosine_phi) / (n1 * cosine_omega - n2 * cosine_phi);
        float reflection_ratio = 0.5f * ( r_parallel * r_parallel + r_ortho * r_ortho);

        // compute transmission ratio 
        float refraction_ratio = 1 - reflection_ratio; 
        Ray reflection_ray( hit_point + Wr * 0.005f , hit_point + Wr  ); // real WR 
        Ray refraction_ray( hit_point + wt * 0.005f , hit_point + wt ); //real Wt 

        parser::Vec3f reflection_color(0.0f , 0.0f , 0.0f );
        parser::Vec3f refraction_color(0.0f , 0.0f , 0.0f );
        if( max_recursion_depth > 0 )
        {
            max_recursion_depth -= 1; 
            //calculate second hit 
            parser::Vec3f second_hit_point(0.0f , 0.0f , 0.0f );
            parser::Vec3f second_normal(0.0f ,0.0f , 0.0f );
            parser::Face second_face;

            // calculate attenuation
            bool is_hit = calculate_second_hitpoint_in_same_object( scene , refraction_ray , hit_point , normal , object_id  , second_hit_point , second_normal ,second_face ,  is_shadow_rays_active);
            if ( is_hit  && cosine_phi > 0 ) // if cosine phi less than 0 abort 
            {
                cosine_phi = std::sqrt(cosine_phi);
                float distance = sqrt( ( hit_point.x - second_hit_point.x) * ( hit_point.x - second_hit_point.x)  + ( hit_point.y - second_hit_point.y) * ( hit_point.y - second_hit_point.y) + ( hit_point.z - second_hit_point.z) * ( hit_point.z - second_hit_point.z)    ); 
                // now exit the medium 
                // now n1 and n2 are swaped
                float temp_n = n2;
                n2 = n1;
                n1 = temp_n;

                second_normal = second_normal * -1; //negate normal
                parser::Vec3f W0 =  refraction_ray.direction * -1;
                float cosine_omega  = parser::dot( W0 , second_normal ) / ( parser::length( second_normal ) * parser::length(W0));
                float sine_omega = std::sqrt(1 -  ( cosine_omega*cosine_omega) );
                float sine_phi  = n1/n2 *  sine_omega;
                float cosine_phi = (1 -  (n1/n2) * (n1/n2) *(1 -  cosine_omega*cosine_omega) );
                parser::Vec3f wt = ( W0 * -1.0f  + (second_normal * cosine_omega )  ) * ( n1 / n2) -  ( second_normal * cosine_phi ) ;
                Ray out_ray( second_hit_point + (wt * ( 0.005f)) ,  second_hit_point + wt );
                //std::cout << "hit point  " << hit_point.x << " " << hit_point.y << "  " << hit_point.z << std::endl;
                //std::cout << "second hit point  " << second_hit_point.x << " " << second_hit_point.y  << " " << out_ray.origin.z << std::endl; 
                if( cosine_phi  > 0 ) //  else abort the refraction ray 
                {
                    cosine_phi = std::sqrt(cosine_phi);
                    reflection_color =   color_pixel(scene ,reflection_ray) * reflection_ratio * 0.1f;
                    refraction_color = color_pixel(scene , out_ray  ) * refraction_ratio; // L(0)
                    // attenuation 
                
                    parser::Vec3f c  = parser::Vec3f( material.absorptionCoefficient.x ,  material.absorptionCoefficient.y  ,  material.absorptionCoefficient.z )  ;
                    float e = 2.7182818f;
                    refraction_color.x = refraction_color.x * powf32(e , -1 * c.x * distance  );
                    refraction_color.y = refraction_color.y * powf32(e , -1 * c.y * distance  );
                    refraction_color.z = refraction_color.z * powf32(e , -1 * c.z * distance  );
                    //std::cout << "refraction color   " << reflection_color.x << " " << reflection_color.y  << " " << reflection_color.z << std::endl; 
                }
                else
                {
                    reflection_color =   color_pixel(scene ,reflection_ray) ;
                }

            }
            else
            {
                    reflection_color =   color_pixel(scene ,reflection_ray) ;
            }
            //refraction_color =  color_pixel(scene ,refraction_ray) * refraction_ratio ;
        }
        
        color = color  + refraction_color + reflection_color;

    }
    else if( material.is_conductor)
    {
        //for wt however...
        // assume rfeeractive index of air is 1
        float k = material.absorption_index ; 
        float n1 = 1.0f; 
        float n2 = material.refraction_index;
        parser::Vec3f view_pos = parser::Vec3f( ray.direction.x  * -1 , ray.direction.y  * -1 , ray.direction.z  * -1  ); 
        float cosine_omega  = parser::dot( view_pos , normal ) / ( parser::length( normal ) * parser::length(view_pos));
        parser::Vec3f vr =  view_pos * -1  +  normal * cosine_omega * 2.0f ;  
        float sine_omega = std::sqrt(1 -  ( cosine_omega*cosine_omega) );
        float sine_phi = n1/n2 *  sine_omega;
        float cosine_phi = std::sqrt(1 -  ( sine_phi*sine_phi) );
        parser::Vec3f wt = (ray.direction +  normal * cosine_omega  ) * ( n1 / n2) -  ( normal * cosine_phi ) ;
        // compute reflection ratio
        float r_s  =  ( ( n2 * n2 + k * k ) - 2 * n2 * cosine_omega + (cosine_omega * cosine_omega ) )  / ( ( n2 * n2 + k * k ) + 2 * n2 * cosine_omega + (cosine_omega * cosine_omega ));  
        float r_p =   ( ( n2 * n2 + k * k ) * (cosine_omega * cosine_omega )  - 2 * n2 * cosine_omega + 1  ) / (( n2 * n2 + k * k ) * (cosine_omega * cosine_omega )  + 2 * n2 * cosine_omega + 1 );
        float reflection_ratio = 0.5f * ( r_s + r_p);

        // compute transmission ratio 
        if( max_recursion_depth > 0 )
        {
            if( material.roughness != 0 )
            {
                parser::Vec3f Wr_ = vr; 
                //calculate orthonormal bases
                // n'
                if( std::abs(vr.x) < std::abs(vr.y) && std::abs(vr.x) < std::abs(vr.z))
                {
                    // x is smallest
                    Wr_.x =  1.0f; 
                }
                else if( std::abs(vr.y) < std::abs(vr.x) && std::abs(vr.y) < std::abs(vr.z))
                {
                    // y is smallest
                    Wr_.y =  1.0f; 
                }
                else
                {
                    // z is smallest
                    Wr_.z =  1.0f; 
                }
                //compute u
                parser::Vec3f u = parser::normalize(parser::cross(Wr_ , vr ));
                parser::Vec3f v = parser::normalize(parser::cross(vr , u ));
                float random_num_u = material.roughness * dis_area_light(gen_area);
                float random_num_v = material.roughness * dis_area_light(gen_area);
                Wr_ = vr + u * random_num_u + v * random_num_v; 
                vr = Wr_; 
            }
            max_recursion_depth -= 1; 
            Ray reflection_ray( hit_point +  ( vr * (0.01f) )  , hit_point + vr  );
            reflection_ray.generate_time();
            parser::Vec3f reflection_color =   color_pixel(scene ,reflection_ray) * reflection_ratio ;
            reflection_color = parser::Vec3f( reflection_color.x * material.mirror.x , reflection_color.y * material.mirror.y , reflection_color.z * material.mirror.z);
            color = color + reflection_color; 
        }
       

    }
    //clamp color
    color.x = std::min( color.x , 255.0f );
    color.y = std::min( color.y , 255.0f );
    color.z = std::min( color.z , 255.0f );

    //clamp to nearest integer
    color.x = (float) ( (int ) color.x );
    color.y = (float) ( (int ) color.y );
    color.z = (float) ( (int ) color.z );

    return color; 

    

}

/*
static parser::Vec3f color_pixel_BVH(parser::Scene& scene , Ray & ray, BVH::BoundingBoxTree *&  node  )
{
     
    parser::Vec3f color(0.0f , 0.0f , 0.0f );
    
    parser::Vec3f hit_point; 
    parser::Vec3f normal;
    parser::Material material; 
    material.is_conductor = false; 
    material.is_dielectric = false; 
    material.is_mirror = false; 


    //find the hitpoint 
    bool is_shadow_rays_active = false;
    int object_id = -1; 
    std::vector<std::pair<int,float> >  object_no_t_pair; 
    BVH::BoundingBox box = BVH::traverse_tree(node , ray , object_no_t_pair );
    if( box.object_id == -1 ) // not hit 
    {
        return parser::Vec3f(scene.background_color.x, scene.background_color.y , scene.background_color.z );
    }
    //else
    //check all ofthe objects it has hit 
    std::pair<int , float > currentpair;
    currentpair.second = INFINITY; 
    bool  is_object_hit = false; 
    
    parser::Vec3f temp_normal; 
    parser::Vec3f temp_hit_point; 
    parser::Material temp_material; 

    for (size_t i = 0; i < object_no_t_pair.size(); i++)
    {
        bool is_object_hit_loop = ray_object_intersection( ray , scene ,  temp_hit_point , temp_normal  , temp_material   , object_id ,  is_shadow_rays_active , object_no_t_pair[i].first);
        if( is_object_hit_loop )
        {
            is_object_hit = is_object_hit_loop;
            if(object_no_t_pair[i].second < currentpair.second  )
            {
                currentpair.second = object_no_t_pair[i].second;
                currentpair.first = object_no_t_pair[i].first;
                hit_point = temp_hit_point;
                normal = temp_normal;
                material = temp_material;
            }
        }:
    
     
    
    if( !is_object_hit)
    {
        return parser::Vec3f(scene.background_color.x, scene.background_color.y , scene.background_color.z );
    }
   
    
    //normalize normal 
    //now we got the  nearest hitpoint and normal of that hitpoint. we can calculate color and cast shadow rays
    
    //we can add ambient shading
    color = parser::Vec3f(scene.ambient_light.x , scene.ambient_light.y , scene.ambient_light.z);
    //  for each light
    is_shadow_rays_active = true; 
    for (size_t i = 0; i < scene.point_lights.size(); i++) 
    {
        //each light pos 
        parser::Vec3f light_pos = parser::Vec3f(scene.point_lights[i].position.x , scene.point_lights[i].position.y ,scene.point_lights[i].position.z);
        
        //add shadow ray epsilon in direction of normal 
        hit_point = hit_point +  normal * scene.shadow_ray_epsilon;
        // cast  shadow ray
        Ray shadow_ray(hit_point , light_pos );

        
        // else you directly went to a light
        
        parser::Vec3f hit_point_temp; 
        parser::Vec3f normal_temp;
        parser::Material material_temp;  
         bool  is_object_hit = ray_object_intersection( shadow_ray , scene ,  hit_point_temp , normal_temp , material_temp , object_id , is_shadow_rays_active);
        
        if( is_object_hit ) // it is in shadow no contribution from light
        {
            if( parser::distance(hit_point_temp , hit_point) < parser::distance(hit_point, light_pos) ) //before light source 
            {
                continue; 
            }
        }
        
        
        // diffuse 
        parser::Vec3f diffuse = parser::Vec3f( material.diffuse.x , material.diffuse.y, material.diffuse.z);
        float cosine_alpha = parser::dot(shadow_ray.direction  ,  normal ) / (parser::length(shadow_ray.direction) * parser::length(normal)  ) ; //normal vectors
        //clamp cosine alpha
        cosine_alpha = std::max(0.0f , cosine_alpha);
        diffuse = diffuse * cosine_alpha; 
        diffuse = parser::Vec3f( diffuse.x * scene.point_lights[i].intensity.x , diffuse.y *  scene.point_lights[i].intensity.y ,diffuse.z *   scene.point_lights[i].intensity.z) / (parser::distance(light_pos , hit_point )*parser::distance(light_pos, hit_point)  );
        //clamp diffuse
        diffuse.x = std::min(255.0f , diffuse.x );
        diffuse.y = std::min(255.0f , diffuse.y );
        diffuse.z = std::min(255.0f , diffuse.z );

        // specular 
        parser::Vec3f specular = parser::Vec3f( material.specular.x , material.specular.y, material.specular.z);
        float phong_exponent = material.phong_exponent; 

        parser::Vec3f viewpos = parser::normalize(ray.origin - hit_point); 
        parser::Vec3f h =  parser::normalize( shadow_ray.direction + viewpos ); // this might be flawed 
        float cosine_h_n = parser::dot(h , normal) /  ( parser::length(h) * parser::length(normal) ); // they re both normalised but just in case
        //clamp cosine
        cosine_h_n = std::max(0.0f , cosine_h_n );

        specular =    parser::Vec3f( specular.x * scene.point_lights[i].intensity.x , specular.y * scene.point_lights[i].intensity.y  , specular.z * scene.point_lights[i].intensity.z ) / (parser::distance(light_pos , hit_point )*parser::distance(light_pos, hit_point)  ) * std::pow( cosine_h_n , phong_exponent   );  
        //clamp diffuse
        specular.x = std::min(255.0f , specular.x );
        specular.y = std::min(255.0f , specular.y );
        specular.z = std::min(255.0f , specular.z );
        
        if( !material.is_mirror)
        {
        color =  color + diffuse  + specular  ;   
        }


    }
    if(material.is_mirror ) // do recursion
    {
        if( max_recursion_depth > 0 )
        {
            max_recursion_depth -= 1; // max recursion depth is reduced 
            parser::Vec3f W0 = parser::Vec3f( ray.direction.x  * -1 , ray.direction.y  * -1 , ray.direction.z  * -1  ); 
            float cosine_omega  = parser::dot( W0 , normal ) / ( parser::length( normal ) * parser::length(W0));
            parser::Vec3f Wr = W0 * -1.0f +  normal * cosine_omega * 2.0f ;  
            Ray recursion_ray( hit_point + Wr * 0.001f , hit_point + Wr );
            parser::Vec3f recursion_color = color_pixel(scene , recursion_ray );
            color =  color + parser::Vec3f( recursion_color.x  * material.mirror.x ,  recursion_color.y  * material.mirror.y  ,  recursion_color.z  * material.mirror.z ); 
        }
        
    }
    else if( material.is_dielectric)
    {
        //for wt however...
        // assume refractive index of air is 1
        float n1 = 1.0f; 
        float n2 = material.refraction_index;
        parser::Vec3f W0 = parser::Vec3f( ray.direction.x  * -1 , ray.direction.y  * -1 , ray.direction.z  * -1  ); // wr we want
        float cosine_omega  = parser::dot( W0 , normal ) / ( parser::length( normal ) * parser::length(W0));
        parser::Vec3f Wr = W0 *-1  +  normal * cosine_omega * 2.0f ;  
        float sine_omega = std::sqrt(1 -  ( cosine_omega*cosine_omega) );
        float sine_phi = n1/n2 *  sine_omega;
        float cosine_phi = (1 -   ( ( (n1/n2) * (n1/n2) )  *(1 -  (cosine_omega*cosine_omega )) )  );
        parser::Vec3f wt =  ((W0 * -1.0f  + ( normal * cosine_omega )   ) * ( n1 / n2) )  -  ( normal * cosine_phi ) ;
        // compute reflection ratio
        float r_parallel  = (n2 * cosine_omega - n1 * cosine_phi) / ( n2* cosine_omega + n1 * cosine_phi); 
        float r_ortho =  (n1 * cosine_omega - n2 * cosine_phi) / (n1 * cosine_omega - n2 * cosine_phi);
        float reflection_ratio = 0.5f * ( r_parallel * r_parallel + r_ortho * r_ortho);

        // compute transmission ratio 
        float refraction_ratio = 1 - reflection_ratio; 
        Ray reflection_ray( hit_point + Wr * 0.005f , hit_point + Wr  ); // real WR 
        Ray refraction_ray( hit_point + wt * 0.005f , hit_point + wt ); //real Wt 

        parser::Vec3f reflection_color(0.0f , 0.0f , 0.0f );
        parser::Vec3f refraction_color(0.0f , 0.0f , 0.0f );
        if( max_recursion_depth > 0 )
        {
            max_recursion_depth -= 1; 
            //calculate second hit 
            parser::Vec3f second_hit_point(0.0f , 0.0f , 0.0f );
            parser::Vec3f second_normal(0.0f ,0.0f , 0.0f );
            // calculate attenuation
            bool is_hit = calculate_second_hitpoint_in_same_object( scene , refraction_ray , hit_point , normal , object_id  , second_hit_point , second_normal , is_shadow_rays_active);
            if ( is_hit  && cosine_phi > 0 ) // if cosine phi less than 0 abort 
            {
                cosine_phi = std::sqrt(cosine_phi);
                float distance = sqrt( ( hit_point.x - second_hit_point.x) * ( hit_point.x - second_hit_point.x)  + ( hit_point.y - second_hit_point.y) * ( hit_point.y - second_hit_point.y) + ( hit_point.z - second_hit_point.z) * ( hit_point.z - second_hit_point.z)    ); 
                // now exit the medium 
                // now n1 and n2 are swaped
                float temp_n = n2;
                n2 = n1;
                n1 = temp_n;

                second_normal = second_normal * -1; //negate normal
                parser::Vec3f W0 =  refraction_ray.direction * -1;
                float cosine_omega  = parser::dot( W0 , second_normal ) / ( parser::length( second_normal ) * parser::length(W0));
                float sine_omega = std::sqrt(1 -  ( cosine_omega*cosine_omega) );
                float sine_phi  = n1/n2 *  sine_omega;
                float cosine_phi = (1 -  (n1/n2) * (n1/n2) *(1 -  cosine_omega*cosine_omega) );
                parser::Vec3f wt = ( W0 * -1.0f  + (second_normal * cosine_omega )  ) * ( n1 / n2) -  ( second_normal * cosine_phi ) ;
                Ray out_ray( second_hit_point + (wt * ( 0.005f)) ,  second_hit_point + wt );
                //std::cout << "hit point  " << hit_point.x << " " << hit_point.y << "  " << hit_point.z << std::endl;
                //std::cout << "second hit point  " << second_hit_point.x << " " << second_hit_point.y  << " " << out_ray.origin.z << std::endl; 
                if( cosine_phi  > 0 ) //  else abort the refraction ray 
                {
                    cosine_phi = std::sqrt(cosine_phi);
                    reflection_color =   color_pixel(scene ,reflection_ray) * reflection_ratio * 0.1f;
                    refraction_color = color_pixel(scene , out_ray  ) * refraction_ratio; // L(0)
                    // attenuation 
                
                    parser::Vec3f c  = parser::Vec3f( material.absorptionCoefficient.x ,  material.absorptionCoefficient.y  ,  material.absorptionCoefficient.z )  ;
                    float e = 2.7182818f;
                    refraction_color.x = refraction_color.x * powf32(e , -1 * c.x * distance  );
                    refraction_color.y = refraction_color.y * powf32(e , -1 * c.y * distance  );
                    refraction_color.z = refraction_color.z * powf32(e , -1 * c.z * distance  );
                    //std::cout << "refraction color   " << reflection_color.x << " " << reflection_color.y  << " " << reflection_color.z << std::endl; 
                }
                else
                {
                    reflection_color =   color_pixel(scene ,reflection_ray) ;
                }

            }
            else
            {
                    reflection_color =   color_pixel(scene ,reflection_ray) ;
            }
            //refraction_color =  color_pixel(scene ,refraction_ray) * refraction_ratio ;
        }
        
        color = color  + refraction_color + reflection_color;

    }
    else if( material.is_conductor)
    {
        //for wt however...
        // assume rfeeractive index of air is 1
        float k = material.absorption_index ; 
        float n1 = 1.0f; 
        float n2 = material.refraction_index;
        parser::Vec3f view_pos = parser::Vec3f( ray.direction.x  * -1 , ray.direction.y  * -1 , ray.direction.z  * -1  ); 
        float cosine_omega  = parser::dot( view_pos , normal ) / ( parser::length( normal ) * parser::length(view_pos));
        parser::Vec3f vr =  view_pos * -1  +  normal * cosine_omega * 2.0f ;  
        float sine_omega = std::sqrt(1 -  ( cosine_omega*cosine_omega) );
        float sine_phi = n1/n2 *  sine_omega;
        float cosine_phi = std::sqrt(1 -  ( sine_phi*sine_phi) );
        parser::Vec3f wt = (ray.direction +  normal * cosine_omega  ) * ( n1 / n2) -  ( normal * cosine_phi ) ;
        // compute reflection ratio
        float r_s  =  ( ( n2 * n2 + k * k ) - 2 * n2 * cosine_omega + (cosine_omega * cosine_omega ) )  / ( ( n2 * n2 + k * k ) + 2 * n2 * cosine_omega + (cosine_omega * cosine_omega ));  
        float r_p =   ( ( n2 * n2 + k * k ) * (cosine_omega * cosine_omega )  - 2 * n2 * cosine_omega + 1  ) / (( n2 * n2 + k * k ) * (cosine_omega * cosine_omega )  + 2 * n2 * cosine_omega + 1 );
        float reflection_ratio = 0.5f * ( r_s + r_p);

        // compute transmission ratio 
        if( max_recursion_depth > 0 )
        {
            max_recursion_depth -= 1; 
            Ray reflection_ray( hit_point +  ( vr * (0.01f) )  , hit_point + vr  );
            parser::Vec3f reflection_color =   color_pixel(scene ,reflection_ray) * reflection_ratio ;
            reflection_color = parser::Vec3f( reflection_color.x * material.mirror.x , reflection_color.y * material.mirror.y , reflection_color.z * material.mirror.z);
            color = color + reflection_color; 
        }
       

    }
    //clamp color
    color.x = std::min( color.x , 255.0f );
    color.y = std::min( color.y , 255.0f );
    color.z = std::min( color.z , 255.0f );

    //clamp to nearest integer
    color.x = (float) ( (int ) color.x );
    color.y = (float) ( (int ) color.y );
    color.z = (float) ( (int ) color.z );

    return color; 

    

}
*/