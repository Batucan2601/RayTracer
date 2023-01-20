#pragma once
#define STB_IMAGE_WRITE_IMPLEMENTATION // stb 
#include "../System_Files/support_files/stb_image_write.h"
#include "../System_Files/support_files/parser.h"
#include "../Header/Prototypes.h"
#include "../Header/Ray.h"
#include "../Header/BVH.h"
#include "../Header/Tonemap.h"
#include "../Header/EnvironmentalLight.h"
#include "../Header/BRDF.h"
#include "../Header/PathTracing.h"


#include <random>
static parser::Vec3f color_pixel_BVH(parser::Scene& scene , Ray & ray, BVH::BoundingBoxTree *&  node  );

// the basic ray tracing block 
int max_recursion_depth; 
// uniform dist for area light
std::random_device rd_area_light;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_area(rd_area_light()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_area_light(-0.5 , 0.5);

//global intevals
parser::Vec3f global_interval_row; 
parser::Vec3f global_interval_col; 

//bounding boxes
std::vector< BVH::BoundingBox > bounding_boxes;
static void generate_image(parser::Scene & scene)
{
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0 , 1.0);
    
    std::cout <<  "camera size  " << scene.cameras.size()  << std::endl; 
    init_perlin_noise();
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
        float* tonemap_image = new  float[width * height * 3 ];
        std::vector<float> luminances; //100 luminances
        //HW2
        //create bvh tree
        //BVH::BoundingBoxTree boxtree;
       //BVH::BoundingBoxTree* head_node = &boxtree; 
        //BVH::generate_BVH_tree(scene ,head_node );
        bounding_boxes =  BVH::generate_bounding_boxes( scene  );

        if(scene.spheredir_lights.size() > 0 )
        {
            sphere_env_open_hdr( scene.images[scene.spheredir_lights[0].image_id].path );
        }
        for (size_t y = 0; y < height; y++)
        {
            for (size_t x = 0; x < width; x++)
            {
                parser::Vec3f  color;
                color.x = 0.0f;
                color.y = 0.0f;
                color.z = 0.0f;
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
                
                if( color.x == -1.0f ) //replace background 
                {
                   color.x = scene.background.image->image[i];
                   color.y =  scene.background.image->image[i + 1] ;
                   color.z =  scene.background.image->image[i + 2] ;
                }
                else //general case 
                {
                    color = color / num_samples ;
                    if( current_camera.tonemap.tmo == "")
                    {
                        color.x = std::min( 249.0f , color.x );
                        color.y = std::min( 249.0f , color.y );
                        color.z = std::min( 249.0f , color.z );
                    }
                }
                if( current_camera.tonemap.tmo == "")
                {
                    image[i ++ ] = (unsigned char) (color.x);
                    image[i ++ ] = (unsigned char) (color.y);
                    image[i ++ ] = (unsigned char) (color.z);
                }
                else
                {

                    tonemap_image[i++] = color.x;
                    tonemap_image[i++] = color.y;
                    tonemap_image[i++] = color.z;

                }
               /*  //initilaize recursion depth
                max_recursion_depth = scene.max_recursion_depth;
                // initialize the ray
                parser::Vec3f current_pixel_world_space =  (starting_p114oint_parser +  ( parser::Vec3f( interval_row.x ,  interval_row.y ,  interval_row.z) * (float)x )  )  +   ( parser::Vec3f( interval_col.x ,  interval_col.y ,  interval_col.z) * (float )y  )  ; 
                Ray ray(parser::Vec3f( current_camera.position.x , current_camera.position.y , current_camera.position.z  )  , current_pixel_world_space );
                parser::Vec3f  color = color_pixel(scene  , ray);
                // color cast 
                image[i++] = (unsigned char) (cstd::vector< BVH::BoundingBox >olor.x);
                image[i++] = (unsigned char) (color.y);
                image[i++] = (unsigned char) (color.z);
               */
                // color cast 
            
            }
        }
        //BVH::delete_boxes(head_node);
        //save the written image
        //stbi_write_png(current_camera.image_name.c_str() , width * num_samples_sqrt, height * num_samples_sqrt,3 , image , width * num_samples_sqrt *  3  );
        if( current_camera.tonemap.tmo == "")
        {
            stbi_write_png(current_camera.image_name.c_str() , width , height,3 , image , width  *  3  );
        }
        else //tonemap
        {
            // 1 - calculate average world luminance
            int i = 0; 
            double average_world_lum = average_world_luminance(tonemap_image , width , height , luminances);
            for (size_t y = 0; y < height; y++)
            {
                for (size_t x = 0; x < width; x++)
                {
                    parser::Vec3f color( tonemap_image[i ] , tonemap_image[i + 1] , tonemap_image[i + 3 ]);
                    parser::Vec3f tonemapped_color = apply_tonemap(  color ,current_camera  ,  average_world_lum  , luminances );
                    image[i++] = (unsigned char)tonemapped_color.x;
                    image[i++] = (unsigned char)tonemapped_color.y;
                    image[i++] = (unsigned char)tonemapped_color.z;
                }
            }

            std::string png = current_camera.image_name; 
            png = png.substr(0 , png.size() - 4 );
            png = png + ".png";
            stbi_write_png( png.c_str() , width , height,3 , image , width  *  3  );
            
            const char** temp_err; 
            SaveEXR(tonemap_image , width , height , 3 , 0 , current_camera.image_name.c_str()  , temp_err  );

            
        }
        
        //in the end just reallocate the image
        if(scene.spheredir_lights.size() > 0 )
        {
            release_sphere_env();
        }
        delete[] image; 
        delete[] tonemap_image; 
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

    global_interval_row = interval_row;
    global_interval_col = interval_col;

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
    bool is_light_object_intersected = false; 
    int object_id = -1; 
    bool  is_object_hit = ray_object_intersection( ray , scene ,  hit_point , normal  , material   , object_id , hit_face, is_shadow_rays_active , is_light_object_intersected);
    if( !is_object_hit)
    {
        if( scene.spheredir_lights.size() > 0 ) // if env spherical 
        {
            return spherical_background(ray.direction);
        }   
        else if(!scene.is_texture_background)
        {
            return parser::Vec3f(scene.background_color.x, scene.background_color.y , scene.background_color.z );
        }
        else
        {
            //return parse_background( current_camera , ray , global_interval_row , global_interval_col    );
            return parser::Vec3f(-1.0f ,0.0f ,0.0f );
        }
    }
    if( is_light_object_intersected && object_id >= scene.meshes.size() + scene.spheres.size() + scene.triangles.size() + scene.mesh_instances.size() )
    {
        if( object_id <  scene.meshes.size() + scene.spheres.size() + scene.triangles.size() + scene.mesh_instances.size() + scene.light_meshes.size() )
        {
            //light mesh get object
            parser::LightMesh* light_mesh_ptr =  &scene.light_meshes[object_id - ( scene.meshes.size() + scene.spheres.size() + scene.triangles.size() + scene.mesh_instances.size()) ];
            return light_mesh_ptr->radiance / ( parser::distance( ray.origin , hit_point));
        }
        else if( object_id <  scene.meshes.size() + scene.spheres.size() + scene.triangles.size() + scene.mesh_instances.size() + scene.light_meshes.size() + scene.light_spheres.size() )
        {
            parser::LightSphere* light_sphere_ptr =  &scene.light_spheres[object_id - ( scene.meshes.size() + scene.spheres.size() + scene.triangles.size() + scene.mesh_instances.size() + scene.light_meshes.size()) ];
            return light_sphere_ptr->radiance / ( parser::distance( ray.origin , hit_point));
        }
    }
    //check texture
    parser::Vec3f texture_color; 
    parser::Material texture_material; 
    if( current_camera.is_importance_sampling ) 
    {
        //sample hemisphere
        float eps1 = dis_importance_sampling_eps1(gen_importance_eps1);
        float eps2 = dis_importance_sampling_eps2(gen_importance_eps2);

        float phi =  2 * M_PI * eps1;
        float theta =  std::asin(std::sqrt(eps1));
        
        //1 -  rotate phi around + y axis
        //rotate phi degrees around 0
        parser::Vec4f unit_vec;
        unit_vec.x = 0.0f;
        unit_vec.y = 0.0f;
        unit_vec.z = -1.0f;
        unit_vec.w = 1;

        float Rx = 0.0f;
        float Ry = 1.0f;
        float Rz = 0.0f;

        parser::Matrix s;
        //rotation matrix 
        //  x is degree 
        // y z w arbitrary axis
        float cosine_theta = std::cos( phi  );
        float sine_theta = std::sin(   phi );
        s.set(0, 0 , cosine_theta + Rx * Rx * (1 - cosine_theta )  );
        s.set(0, 1 , Rx * Ry * (1 - cosine_theta ) - Rz * sine_theta  );  
        s.set(0, 2 , Rx * Rz * (1 -cosine_theta ) + Ry *sine_theta );  
        s.set(0, 3 , 0 );

        s.set(1, 0 , Ry*Rx*(1-cosine_theta) + Rz*sine_theta);
        s.set(1, 1 , cosine_theta + Ry*Ry*(1-cosine_theta) );
        s.set(1, 2 , Ry*Rz*(1-cosine_theta)-Rx*sine_theta );
        s.set(1, 3 , 0 );

        s.set(2, 0 , Ry*Rx*(1-cosine_theta) - Ry*sine_theta);
        s.set(2, 1 ,Rz*Ry*(1-cosine_theta) + Rx*sine_theta );
        s.set(2, 2 , cosine_theta + Rz*Rz*(1-cosine_theta) );
        s.set(2, 3 , 0 );

        s.set(3,0,0);
        s.set(3,1,0);
        s.set(3,2,0);
        s.set(3,3,1);  
        
        unit_vec = s * unit_vec;

        parser::Matrix s2;
        // theta around + x 
        Rx = 1.0f;
        Ry = 0.0f;
        Rz = 0.0f;

        cosine_theta = std::cos( theta );
        sine_theta = std::sin( theta );
        s2.set(0, 0 , cosine_theta + Rx * Rx * (1 - cosine_theta )  );
        s2.set(0, 1 , Rx * Ry * (1 - cosine_theta ) - Rz * sine_theta  );  
        s2.set(0, 2 , Rx * Rz * (1 -cosine_theta ) + Ry *sine_theta );  
        s2.set(0, 3 , 0 );

        s2.set(1, 0 , Ry*Rx*(1-cosine_theta) + Rz*sine_theta);
        s2.set(1, 1 , cosine_theta + Ry*Ry*(1-cosine_theta) );
        s2.set(1, 2 , Ry*Rz*(1-cosine_theta)-Rx*sine_theta );
        s2.set(1, 3 , 0 );

        s2.set(2, 0 , Ry*Rx*(1-cosine_theta) - Ry*sine_theta);
        s2.set(2, 1 ,Rz*Ry*(1-cosine_theta) + Rx*sine_theta );
        s2.set(2, 2 , cosine_theta + Rz*Rz*(1-cosine_theta) );
        s2.set(2, 3 , 0 );

        s2.set(3,0,0);
        s2.set(3,1,0);
        s2.set(3,2,0);
        s2.set(3,3,1);  

        unit_vec = s2 * unit_vec;
        
        unit_vec.x =   (unit_vec.x/unit_vec.w);
        unit_vec.y =   (unit_vec.y/unit_vec.w);
        unit_vec.z =   (unit_vec.z/unit_vec.w);


        Ray importance_ray( hit_point + parser::Vec3f(unit_vec.x , unit_vec.y ,unit_vec.z ) * 0.01f, hit_point + parser::Vec3f(unit_vec.x , unit_vec.y ,unit_vec.z ) );
        
        //toggle importance sampling in order to not have stack overflow 
        current_camera.is_importance_sampling = false; 
        color = color + color_pixel(scene , importance_ray);
        current_camera.is_importance_sampling = true; 

    }
    if( current_camera.is_next_event_estimation)
    {

    }
    
    bool is_intersection_textured = is_texture_present(  scene ,   object_id ,  hit_point  ,  normal , hit_face ,  texture_color ,  texture_material);
    //we can add ambient shading
    color =  color + parser::Vec3f(scene.ambient_light.x , scene.ambient_light.y , scene.ambient_light.z);
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
        bool is_light_object_intersected  = false; 
        bool  is_object_hit = ray_object_intersection( shadow_ray , scene ,  hit_point_temp , normal_temp , material_temp , object_id , hit_face_temp ,  is_shadow_rays_active , is_light_object_intersected);
        
        if( is_object_hit ) // it is in shadow no contribution from light
        {
            if( parser::distance(hit_point_temp , hit_point) < parser::distance(hit_point, light_pos) ) //before light source 
            {
                continue; 
            }
        }
        
        if( material.is_BRDF == false  ) // no brdf 
        {
                // diffuse 
            parser::Vec3f diffuse = parser::Vec3f( material.diffuse.x , material.diffuse.y, material.diffuse.z);
            float cosine_alpha = parser::dot(shadow_ray.direction  ,  normal ) / (parser::length(shadow_ray.direction) * parser::length(normal)  ) ; //normal vectors
            
            if( is_intersection_textured )
            {
                diffuse = parser::Vec3f( texture_material.diffuse.x * scene.point_lights[i].intensity.x , texture_material.diffuse.y *  scene.point_lights[i].intensity.y , texture_material.diffuse.z *   scene.point_lights[i].intensity.z) / (parser::distance(light_pos , hit_point )*parser::distance(light_pos, hit_point)  );
                //clamp cosine alpha
                cosine_alpha = std::max(0.0f , cosine_alpha);
                diffuse = diffuse * cosine_alpha; 

            }
            else
            {
                diffuse = parser::Vec3f( diffuse.x * scene.point_lights[i].intensity.x , diffuse.y *  scene.point_lights[i].intensity.y ,diffuse.z *   scene.point_lights[i].intensity.z) / (parser::distance(light_pos , hit_point )*parser::distance(light_pos, hit_point)  );
                //clamp cosine alpha
                cosine_alpha = std::max(0.0f , cosine_alpha);
                diffuse = diffuse * cosine_alpha; 
            }
            if( current_camera.tonemap.tmo == "") //empty 
            {//clamp diffuse
                diffuse.x = std::min(255.0f , diffuse.x );
                diffuse.y = std::min(255.0f , diffuse.y );
                diffuse.z = std::min(255.0f , diffuse.z );
            }


            // specular 
            // specular 
            parser::Vec3f specular = parser::Vec3f( material.specular.x , material.specular.y, material.specular.z);
            float phong_exponent = material.phong_exponent; 

            parser::Vec3f viewpos = parser::normalize(ray.origin - hit_point); 
            parser::Vec3f h =  parser::normalize( shadow_ray.direction + viewpos ); // this might be flawed 
            float cosine_h_n = parser::dot(h , normal) /  ( parser::length(h) * parser::length(normal) ); // they re both normalised but just in case
            //clamp cosine
            cosine_h_n = std::max(0.0f , cosine_h_n );
            if( is_intersection_textured )
            {
                specular =    parser::Vec3f(  texture_material.specular.x * scene.point_lights[i].intensity.x , texture_material.specular.y * scene.point_lights[i].intensity.y  , texture_material.specular.z * scene.point_lights[i].intensity.z ) / (parser::distance(light_pos , hit_point )*parser::distance(light_pos, hit_point)  ) * std::pow( cosine_h_n , phong_exponent   );  
            }
            else
            {
                specular =    parser::Vec3f( specular.x * scene.point_lights[i].intensity.x , specular.y * scene.point_lights[i].intensity.y  , specular.z * scene.point_lights[i].intensity.z ) / (parser::distance(light_pos , hit_point )*parser::distance(light_pos, hit_point)  ) * std::pow( cosine_h_n , phong_exponent   );  
            }
            if( current_camera.tonemap.tmo == "") //empty 
            {
                //clamp diffuse
                specular.x = std::min(255.0f , specular.x );
                specular.y = std::min(255.0f , specular.y );
                specular.z = std::min(255.0f , specular.z );
            }
            // replace_all u ekle 
            if( texture_material.diffuse.x == -1.0f && texture_material.specular.x == -1.0f )
            {
                //replace all 
                color = color + texture_color; // sphere_point_hdr icin
            }
            else
            {
                //normal case 
                color =  color + diffuse  + specular  ;   
            }
        }
        else // material.isbrdf
        {
            parser::Material material_temp; 
            parser::Vec3f diffuse;
            parser::Vec3f specular;

            
            diffuse = parser::Vec3f( material.diffuse.x , material.diffuse.y, material.diffuse.z);
            float cosine_theta_value = parser::dot( normal , shadow_ray.direction ) / ( parser::length(shadow_ray.direction) * parser::length(normal) ) ; 
            cosine_theta_value = std::max(0.0f , cosine_theta_value);

            parser::Vec3f W0 = parser::Vec3f( shadow_ray.direction.x  * -1 , shadow_ray.direction.y  * -1 , shadow_ray.direction.z  * -1  ); 
            float cosine_omega  = parser::dot( W0 , normal ) / ( parser::length( normal ) * parser::length(W0));
            parser::Vec3f Wr = W0 * -1.0f +  normal * cosine_omega * 2.0f ;  // perfect directed vector 

            float cosine_alpha_value = parser::dot(Wr , ray.direction * -1) / (parser::length(Wr) * parser::length(ray.direction));

            parser::Vec3f viewpos = parser::normalize(ray.origin - hit_point); 
            parser::Vec3f h =  parser::normalize( shadow_ray.direction + viewpos ); // this might be flawed 
            float cosine_h_n = parser::dot(h , normal) /  ( parser::length(h) * parser::length(normal) ); // they re both normalised but just in case
            //clamp cosine
            cosine_h_n = std::max(0.0f , cosine_h_n );


            if( is_intersection_textured )
            {
                diffuse = parser::Vec3f( texture_material.diffuse.x , texture_material.diffuse.y, texture_material.diffuse.z);
                specular = parser::Vec3f( texture_material.specular.x , texture_material.specular.y, texture_material.specular.z);
            }
            else
            {
                diffuse = parser::Vec3f( material.diffuse.x , material.diffuse.y, material.diffuse.z);
                specular = parser::Vec3f( material.specular.x , material.specular.y, material.specular.z);
            }
            material_temp = material; 
            material_temp.diffuse = diffuse; 
            material_temp.specular = specular; 
            parser::Vec3f brdf_material  = get_BRDF(scene , material_temp ,  cosine_alpha_value , cosine_theta_value,  cosine_h_n , current_camera  , shadow_ray.direction , ray.direction * -1 , h , normal );
            
            parser::Vec3f new_color(brdf_material.x * scene.point_lights[i].intensity.x , brdf_material.y *  scene.point_lights[i].intensity.y ,brdf_material.z *   scene.point_lights[i].intensity.z);
            new_color = new_color / (parser::distance(light_pos , hit_point )*parser::distance(light_pos, hit_point)  );
            new_color = new_color * cosine_theta_value;
            
            
            if(current_camera.is_russian_roulette)
            {
                float eps1 = dis_russian_roulette_eps(gen_russian_roulette_eps);
                // use average of the 3 brdf channels ? 
                float brdf = (brdf_material.x + brdf_material.y + brdf_material.z) / 3.0f ;

                color = color +  new_color * ray.throughput ;
                ray.throughput = ray.throughput * brdf;

                if( !(eps1 >=  ray.throughput) ) //continue recursion
                {
                    color = color + color_pixel(scene , ray);
                }
                else
                {
                    color = color +  new_color / ray.throughput ;
                }
            }
            else
            {
                color = color + new_color; //plain addition
                if( new_color.x > 100 )
                int a = 1; 
            }
        }
        

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
        
        bool is_light_object_intersected  = false; 
        bool  is_object_hit = ray_object_intersection( shadow_ray , scene ,  hit_point_temp , normal_temp , material_temp , object_id , hit_face_temp , is_shadow_rays_active , is_light_object_intersected);
        
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
        //clamp cosine135
        cosine_h_n = std::max(0.0f , cosine_h_n );

        specular =    parser::Vec3f( specular.x * L.x, specular.y * L.y  , specular.z * L.z )  * std::pow( cosine_h_n , phong_exponent   );  
        //clamp diffuse
        specular.x = std::min(255.0f , specular.x );
        specular.y = std::min(255.0f , specular.y );
        specular.z = std::min(255.0f , specular.z );
        
         
        color =  color + diffuse  + specular  ; 

    }
    //  for each directional light 
    is_shadow_rays_active = true; 
    for (size_t i = 0; i < scene.directional_lights.size(); i++)
    {
         // cast  shadow ray
        Ray shadow_ray(hit_point , hit_point - scene.directional_lights[i].direction );
        shadow_ray.time_parameter = ray.time_parameter;
        
        // else you directly went to a light
        
        parser::Vec3f hit_point_temp; 
        parser::Face hit_face_temp; 
        parser::Vec3f normal_temp;
        parser::Material material_temp;  
        bool is_light_object_intersected  = false; 

        bool  is_object_hit = ray_object_intersection( shadow_ray , scene ,  hit_point_temp , normal_temp , material_temp , object_id , hit_face_temp ,  is_shadow_rays_active , is_light_object_intersected );
        
        if( is_object_hit ) // it is in shadow no contribution from light
        {
                continue; 
        }
        
        // diffuse 
        parser::Vec3f diffuse = parser::Vec3f( material.diffuse.x , material.diffuse.y, material.diffuse.z);
        float cosine_alpha = parser::dot(shadow_ray.direction  ,  normal ) / (parser::length(shadow_ray.direction) * parser::length(normal)  ) ; //normal vectors
        
        if( is_intersection_textured )
        {
            diffuse = parser::Vec3f( texture_material.diffuse.x * scene.directional_lights[i].radiance.x , texture_material.diffuse.y *  scene.directional_lights[i].radiance.y , texture_material.diffuse.z *   scene.directional_lights[i].radiance.z) ;
            //clamp cosine alpha
            cosine_alpha = std::max(0.0f , cosine_alpha);
            diffuse = diffuse * cosine_alpha; 

        }
        else
        {
            diffuse = parser::Vec3f( diffuse.x * scene.directional_lights[i].radiance.x , diffuse.y *  scene.directional_lights[i].radiance.y ,diffuse.z *   scene.directional_lights[i].radiance.z) ;
             //clamp cosine alpha
            cosine_alpha = std::max(0.0f , cosine_alpha);
            diffuse = diffuse * cosine_alpha; 
        }
        if( current_camera.tonemap.tmo == "") //empty 
            {//clamp diffuse
            diffuse.x = std::min(255.0f , diffuse.x );
            diffuse.y = std::min(255.0f , diffuse.y );
            diffuse.z = std::min(255.0f , diffuse.z );
        }


        // specular 
        // specular 
        parser::Vec3f specular = parser::Vec3f( material.specular.x , material.specular.y, material.specular.z);
        float phong_exponent = material.phong_exponent; 

        parser::Vec3f viewpos = parser::normalize(ray.origin - hit_point); 
        parser::Vec3f h =  parser::normalize( shadow_ray.direction + viewpos ); // this might be flawed 
        float cosine_h_n = parser::dot(h , normal) /  ( parser::length(h) * parser::length(normal) ); // they re both normalised but just in case
        //clamp cosine
        cosine_h_n = std::max(0.0f , cosine_h_n );
        if( is_intersection_textured )
        {
            specular =    parser::Vec3f(  texture_material.specular.x * scene.directional_lights[i].radiance.x , texture_material.specular.y * scene.directional_lights[i].radiance.y  , texture_material.specular.z * scene.directional_lights[i].radiance.z ) * std::pow( cosine_h_n , phong_exponent   );  
        }
        else
        {
            specular =    parser::Vec3f( specular.x * scene.directional_lights[i].radiance.x , specular.y * scene.directional_lights[i].radiance.y  , specular.z * scene.directional_lights[i].radiance.z )  * std::pow( cosine_h_n , phong_exponent   );  
        }
        if( current_camera.tonemap.tmo == "") //empty 
        {
            //clamp specular
            specular.x = std::min(255.0f , specular.x );
            specular.y = std::min(255.0f , specular.y );
            specular.z = std::min(255.0f , specular.z );
        }
         
        color =  color + diffuse  + specular  ;   
    }
    is_shadow_rays_active = true; 
    for (size_t i = 0; i < scene.spot_lights.size(); i++)
    {
         // cast  shadow ray
        Ray shadow_ray(hit_point , scene.spot_lights[i].position );
        shadow_ray.time_parameter = ray.time_parameter;
        
        // else you directly went to a light
        
        parser::Vec3f hit_point_temp; 
        parser::Face hit_face_temp; 
        parser::Vec3f normal_temp;
        parser::Material material_temp;  
        bool is_light_object_intersected  = false; 
        bool  is_object_hit = ray_object_intersection( shadow_ray , scene ,  hit_point_temp , normal_temp , material_temp , object_id , hit_face_temp ,  is_shadow_rays_active , is_light_object_intersected);
        
        if( is_object_hit ) // it is in shadow no contribution from light
        {
            if( parser::distance(hit_point_temp , hit_point) < parser::distance(hit_point, scene.spot_lights[i].position) ) //before light source 
            {
                continue; 
            }
        }
        // calculate the angle between hitpoint and spot
        parser::Vec3f spot_to_hit = hit_point - scene.spot_lights[i].position;
        float alpha_for_spot = parser::dot(spot_to_hit ,scene.spot_lights[i].direction ) / ( parser::length(spot_to_hit ) * parser::length(scene.spot_lights[i].direction ));
        //now calculate if L1 L2 or L3
        float alpha_for_spot_value = std::acos(alpha_for_spot) * 180/M_PI;
        parser::Vec3f irradiance;
        if(  alpha_for_spot_value <= scene.spot_lights[i].fallofAngle/2  && alpha_for_spot_value >= -scene.spot_lights[i].fallofAngle/2   ) //L1 
        {
            //L1
            irradiance = scene.spot_lights[i].intensity / parser::length(spot_to_hit);
        }
        else if( alpha_for_spot_value <= scene.spot_lights[i].coverageAngle/2 && alpha_for_spot_value >= scene.spot_lights[i].coverageAngle/2 ) 
        {
            float cos_fallof_2 = std::cos( M_PI/180 * scene.spot_lights[i].fallofAngle/2  );
            float cos_cover_2 = std::cos( M_PI/180 * scene.spot_lights[i].coverageAngle/2  );

            float s = std::pow( (alpha_for_spot - cos_fallof_2 ) / (cos_fallof_2 - cos_cover_2)    , 4 );
            //L2
            irradiance =  scene.spot_lights[i].intensity / parser::length(spot_to_hit) * s;
        }
        else 
        {
            //L3
            irradiance.x = 0.0f;
            irradiance.y = 0.0f;
            irradiance.z = 0.0f;

        }

        // diffuse 
        parser::Vec3f diffuse = parser::Vec3f( material.diffuse.x , material.diffuse.y, material.diffuse.z);
        float cosine_alpha = parser::dot(shadow_ray.direction  ,  normal ) / (parser::length(shadow_ray.direction) * parser::length(normal)  ) ; //normal vectors
        
        if( is_intersection_textured )
        {
            diffuse = parser::Vec3f( texture_material.diffuse.x * irradiance.x , texture_material.diffuse.y *  irradiance.y , texture_material.diffuse.z *   irradiance.z) ;
            //clamp cosine alpha
            cosine_alpha = std::max(0.0f , cosine_alpha);
            diffuse = diffuse * cosine_alpha; 

        }
        else
        {
            diffuse = parser::Vec3f( diffuse.x * irradiance.x , diffuse.y *  irradiance.y ,diffuse.z *   irradiance.z) ;
             //clamp cosine alpha
            cosine_alpha = std::max(0.0f , cosine_alpha);
            diffuse = diffuse * cosine_alpha; 
        }
        if( current_camera.tonemap.tmo == "") //empty 
            {//clamp diffuse
            diffuse.x = std::min(255.0f , diffuse.x );
            diffuse.y = std::min(255.0f , diffuse.y );
            diffuse.z = std::min(255.0f , diffuse.z );
        }


        // specular 
        parser::Vec3f specular = parser::Vec3f( material.specular.x , material.specular.y, material.specular.z);
        float phong_exponent = material.phong_exponent; 

        parser::Vec3f viewpos = parser::normalize(ray.origin - hit_point); 
        parser::Vec3f h =  parser::normalize( shadow_ray.direction + viewpos ); // this might be flawed 
        float cosine_h_n = parser::dot(h , normal) /  ( parser::length(h) * parser::length(normal) ); // they re both normalised but just in case
        //clamp cosine
        cosine_h_n = std::max(0.0f , cosine_h_n );
        if( is_intersection_textured )
        {
            specular =    parser::Vec3f(  texture_material.specular.x * irradiance.x , texture_material.specular.y * irradiance.y  , texture_material.specular.z * irradiance.z ) * std::pow( cosine_h_n , phong_exponent   );  
        }
        else
        {
            specular =    parser::Vec3f( specular.x * irradiance.x , specular.y * irradiance.y  , specular.z * irradiance.z )  * std::pow( cosine_h_n , phong_exponent   );  
        }
        if( current_camera.tonemap.tmo == "") //empty 
        {
            //clamp specular
            specular.x = std::min(255.0f , specular.x );
            specular.y = std::min(255.0f , specular.y );
            specular.z = std::min(255.0f , specular.z );
        }
         
        color =  color + diffuse  + specular  ;   
    }  
    is_shadow_rays_active = true; 
    for (size_t i = 0; i < scene.spheredir_lights.size(); i++)
    {
        //random rejection sampling
        parser::Vec3f random_vec = get_random_vector_rejection_sampling( normal );
        color = color + get_hdr_image_color(random_vec);
    }
    is_shadow_rays_active = true; 
    for (size_t i = 0; i < scene.light_meshes.size(); i++)
    {
        if( current_camera.is_next_event_estimation )
        {
            //do sampling I dont know how
        }
    }
    is_shadow_rays_active = true; 
    for (size_t i = 0; i < scene.light_spheres.size(); i++)
    {
        if( current_camera.is_next_event_estimation )
        {
            //do sampling I dont know how
        }
    }
    
    if(material.is_mirror ) // do recursion
    {
        if( max_recursion_depth > 0 )
        {
            
            max_recursion_depth -= 1; // max recursion depth is reduced 
            if( current_camera.is_russian_roulette )
            {
                max_recursion_depth += 1; // recursion depth becomes useless 
            }
            
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
            if( scene.spheredir_lights.size()  > 0 )
            {
            color =  parser::Vec3f( recursion_color.x  * material.mirror.x ,  recursion_color.y  * material.mirror.y  ,  recursion_color.z  * material.mirror.z ); 
            }
            else
            {
            color =  color + parser::Vec3f( recursion_color.x  * material.mirror.x ,  recursion_color.y  * material.mirror.y  ,  recursion_color.z  * material.mirror.z ); 
            }

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
            if( current_camera.is_russian_roulette )
            {
                max_recursion_depth += 1; // recursion depth becomes useless 
            }
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
        
        if(scene.spheredir_lights.size() > 0 )
        {
            color =  refraction_color + reflection_color;
        }
        else
        {
            color = color +  refraction_color + reflection_color;

        }

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
            if( current_camera.is_russian_roulette )
            {
                max_recursion_depth += 1; // recursion depth becomes useless 
            }
            Ray reflection_ray( hit_point +  ( vr * (0.01f) )  , hit_point + vr  );
            reflection_ray.generate_time();
            parser::Vec3f reflection_color =   color_pixel(scene ,reflection_ray) * reflection_ratio ;
            reflection_color = parser::Vec3f( reflection_color.x * material.mirror.x , reflection_color.y * material.mirror.y , reflection_color.z * material.mirror.z);
            color = color + reflection_color; 
        }
       

    }
    //clamp color
    if (current_camera.tonemap.tmo == "")
    {
        color.x = std::min( color.x , 255.0f );
        color.y = std::min( color.y , 255.0f );
        color.z = std::min( color.z , 255.0f );

            //clamp to nearest integer
        color.x = (float) ( (int ) color.x );
        color.y = (float) ( (int ) color.y );
        color.z = (float) ( (int ) color.z );
    }

    color.x = std::max(color.x , 0.0f );
    color.y = std::max(color.y , 0.0f );
    color.z = std::max(color.z , 0.0f );

   

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