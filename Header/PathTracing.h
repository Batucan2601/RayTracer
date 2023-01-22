#pragma once
#include "../System_Files/support_files/parser.h"
#include "../Header/Ray.h"
#include "../Header/Prototypes.h"
#include <random>
// generator for importance sampling

// epsilon1 
std::random_device rd_importance_sampling_eps1;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_importance_eps1(rd_importance_sampling_eps1()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_importance_sampling_eps1(0 ,1);

//epslion2 
std::random_device rd_importance_sampling_eps2;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_importance_eps2(rd_importance_sampling_eps2()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_importance_sampling_eps2(0 ,1);

//generator for russian roulette 
std::random_device rd_russian_roulette_eps;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_russian_roulette_eps(rd_russian_roulette_eps()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_russian_roulette_eps(0 ,1);

//generaor for next event estimation 
// epsilon1 
std::random_device rd_next_event_eps1;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_next_event_eps1(rd_next_event_eps1()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_next_event_eps1(0 ,1);

//epslion2 
std::random_device rd_next_event_eps2;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_next_event_eps2(rd_next_event_eps2()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_next_event_eps2(0 ,1);

// epsilon for random triangle
std::random_device rd_random_triangle_eps;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_random_triangle_eps(rd_random_triangle_eps()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_random_triangle_eps(0 ,1);

// 2 epsilons for random points inside a random triangle
std::random_device rd_random_triangle_point_eps1;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_random_triangle_point_eps1(rd_random_triangle_point_eps1()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_random_triangle_point_eps1(0 ,1);

std::random_device rd_random_triangle_point_eps2;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_random_triangle_point_eps2(rd_random_triangle_point_eps2()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_random_triangle_point_eps2(0 ,1);

std::random_device rd_next_event;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_next_event(rd_next_event()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_next_event(0 ,1);

std::vector<std::vector<float> > light_mesh_cdf; //prob for each triangle 
std::vector<float> light_mesh_areas; 

std::random_device rd_sphere_eps1 ;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_sphere_eps1(rd_sphere_eps1()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_sphere_eps1(0 ,1);

std::random_device rd_sphere_eps2 ;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_sphere_eps2(rd_sphere_eps2()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_sphere_eps2(0 ,1);

void calculate_total_area_light_mesh(parser::Scene &scene )
{
    if(scene.light_meshes.size()  == 0 )
    {
        return; 
    }

    for (size_t j = 0; j < scene.light_meshes.size(); j++)
    {
        parser::LightMesh light_mesh = scene.light_meshes[j];
        float total_area = 0;
        for (size_t i = 0; i < light_mesh.faces.size(); i++)
        {
            parser::Vec3f p1 = parser::Vec3f( scene.vertex_data[ (light_mesh.faces[i].v0_id-1)    ].x , scene.vertex_data[(light_mesh.faces[i].v0_id-1)    ].y , scene.vertex_data[ (light_mesh.faces[i].v0_id-1)  ].z); 
            parser::Vec3f p2 = parser::Vec3f( scene.vertex_data[ (light_mesh.faces[i].v1_id-1)   ].x , scene.vertex_data[(light_mesh.faces[i].v1_id-1)].y , scene.vertex_data[ (light_mesh.faces[i].v1_id-1)   ].z);
            parser::Vec3f p3 = parser::Vec3f( scene.vertex_data[ (light_mesh.faces[i].v2_id-1)   ].x , scene.vertex_data[(light_mesh.faces[i].v2_id-1)].y , scene.vertex_data[ (light_mesh.faces[i].v2_id-1)   ].z);

            //area
            total_area +=  parser::length( parser::cross( p2 - p1 , p3  - p1 ) ) / 2;
           
        }
        light_mesh_areas.push_back(total_area);
        std::vector<float> cdf_for_single_light_mesh;
        float prob_so_far = 0;  
        //no we know the area we cna convert the triangles to cdf
        for (size_t i = 0; i < light_mesh.faces.size(); i++)
        {
            parser::Vec3f p1 = parser::Vec3f( scene.vertex_data[ (light_mesh.faces[i].v0_id-1)    ].x , scene.vertex_data[(light_mesh.faces[i].v0_id-1)    ].y , scene.vertex_data[ (light_mesh.faces[i].v0_id-1)  ].z); 
            parser::Vec3f p2 = parser::Vec3f( scene.vertex_data[ (light_mesh.faces[i].v1_id-1)   ].x , scene.vertex_data[(light_mesh.faces[i].v1_id-1)].y , scene.vertex_data[ (light_mesh.faces[i].v1_id-1)   ].z);
            parser::Vec3f p3 = parser::Vec3f( scene.vertex_data[ (light_mesh.faces[i].v2_id-1)   ].x , scene.vertex_data[(light_mesh.faces[i].v2_id-1)].y , scene.vertex_data[ (light_mesh.faces[i].v2_id-1)   ].z);

            //probability =  triangle area / total area 
            float prob = parser::length( parser::cross( p2 - p1 , p3  - p1 ) ) / 2  / total_area;
            prob_so_far += prob; 
            cdf_for_single_light_mesh.push_back(prob_so_far);
        }
        light_mesh_cdf.push_back(cdf_for_single_light_mesh);
        
    }
    
       
}
parser::Vec3f calculate_random_point_inside_triangle(  const parser::Scene & scene , const parser::Face &face )
{
    parser::Vec3f p1 = parser::Vec3f( scene.vertex_data[ (face.v0_id-1)    ].x , scene.vertex_data[(face.v0_id-1)    ].y , scene.vertex_data[(face.v0_id-1)  ].z); 
    parser::Vec3f p2 = parser::Vec3f( scene.vertex_data[ (face.v1_id-1)   ].x , scene.vertex_data[(face.v1_id-1)].y , scene.vertex_data[ (face.v1_id-1)   ].z);
    parser::Vec3f p3 = parser::Vec3f( scene.vertex_data[ (face.v2_id-1)   ].x , scene.vertex_data[(face.v2_id-1)].y , scene.vertex_data[ (face.v2_id-1)   ].z);

    //random number between p1 and p2
    float eps1 = dis_random_triangle_point_eps1(gen_random_triangle_point_eps1);
    float eps2 = dis_random_triangle_point_eps2(gen_random_triangle_point_eps2);

    parser::Vec3f p;
    p = p1 * eps1 +  p2 * (1 - eps1);

    parser::Vec3f q; 
    q = p3 * std::sqrt(eps2 ) + p * ( 1 - std::sqrt(1 - eps2 ));


    return q; 

}

parser::Vec3f calculate_random_point_on_sphere(   parser::Scene & scene  , parser::Vec3f &hit_point ,  parser::LightSphere & sphere , float & cosine_theta_max_input )
{
    parser::Vec3f center = parser::Vec3f(scene.vertex_data[sphere.center_vertex_id-1].x , scene.vertex_data[sphere.center_vertex_id-1 ].y  , scene.vertex_data[sphere.center_vertex_id-1].z );
    
    float theta_max = std::asin(  sphere.radius / parser::distance(  center , hit_point) );
    
    float cosine_theta_max = std::sqrt( 1 -  std::sin(theta_max * theta_max));
    cosine_theta_max_input = cosine_theta_max;

    float x = dis_sphere_eps1(gen_sphere_eps1);
    float omega = std::acos(1 - x  + x * cosine_theta_max );
    float phi = 2 * M_PI * dis_sphere_eps2(gen_sphere_eps2);

    parser::Vec3f w = center - hit_point;
    parser::Vec3f u(hit_point.x+ 1.0f, hit_point.y ,hit_point.z);

    parser::Vec3f v = parser::cross( u , w );

    parser::Vec3f l = w * std::cos( omega )  + v *  std::sin(omega) * std::cos(phi) +u * std::sin(omega) * std::sin(phi);
    
    Ray ray( hit_point , hit_point + l );


    parser::Sphere temp_sphere;
    temp_sphere.radius = sphere.radius;
    temp_sphere.center_vertex_id = sphere.center_vertex_id;
    temp_sphere.transformation.set(0,0,1);
    temp_sphere.transformation.set(0,1,0);
    temp_sphere.transformation.set(0,2,0);
    temp_sphere.transformation.set(0,3,0);
    temp_sphere.transformation.set(1,0,0);
    temp_sphere.transformation.set(1,1,1);
    temp_sphere.transformation.set(1,2,0);
    temp_sphere.transformation.set(1,3,0);
    temp_sphere.transformation.set(2,0,0);
    temp_sphere.transformation.set(2,1,0);
    temp_sphere.transformation.set(2,2,1);
    temp_sphere.transformation.set(2,3,0);
    temp_sphere.transformation.set(3,0,0);
    temp_sphere.transformation.set(3,1,0);
    temp_sphere.transformation.set(3,2,0);
    temp_sphere.transformation.set(3,3,1);

    temp_sphere.transformation_inverse.set(0,0,1);
    temp_sphere.transformation_inverse.set(0,1,0);
    temp_sphere.transformation_inverse.set(0,2,0);
    temp_sphere.transformation_inverse.set(0,3,0);
    temp_sphere.transformation_inverse.set(1,0,0);
    temp_sphere.transformation_inverse.set(1,1,1);
    temp_sphere.transformation_inverse.set(1,2,0);
    temp_sphere.transformation_inverse.set(1,3,0);
    temp_sphere.transformation_inverse.set(2,0,0);
    temp_sphere.transformation_inverse.set(2,1,0);
    temp_sphere.transformation_inverse.set(2,2,1);
    temp_sphere.transformation_inverse.set(2,3,0);
    temp_sphere.transformation_inverse.set(3,0,0);
    temp_sphere.transformation_inverse.set(3,1,0);
    temp_sphere.transformation_inverse.set(3,2,0);
    temp_sphere.transformation_inverse.set(3,3,1);

    parser::Vec3f temp_normal;
    parser::Vec3f temp_hitpoint;
    bool success = false; 
    #pragma region ray sphere
     // multiply ray direction and origin with inverse of this
    parser::Vec4f origin_inv; 
    origin_inv.x = ray.origin.x;
    origin_inv.y = ray.origin.y;
    origin_inv.z = ray.origin.z;
    origin_inv.w = 1.0f;
    
    parser::Vec4f dir_inv; 
    dir_inv.x = ray.origin.x + ray.direction.x;;
    dir_inv.y = ray.origin.y + ray.direction.y;;
    dir_inv.z = ray.origin.z + ray.direction.z;;
    dir_inv.w = 1.0f;

    origin_inv =  temp_sphere.transformation_inverse * origin_inv ;
    dir_inv = temp_sphere.transformation_inverse * dir_inv;
    parser::Vec3f origin_inv_3f; 
    parser::Vec3f dir_inv_3f; 

    origin_inv_3f.x = origin_inv.x / origin_inv.w;
    origin_inv_3f.y = origin_inv.y/ origin_inv.w;
    origin_inv_3f.z = origin_inv.z/ origin_inv.w;

    dir_inv_3f.x = dir_inv.x / dir_inv.w;
    dir_inv_3f.y = dir_inv.y/ dir_inv.w;
    dir_inv_3f.z = dir_inv.z/ dir_inv.w;

    Ray inv_ray(origin_inv_3f  , dir_inv_3f  );

    //center translation 

    // at^2 + bT + c = 0
    parser::Vec3f O = inv_ray.origin; 
    parser::Vec3f k = inv_ray.direction; 
    

    /*std::cout << "  " << std::endl;
    std::cout << "3f " <<  radius_vec_3f.x << " " << radius_vec_3f.y << " " << radius_vec_3f.z << std::endl;
    std::cout << "4f " <<  radius_vec_4f.x << " " << radius_vec_4f.y << " " << radius_vec_4f.z << " " << radius_vec_4f.w <<  std::endl;
    std::cout << "  " << std::endl;*/
    
    float a = parser::dot(k , k );
    float b = 2 * parser::dot(k , O-center );
    float c = parser::dot(O-center , O-center) - (sphere.radius  * sphere.radius);
    float discriminant =  b*b - 4 * a * c; 

    #ifdef DEBUG 
    if( isnan(discriminant))
    {
        std::cout << " starrt " << std::endl; 
        std::cout << "center " << center.x << " " << center.y << " " << center.z << std::endl;
        std::cout << "ray dir " << k.x << " " << k.y << " " << k.z  << std::endl;
        std::cout << "ray origin " << O.x << " " << O.y << " " << O.z  << std::endl;
        std::cout << "a " <<  a <<   " b "  << b  << " c " << c <<  std::endl; 
        success =  false;
    }
    #endif 
    if( discriminant < 0 )
    {
        success =  false; 
    }
    
    if( discriminant == 0 )
    {
        float root =  ( -b )  / (2*a);
        if( root < 0 )
        {
            success= false; 
        }
        temp_hitpoint = (parser::Vec3f)O +  (parser::Vec3f)k * root;

        //calculate normal 
        temp_normal = parser::normalize( temp_hitpoint - center );


        //multiply with the transform matrix
        parser::Vec4f hitpoint_4f;
        hitpoint_4f.x= temp_hitpoint.x; 
        hitpoint_4f.y= temp_hitpoint.y; 
        hitpoint_4f.z= temp_hitpoint.z; 
        hitpoint_4f.w= 1.0f; 

        parser::Vec4f normal_inv;
        normal_inv.x= temp_normal.x; 
        normal_inv.y= temp_normal.y; 
        normal_inv.z= temp_normal.z; 
        normal_inv.w= 1.0f; 

        hitpoint_4f = temp_sphere.transformation * hitpoint_4f;

        //deneme
        parser::Matrix trans_inverse_transpose = temp_sphere.transformation.inverse();
        trans_inverse_transpose.Transpose();


        normal_inv = trans_inverse_transpose * normal_inv;

        temp_hitpoint.x = hitpoint_4f.x / hitpoint_4f.w;
        temp_hitpoint.y = hitpoint_4f.y / hitpoint_4f.w;
        temp_hitpoint.z = hitpoint_4f.z / hitpoint_4f.w;

        temp_normal.x = normal_inv.x / normal_inv.w;
        temp_normal.y = normal_inv.y / normal_inv.w;
        temp_normal.z = normal_inv.z / normal_inv.w;
    std::cout << parser::dot( temp_normal , temp_hitpoint ) << std::endl; 

        success =  true;  
    }
    // else there are two roots

    float root_1 =  (-b-std::sqrt(discriminant) )  / (2*a);
    float root_2 = (-b+std::sqrt(discriminant) ) / (2*a);
    if( root_1 < 0 && root_2 < 0  )
    {
        success = false; 
    }
    //two hitpoints
    parser::Vec3f h1 = O + k *root_1;
    parser::Vec3f h2 = O + k *root_2;
    if( root_1 > 0 && root_2 > 0 )
    {
        
        //get the smaller length root * k;
        float d1 = parser::distance( O ,  h1 );
        float d2 = parser::distance( O ,  h2 );
        if( d1 > d2 )
        {
            temp_hitpoint = h2; 
        }
        else
        {
            temp_hitpoint = h1;
        }
    }
    else
    {
        if( root_1 > 0 )
        {
            temp_hitpoint = h1; 
        }
        else if( root_2 > 0 )
        {
            temp_hitpoint = h2; 
        }
    }
    
    if( success )
    {
            //calculate normal
        temp_normal = parser::normalize( temp_hitpoint - center );

        //multiply with the transform matrix
        parser::Vec4f hitpoint_4f;
        hitpoint_4f.x= temp_hitpoint.x; 
        hitpoint_4f.y= temp_hitpoint.y; 
        hitpoint_4f.z= temp_hitpoint.z; 
        hitpoint_4f.w= 1.0f; 

        parser::Vec4f normal_inv;
        normal_inv.x= temp_normal.x; 
        normal_inv.y= temp_normal.y; 
        normal_inv.z= temp_normal.z; 
        normal_inv.w= 1.0f; 

        hitpoint_4f = temp_sphere.transformation * hitpoint_4f;

        //deneme
        parser::Matrix trans_inverse_transpose = temp_sphere.transformation.inverse();
        trans_inverse_transpose.Transpose();
        
        normal_inv = trans_inverse_transpose * normal_inv;

        temp_hitpoint.x = hitpoint_4f.x / hitpoint_4f.w;
        temp_hitpoint.y = hitpoint_4f.y / hitpoint_4f.w;
        temp_hitpoint.z = hitpoint_4f.z / hitpoint_4f.w;

        temp_normal.x = normal_inv.x;
        temp_normal.y = normal_inv.y;
        temp_normal.z = normal_inv.z;

        temp_normal = parser::normalize(temp_normal);
        success = true; 

    }   
 

    #pragma endregion raysphere
    if(success )
    {
        return temp_hitpoint;
    }
    else
    {
        return parser::Vec3f(0.0f,0.0f,0.0f);
    }
}