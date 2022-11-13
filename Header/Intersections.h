#pragma once
#include "Renderer.h"
#include <math.h>       /* isnan, std */
#include "../System_Files/support_files/parser.h"
static bool ray_plane_intersection(const Ray& ray , const parser::Vec3f & p1 ,  const parser::Vec3f & normal , parser::Vec3f & hitpoint  )
{
    
    float n_d = parser::dot(parser::normalize(normal) , parser::normalize(ray.direction ) );
    if( n_d < 1e-8 && n_d > -1e-8  ) //parallel dot 0 
    {
        #ifdef DEBUG
        //std::cout << "nd " <<  n_d   << std::endl; 
        //std::cout <<  " normal " << normal.x << " " <<  normal.y << " " << normal.z << std::endl; 
        //std::cout <<  " ray ir " << ray.direction.x << " " <<  ray.direction.y << " " << ray.direction.z << std::endl; 
        
        #endif

        return false; 
    }
    parser::Vec3f ro_p = (parser::Vec3f)p1 - (parser::Vec3f)ray.origin;
    float t = parser::dot(ro_p , normal ) / n_d; 

    //plug t
    hitpoint = (parser::Vec3f)ray.origin + (parser::Vec3f)ray.direction * t; 
    #ifdef DEBUG

    //std::cout << hitpoint.x << " " << hitpoint.y << "  " << hitpoint.z << std::endl; 
    if( isnan(hitpoint.x) ||isnan(hitpoint.y) || isnan(hitpoint.z) )
    {
        std::cout << "ray origin " << std::endl;
        std::cout << ray.origin.x << " " << ray.origin.y << " "  << ray.origin.z << " " <<std::endl;  
        std::cout << " ray dire" << std::endl; 
        std::cout << ray.direction.x << " " << ray.direction.y << " "  << ray.direction.z << " " <<std::endl;  
        std::cout <<  " t " << std::endl; 
        std::cout << t <<   " " << n_d <<   std::endl;
        std::cout << normal.x << " " << normal.y << " "  << normal.z << " " <<std::endl;  

    } 
    

    #endif
    if( t < 0 ) // line intersection not a ray 
    {

        return false; 
    }

    return true; 

    
}
static bool is_point_in_triangle(const parser::Vec3f & p1 ,const parser::Vec3f & p2 , const parser::Vec3f & p3 , const parser::Vec3f & hitpoint)
{
    // assume same plane

    
    parser::Vec3f normal = parser::cross((parser::Vec3f)p2 - (parser::Vec3f)p1 , (parser::Vec3f)p3 - (parser::Vec3f)p1 );
    
    parser::Vec3f C1 = parser::cross( (parser::Vec3f)p2 - (parser::Vec3f)p1 , (parser::Vec3f)hitpoint - (parser::Vec3f)p1  );
    
    if( parser::dot( C1 , normal) < 0 )
    {
        return false;
    }
    parser::Vec3f C2 = parser::cross( (parser::Vec3f)p3 - (parser::Vec3f)p2 , (parser::Vec3f)hitpoint - (parser::Vec3f)p2  );
    if( parser::dot( C2 , normal) < 0 )
    {

        return false;
    }
    parser::Vec3f C3 = parser::cross( (parser::Vec3f)p1 - (parser::Vec3f)p3 , (parser::Vec3f)hitpoint - (parser::Vec3f)p3  );
    
    if( parser::dot( C3 , normal) < 0 )
    {

        return false;
    }
    return true; 

}
static bool ray_triangle_intersection(const Ray& ray , const parser::Vec3f & p1 ,const  parser::Vec3f & p2 , const parser::Vec3f &  p3 , const parser::Vec3f  & normal , parser::Vec3f & hitpoint )
{
  // I will do the inefficent method since I am used to it
  
  // 1 - generate plane from 3 points
    bool is_plane_hit = ray_plane_intersection(ray , p1 , normal , hitpoint  );
    if( !is_plane_hit)
    {
        return false; 
    }
    bool is_in_triangle = is_point_in_triangle( p1 ,p2 ,p3 , hitpoint);

    
    return is_in_triangle;
}
static bool ray_sphere_intersection(const Ray& ray ,  parser::Sphere& sphere ,  parser::Vec3f & center , parser::Vec3f &normal ,  parser::Vec3f & hitpoint  )
{
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

    origin_inv =  sphere.transformation_inverse * origin_inv ;
    dir_inv = sphere.transformation_inverse * dir_inv;
    parser::Vec3f origin_inv_3f; 
    parser::Vec3f dir_inv_3f; 

    origin_inv_3f.x = origin_inv.x / origin_inv.w;
    origin_inv_3f.y = origin_inv.y/ origin_inv.w;
    origin_inv_3f.z = origin_inv.z/ origin_inv.w;

    dir_inv_3f.x = dir_inv.x / dir_inv.w;
    dir_inv_3f.y = dir_inv.y/ dir_inv.w;
    dir_inv_3f.z = dir_inv.z/ dir_inv.w;

    Ray inv_ray(origin_inv_3f  , dir_inv_3f  );


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
        return false;
    }
    #endif 
    if( discriminant < 0 )
    {
        return false; 
    }
    
    if( discriminant == 0 )
    {
        float root =  ( -b )  / (2*a);
        if( root < 0 )
        {
            return false; 
        }
        hitpoint = (parser::Vec3f)O +  (parser::Vec3f)k * root;

        //calculate normal 
        normal = parser::normalize( hitpoint - center );


        //multiply with the transform matrix
        parser::Vec4f hitpoint_4f;
        hitpoint_4f.x= hitpoint.x; 
        hitpoint_4f.y= hitpoint.y; 
        hitpoint_4f.z= hitpoint.z; 
        hitpoint_4f.w= 1.0f; 

        parser::Vec4f normal_inv;
        normal_inv.x= normal.x; 
        normal_inv.y= normal.y; 
        normal_inv.z= normal.z; 
        normal_inv.w= 1.0f; 

        hitpoint_4f = sphere.transformation * hitpoint_4f;

        //deneme
        parser::Matrix trans_inverse_transpose = sphere.transformation.inverse();
        trans_inverse_transpose.Transpose();


        normal_inv = trans_inverse_transpose * normal_inv;

        hitpoint.x = hitpoint_4f.x / hitpoint_4f.w;
        hitpoint.y = hitpoint_4f.y / hitpoint_4f.w;
        hitpoint.z = hitpoint_4f.z / hitpoint_4f.w;

        normal.x = normal_inv.x / normal_inv.w;
        normal.y = normal_inv.y / normal_inv.w;
        normal.z = normal_inv.z / normal_inv.w;

        return true;  
    }
    // else there are two roots

    float root_1 =  (-b-std::sqrt(discriminant) )  / (2*a);
    float root_2 = (-b+std::sqrt(discriminant) ) / (2*a);
    if( root_1 < 0 && root_2 < 0  )
    {
        return false; 
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
            hitpoint = h2; 
        }
        else
        {
            hitpoint = h1;
        }
    }
    else
    {
        if( root_1 > 0 )
        {
            hitpoint = h1; 
        }
        else if( root_2 > 0 )
        {
            hitpoint = h2; 
        }
    }
    
    
    //calculate normal
    normal = parser::normalize( hitpoint - center );

    //multiply with the transform matrix
    parser::Vec4f hitpoint_4f;
    hitpoint_4f.x= hitpoint.x; 
    hitpoint_4f.y= hitpoint.y; 
    hitpoint_4f.z= hitpoint.z; 
    hitpoint_4f.w= 1.0f; 

    parser::Vec4f normal_inv;
    normal_inv.x= normal.x; 
    normal_inv.y= normal.y; 
    normal_inv.z= normal.z; 
    normal_inv.w= 1.0f; 

    hitpoint_4f = sphere.transformation * hitpoint_4f;


    //deneme
    parser::Matrix trans_inverse_transpose = sphere.transformation.inverse();
    trans_inverse_transpose.Transpose();
    
    normal_inv = trans_inverse_transpose * normal_inv;

    hitpoint.x = hitpoint_4f.x / hitpoint_4f.w;
    hitpoint.y = hitpoint_4f.y / hitpoint_4f.w;
    hitpoint.z = hitpoint_4f.z / hitpoint_4f.w;

    normal.x = normal_inv.x / normal_inv.w;
    normal.y = normal_inv.y / normal_inv.w;
    normal.z = normal_inv.z / normal_inv.w;
    return true; 

}
//a mesh intersection 
static bool calculate_intersection(parser::Scene& scene ,parser::Mesh& object ,const Ray & ray , parser::Vec3f & intersection_normal ,  parser::Vec3f & intersection_point  )
{
    std::vector<parser::Vec3f> hit_points; // there can be multiple hitpoints for an object 
    std::vector<parser::Vec3f> normals; 
    std::vector<parser::Face> faces; 

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

    origin_inv =  object.transformation_inverse * origin_inv ;
    dir_inv = object.transformation_inverse * dir_inv;
    parser::Vec3f origin_inv_3f; 
    parser::Vec3f dir_inv_3f; 

    origin_inv_3f.x = origin_inv.x / origin_inv.w;
    origin_inv_3f.y = origin_inv.y/ origin_inv.w;
    origin_inv_3f.z = origin_inv.z/ origin_inv.w;

    dir_inv_3f.x = dir_inv.x / dir_inv.w;
    dir_inv_3f.y = dir_inv.y/ dir_inv.w;
    dir_inv_3f.z = dir_inv.z/ dir_inv.w;

    Ray inv_ray(origin_inv_3f  , dir_inv_3f  );
        
    for (size_t i = 0; i < object.faces.size(); i++)
    {
        //get points
        parser::Vec3f p1 = parser::Vec3f( scene.vertex_data[ (object.faces[i].v0_id-1)    ].x , scene.vertex_data[ (object.faces[i].v0_id-1)    ].y , scene.vertex_data[ (object.faces[i].v0_id-1)  ].z); 
        parser::Vec3f p2 = parser::Vec3f( scene.vertex_data[ (object.faces[i].v1_id-1)   ].x , scene.vertex_data[ (object.faces[i].v1_id-1)  ].y , scene.vertex_data[ (object.faces[i].v1_id-1)   ].z);
        parser::Vec3f p3 = parser::Vec3f( scene.vertex_data[ (object.faces[i].v2_id-1)   ].x , scene.vertex_data[ (object.faces[i].v2_id-1)    ].y , scene.vertex_data[ (object.faces[i].v2_id-1)   ].z);

        parser::Vec3f normal = parser::normalize(parser::cross((p2-p1) , (p3-p1)) );
        #ifdef DEBUG 
        if (isnan(normal.x) )
        {
            std::cout << " ray dir " << ray.direction.x << " " << ray.direction.y << " " << ray.direction.z << std::endl; 
             //calculate normal
            std::cout << " new trianglee " << std::endl; 
            std::cout << p1.x << "  " << p1.y << " "  << p1.z << std::endl;
            std::cout << p2.x << "  " << p2.y << " "  << p2.z << std::endl;
            std::cout << p3.x << "  " << p3.y << " "  << p3.z << std::endl;
            std::cout << "indices " << std::endl; 
        }
        #endif

       
       
        
        parser::Vec3f temp_intersection_point(0.0f , 0.0f , 0.f );
        if( ray_triangle_intersection(inv_ray , p1 , p2 , p3 , normal , temp_intersection_point) )
        {
            /*std::cout << " P 1 " << p1.x << " " << p1.y << " " << p1.z << std::endl;
            std::cout << " P 2 " << p2.x << " " << p2.y << " " << p2.z << std::endl;
            std::cout << " P 3 " << p3.x << " " << p3.y << " " << p3.z << std::endl;*/

            //if hit  multiply bacck with transformation 
            parser::Vec4f temp_intersection_real_coord;
            temp_intersection_real_coord.x = temp_intersection_point.x;
            temp_intersection_real_coord.y =  temp_intersection_point.y;
            temp_intersection_real_coord.z =  temp_intersection_point.z;
            temp_intersection_real_coord.w =  1.0f;

            temp_intersection_real_coord = object.transformation * temp_intersection_real_coord;
            temp_intersection_point.x = temp_intersection_real_coord.x / temp_intersection_real_coord.w;
            temp_intersection_point.y = temp_intersection_real_coord.y / temp_intersection_real_coord.w;
            temp_intersection_point.z = temp_intersection_real_coord.z / temp_intersection_real_coord.w;
            // if hit multiply back with transformation 
            parser::Vec4f normal_real_cord;
            normal_real_cord.x = normal.x;
            normal_real_cord.y =  normal.y;
            normal_real_cord.z =  normal.z;
            normal_real_cord.w =  1.0f;

            normal_real_cord =  object.transformation * normal_real_cord ;
            normal.x = normal_real_cord.x / normal_real_cord.w;
            normal.y = normal_real_cord.y / normal_real_cord.w;
            normal.z = normal_real_cord.z / normal_real_cord.w;

            hit_points.push_back(temp_intersection_point);
            normals.push_back(normal);
        }
    }
    if( hit_points.size() == 0 )
    {
        return false; 
    }
    if( hit_points.size() == 1 )
    {
        intersection_point = hit_points[0];
        intersection_normal = normals[0];
        
        return true; 
    }
    std::vector<float> distances;
    for (size_t i = 0; i < hit_points.size(); i++)
    {
        distances.push_back(parser::distance(hit_points[i] , ray.origin) ) ;
    }
    float smallest_length = 99999.0f; 
    int smallest_index = 0; 
    for (size_t i = 0; i < distances.size(); i++)
    {
        if( smallest_length > distances[i])
        {
            smallest_index = i;
            smallest_length = distances[i];
        }
    }
    // fill intersection and normal
    intersection_point = hit_points[smallest_index];
    intersection_normal = normals[smallest_index];

    return true; 
        
    
}
// a sphere intersection 
static bool calculate_intersection(parser::Scene& scene ,parser::Sphere& object , const Ray & ray , parser::Vec3f & intersection_normal ,  parser::Vec3f & intersection_point  )
{
    parser::Vec3f center = parser::Vec3f(scene.vertex_data[object.center_vertex_id-1].x , scene.vertex_data[object.center_vertex_id-1 ].y , scene.vertex_data[object.center_vertex_id-1].z );
    bool is_ray_intersected = ray_sphere_intersection( ray , object , center ,  intersection_normal ,   intersection_point );
    return is_ray_intersected; 
}
// a trianlge intersection 
static bool calculate_intersection(parser::Scene& scene ,parser::Triangle& object ,const Ray & ray , parser::Vec3f & intersection_normal ,  parser::Vec3f & intersection_point  )
{
     //get points
    parser::Vec3f p1 = parser::Vec3f( scene.vertex_data[ (object.indices.v0_id-1)    ].x , scene.vertex_data[ (object.indices.v0_id-1)    ].y , scene.vertex_data[ (object.indices.v0_id-1)  ].z); 
    parser::Vec3f p2 = parser::Vec3f( scene.vertex_data[ (object.indices.v1_id-1)   ].x , scene.vertex_data[ (object.indices.v1_id-1)  ].y , scene.vertex_data[ (object.indices.v1_id-1)   ].z);
    parser::Vec3f p3 = parser::Vec3f( scene.vertex_data[ (object.indices.v2_id-1)   ].x , scene.vertex_data[ (object.indices.v2_id-1)    ].y , scene.vertex_data[ (object.indices.v2_id-1)   ].z);
    
    parser::Vec3f normal = parser::normalize(parser::cross((p2-p1) , (p3-p1)) );

    parser::Vec3f temp_intersection_point(0.0f , 0.0f , 0.f );

    // multiply ray direction and origin with inverse of this
        parser::Vec4f origin_inv; 
        origin_inv.x = ray.origin.x;
        origin_inv.y = ray.origin.y;
        origin_inv.z = ray.origin.z;
        origin_inv.w = 1.0f;
        
        parser::Vec4f dir_inv; 
        dir_inv.x = ray.origin.x + ray.direction.x;
        dir_inv.y = ray.origin.y + ray.direction.y;
        dir_inv.z = ray.origin.z + ray.direction.z;
        dir_inv.w = 1.0f;

        origin_inv = object.transformation_inverse * origin_inv ;
        dir_inv = object.transformation_inverse * dir_inv ;

        parser::Vec3f origin_inv_3f; 
        parser::Vec3f dir_inv_3f; 

        origin_inv_3f.x = origin_inv.x / origin_inv.w;
        origin_inv_3f.y = origin_inv.y/ origin_inv.w;
        origin_inv_3f.z = origin_inv.z/ origin_inv.w;

        dir_inv_3f.x = dir_inv.x / dir_inv.w;
        dir_inv_3f.y = dir_inv.y/ dir_inv.w;
        dir_inv_3f.z = dir_inv.z/ dir_inv.w;

        Ray inv_ray(origin_inv_3f  ,  dir_inv_3f   );
    if( ray_triangle_intersection(inv_ray , p1 , p2 , p3 , normal , temp_intersection_point) )
    {
        //if hit  multiply bacck with transformation 
            parser::Vec4f temp_intersection_real_coord;
            temp_intersection_real_coord.x = temp_intersection_point.x;
            temp_intersection_real_coord.y =  temp_intersection_point.y;
            temp_intersection_real_coord.z =  temp_intersection_point.z;
            temp_intersection_real_coord.w =  1.0f;

            temp_intersection_real_coord = object.transformation * temp_intersection_real_coord ;
            temp_intersection_point.x = temp_intersection_real_coord.x / temp_intersection_real_coord.w;
            temp_intersection_point.y = temp_intersection_real_coord.y / temp_intersection_real_coord.w;
            temp_intersection_point.z = temp_intersection_real_coord.z / temp_intersection_real_coord.w;
            // if hit multiply back with transformation 
            parser::Vec4f normal_real_cord;
            normal_real_cord.x = normal.x;
            normal_real_cord.y =  normal.y;
            normal_real_cord.z =  normal.z;
            normal_real_cord.w =  1.0f;

            normal_real_cord = object.transformation * normal_real_cord ;
            normal.x = normal_real_cord.x / normal_real_cord.w;
            normal.y = normal_real_cord.y / normal_real_cord.w;
            normal.z = normal_real_cord.z / normal_real_cord.w;

            

            intersection_normal = parser::normalize(normal); ;
            intersection_point = temp_intersection_point; 
            return true;
    }
    //false otherwise
    return false;
}
//a mesh instance intersection 
static bool calculate_intersection(parser::Scene& scene ,parser::MeshInstance& object ,const Ray & ray , parser::Vec3f & intersection_normal ,  parser::Vec3f & intersection_point  )
{
    std::vector<parser::Vec3f> hit_points; // there can be multiple hitpoints for an object 
    std::vector<parser::Vec3f> normals; 
    std::vector<parser::Face> faces; 
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

    origin_inv =  object.transformation_inverse * origin_inv ;
    dir_inv = object.transformation_inverse * dir_inv;
    parser::Vec3f origin_inv_3f; 
    parser::Vec3f dir_inv_3f; 

    origin_inv_3f.x = origin_inv.x / origin_inv.w;
    origin_inv_3f.y = origin_inv.y/ origin_inv.w;
    origin_inv_3f.z = origin_inv.z/ origin_inv.w;

    dir_inv_3f.x = dir_inv.x / dir_inv.w;
    dir_inv_3f.y = dir_inv.y/ dir_inv.w;
    dir_inv_3f.z = dir_inv.z/ dir_inv.w;

    Ray inv_ray(origin_inv_3f  , dir_inv_3f  );
    for (size_t i = 0; i < object.mesh_ptr->faces.size(); i++)
    {
        //get points
        parser::Vec3f p1 = parser::Vec3f( scene.vertex_data[ (object.mesh_ptr->faces[i].v0_id-1)    ].x , scene.vertex_data[ (object.mesh_ptr->faces[i].v0_id-1)    ].y , scene.vertex_data[ (object.mesh_ptr->faces[i].v0_id-1)  ].z); 
        parser::Vec3f p2 = parser::Vec3f( scene.vertex_data[ (object.mesh_ptr->faces[i].v1_id-1)   ].x , scene.vertex_data[ (object.mesh_ptr->faces[i].v1_id-1)  ].y , scene.vertex_data[ (object.mesh_ptr->faces[i].v1_id-1)   ].z);
        parser::Vec3f p3 = parser::Vec3f( scene.vertex_data[ (object.mesh_ptr->faces[i].v2_id-1)   ].x , scene.vertex_data[ (object.mesh_ptr->faces[i].v2_id-1)    ].y , scene.vertex_data[ (object.mesh_ptr->faces[i].v2_id-1)   ].z);

        parser::Vec3f normal = parser::normalize(parser::cross((p2-p1) , (p3-p1)) );
        #ifdef DEBUG 
        if (isnan(normal.x) )
        {
            std::cout << " ray dir " << ray.direction.x << " " << ray.direction.y << " " << ray.direction.z << std::endl; 
             //calculate normal
            std::cout << " new trianglee " << std::endl; 
            std::cout << p1.x << "  " << p1.y << " "  << p1.z << std::endl;
            std::cout << p2.x << "  " << p2.y << " "  << p2.z << std::endl;
            std::cout << p3.x << "  " << p3.y << " "  << p3.z << std::endl;
            std::cout << "indices " << std::endl; 
        }
        #endif
        
        
        parser::Vec3f temp_intersection_point(0.0f , 0.0f , 0.f );
        if( ray_triangle_intersection(inv_ray , p1 , p2 , p3 , normal , temp_intersection_point) )
        {
            //if hit  multiply bacck with transformation 
            parser::Vec4f temp_intersection_real_coord;
            temp_intersection_real_coord.x = temp_intersection_point.x;
            temp_intersection_real_coord.y =  temp_intersection_point.y;
            temp_intersection_real_coord.z =  temp_intersection_point.z;
            temp_intersection_real_coord.w =  1.0f;

            temp_intersection_real_coord = object.transformation * temp_intersection_real_coord;
            temp_intersection_point.x = temp_intersection_real_coord.x / temp_intersection_real_coord.w;
            temp_intersection_point.y = temp_intersection_real_coord.y / temp_intersection_real_coord.w;
            temp_intersection_point.z = temp_intersection_real_coord.z / temp_intersection_real_coord.w;
            // if hit multiply back with transformation 
            parser::Vec4f normal_real_cord;
            normal_real_cord.x = normal.x;
            normal_real_cord.y =  normal.y;
            normal_real_cord.z =  normal.z;
            normal_real_cord.w =  1.0f;

            normal_real_cord =  object.transformation * normal_real_cord ;
            normal.x = normal_real_cord.x / normal_real_cord.w;
            normal.y = normal_real_cord.y / normal_real_cord.w;
            normal.z = normal_real_cord.z / normal_real_cord.w;
            hit_points.push_back(temp_intersection_point);
            normals.push_back(normal);
            faces.push_back(object.mesh_ptr->faces[i]);
        }
    }
    if( hit_points.size() == 0 )
    {
        return false; 
    }
    if( hit_points.size() == 1 )
    {
        intersection_point = hit_points[0];
        intersection_normal = normals[0];
        
        return true; 
    }
    std::vector<float> distances;
    for (size_t i = 0; i < hit_points.size(); i++)
    {
        distances.push_back(parser::distance(hit_points[i] , ray.origin) ) ;
    }
    float smallest_length = 99999.0f; 
    int smallest_index = 0; 
    for (size_t i = 0; i < distances.size(); i++)
    {
        if( smallest_length > distances[i])
        {
            smallest_index = i;
            smallest_length = distances[i];
        }
    }
    // fill intersection and normal
    intersection_point = hit_points[smallest_index];
    intersection_normal = normals[smallest_index];

    return true; 
        
    
}
static bool ray_object_intersection( const Ray & ray , parser::Scene & scene , parser::Vec3f &hitpoint , parser::Vec3f & normal , parser::Material &material ,   int &  prev_object_id  , bool is_shadow_rays_active  )
{
    std::vector<parser::Vec3f> hit_points;
    std::vector<parser::Vec3f> normals;
    std::vector<parser::Material> materials; 
    std::vector<float> distances;
    std::vector<int> object_id_list;

    int mesh_count = 0; 
    int triangle_count = 0;  
    int sphere_count = 0;
    int mesh_instances_count = 0;

    int object_id = 0; 

    while( true )
    {
        // if all done s
        
        if( mesh_count == scene.meshes.size() && sphere_count == scene.spheres.size() && triangle_count == scene.triangles.size() && mesh_instances_count == scene.mesh_instances.size() )
        {
            break; 
        }
        parser::Vec3f intersection_point(0.0f , 0.0f ,0.0f);  
        parser::Vec3f intersection_normal(0.0f , 0.0f ,0.0f);  
        bool is_interected_with_this_triangle = false; 
        bool is_intersected_with_this_object = false; 
        bool is_intersected_with_this_sphere = false; 

        while( mesh_count < scene.meshes.size())
        {
            object_id = mesh_count + sphere_count + triangle_count; // give every object unique id 
            is_intersected_with_this_object = calculate_intersection(scene , scene.meshes[mesh_count] , ray , intersection_normal ,   intersection_point );
            if( is_intersected_with_this_object && object_id != prev_object_id ) // the second clause is for shadow rays in order not the intersect itself
            {
                hit_points.push_back(intersection_point ); 
                normals.push_back(intersection_normal );
                materials.push_back(scene.materials[scene.meshes[mesh_count].material_id - 1] );
                object_id_list.push_back(object_id);
            }
            
            mesh_count += 1; 
        }
        while( sphere_count < scene.spheres.size())
        {
            is_intersected_with_this_sphere = calculate_intersection(scene , scene.spheres[sphere_count] ,  ray ,   intersection_normal , intersection_point );
            object_id = mesh_count + sphere_count + triangle_count; // give every object unique id 
            if( is_intersected_with_this_sphere && object_id != prev_object_id  ) // the second clause is for shadow rays in order not the intersect itself
            {
                hit_points.push_back(intersection_point );
                normals.push_back(intersection_normal );
                materials.push_back(scene.materials[scene.spheres[sphere_count].material_id - 1] );
                object_id_list.push_back(object_id);
                
            }
            sphere_count += 1; 

        }

        while( triangle_count < scene.triangles.size() )
        {

            is_interected_with_this_triangle = calculate_intersection(scene , scene.triangles[triangle_count] ,  ray ,   intersection_normal , intersection_point );
            object_id = mesh_count + sphere_count + triangle_count; // give every object unique id 
            if( is_interected_with_this_triangle && object_id != prev_object_id ) // the second clause is for shadow rays in order not the intersect itself
            {
                hit_points.push_back(intersection_point );
                normals.push_back(intersection_normal );
                materials.push_back(scene.materials[scene.triangles[triangle_count].material_id - 1] );
                object_id_list.push_back(object_id);

            }
            triangle_count += 1; 
        }


        while( mesh_instances_count < scene.mesh_instances.size())
        {
            object_id = mesh_count + sphere_count + triangle_count + mesh_instances_count; // give every object unique id 
            is_intersected_with_this_object = calculate_intersection(scene , scene.mesh_instances[mesh_instances_count] , ray , intersection_normal ,   intersection_point );
            if( is_intersected_with_this_object && object_id != prev_object_id ) // the second clause is for shadow rays in order not the intersect itself
            {
                hit_points.push_back(intersection_point ); 
                normals.push_back(intersection_normal );
                materials.push_back(scene.materials[scene.mesh_instances[mesh_instances_count].material_id - 1] );
                object_id_list.push_back(object_id);
            }
            
            mesh_instances_count += 1; 
        }

    }
 
    //now we found the hitpoints get the smallest of them and return black if nothing found
    if( hit_points.size() == 0 )
    {
        return false; 
    }
    for (size_t i = 0; i < hit_points.size(); i++)
    {
        distances.push_back(parser::distance(hit_points[i] , ray.origin));
    }

    float smallest_length = 999999;
    int smallest_index = 0;  
    for (size_t i = 0; i < distances.size(); i++)
    {
        if( smallest_length > distances[i])
        {
            smallest_index = i; 
            smallest_length = distances[i]; 
        }
    }
    hitpoint = hit_points[smallest_index];
    normal = normals[smallest_index];
    material = materials[smallest_index];
    
    if( !is_shadow_rays_active) // if not shadow rays active  set prev_object_no_to object
    {
        prev_object_id = object_id_list[smallest_index];
    }
    return true; 

}

static bool calculate_second_hitpoint_in_same_object( parser::Scene & scene , const Ray & refracted_ray ,  parser::Vec3f & hit_point , parser::Vec3f & normal  , int & object_id , parser::Vec3f & second_hit_point  , parser::Vec3f & second_normal   )
{
    // fetch the object with 
    // 1 - meshes
    // 2 - spheres
    // 3 - triangles
    parser::Mesh mesh;
    parser::Sphere sphere;
    parser::Triangle triangle;
    bool if_mesh = false; 
    bool if_sphere = false; 
    bool if_triangle = false; 

    if( object_id < scene.meshes.size()   )
    {
        if_mesh = true; 
        mesh = scene.meshes[object_id];
    }
    else if( object_id < scene.meshes.size() +  scene.spheres.size() )
    {
        if_sphere = true; 
        sphere = scene.spheres[object_id - scene.meshes.size() ];
    }
    else if( object_id < scene.meshes.size() +  scene.spheres.size() + scene.triangles.size() )
    {
        if_triangle = true; 
        triangle = scene.triangles[object_id - scene.meshes.size()   - scene.triangles.size()  ];
    }
    
    if( if_mesh )
    {

        bool is_intersection = calculate_intersection(scene, mesh , refracted_ray , second_normal , second_hit_point);
        if(!is_intersection)
        {
            //std::cout << " cannot found second hit_point in mesh" << std::endl;
            return false;
        }
        return true;  
    }
    else if( if_sphere )
    {

        bool is_intersection = calculate_intersection(scene, sphere , refracted_ray , second_normal , second_hit_point);
        if(!is_intersection)
        {
            //std::cout << " cannot found second hit_point in sphere" << std::endl;
            return false;

        }
        return true;  
    }
    else if( if_triangle )
    {
        bool is_intersection = calculate_intersection(scene, triangle , refracted_ray , second_normal , second_hit_point);
        if(!is_intersection)
        {
            //std::cout << " cannot found second hit_point in triangle" << std::endl;
            return false;

        }
        return true;  
    }
    else
    {
        std::cout <<  " a problem  has occured " <<  std::endl; 
        return false ;
    }

}


//for bvh
static bool ray_object_intersection( const Ray & ray , parser::Scene & scene , parser::Vec3f &hitpoint , parser::Vec3f & normal , parser::Material &material ,   int &  prev_object_id  , bool is_shadow_rays_active, int & object_id  )
{
    std::vector<parser::Vec3f> hit_points;
    std::vector<parser::Vec3f> normals;
    std::vector<parser::Material> materials; 
    std::vector<float> distances;
    std::vector<int> object_id_list;
    
    bool is_intersected_with_this_object = false; 
    if( scene.meshes.size() > object_id)
    {
        is_intersected_with_this_object = calculate_intersection(scene , scene.meshes[object_id ] , ray , hitpoint ,   normal );
        material = scene.materials[scene.meshes[object_id].material_id - 1 ];
    }
    else if(scene.spheres.size() + scene.meshes.size() > object_id  )
    {
        is_intersected_with_this_object = calculate_intersection(scene , scene.spheres[object_id - scene.meshes.size()] , ray , hitpoint ,   normal );
        material = scene.materials[scene.spheres[object_id - scene.meshes.size()].material_id  - 1 ];

    }
    else if( scene.triangles.size() + scene.spheres.size() + scene.meshes.size() > object_id )
    {
        is_intersected_with_this_object = calculate_intersection(scene , scene.triangles[object_id  - scene.meshes.size() - scene.spheres.size() ] , ray , hitpoint ,   normal );
        material = scene.materials[scene.triangles[object_id  - scene.meshes.size() - scene.spheres.size()].material_id  - 1 ];

    }
    else if( scene.triangles.size() + scene.spheres.size() + scene.meshes.size() + scene.mesh_instances.size() > object_id )
    {
        is_intersected_with_this_object = calculate_intersection(scene , scene.mesh_instances[object_id  - scene.meshes.size()- scene.spheres.size() - scene.mesh_instances.size() ] , ray , hitpoint ,   normal );
        material = scene.materials[scene.mesh_instances[object_id  - scene.meshes.size()- scene.spheres.size() - scene.mesh_instances.size()].material_id  - 1 ];

    }
    if( !is_intersected_with_this_object)
    {
        return false; 
    }
    if( !is_shadow_rays_active) // if not shadow rays active  set prev_object_no_to object
    {
        prev_object_id = object_id;
    }
    
    return true; 
}