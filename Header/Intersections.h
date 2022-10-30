#include "Renderer.h"
#include <math.h>       /* isnan, std */

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

    /*parser::Vec3f AB  = p2 - p1; 
    parser::Vec3f BC  = p3 - p2;
    parser::Vec3f CA  = p1 - p3;
    parser::Vec3f AP  = hitpoint - p1;
    parser::Vec3f BP  = hitpoint - p2;
    parser::Vec3f CP  = hitpoint - p3;


    parser::Vec3f ABPA = parser::cross(AB , AP);
    parser::Vec3f BCPB = parser::cross(BC , BP);
    parser::Vec3f CAPC = parser::cross(CA , CP );

    if( ! ( (ABPA.x >= 0.0f && BCPB.x >= 0.0calculate_intersection
       // all not same sign means it is 
       return false; 
    }
    if( ! ( (ABPA.y >= 0.0f && BCPB.y >= 0.0f && CAPC.y >= 0.0f ) || (ABPA.y <= 0.0f && BCPB.y <= 0.0f && CAPC.y <= 0.0f )  ) )
    {
       // all not same sign means it is 
       return false; 
    }
    if( ! ( (ABPA.z >= 0.0f && BCPB.z >= 0.0f && CAPC.z >= 0.0f ) || (ABPA.z <= 0.0f && BCPB.z <= 0.0f && CAPC.z <= 0.0f )  ) )
    {
       // all not same sign means it is 
       return false; 
    }
    return true;*/
    
    parser::Vec3f normal = parser::cross((parser::Vec3f)p2 - (parser::Vec3f)p1 , (parser::Vec3f)p3 - (parser::Vec3f)p1 );
    
    parser::Vec3f C1 = parser::cross( (parser::Vec3f)p2 - (parser::Vec3f)p1 , (parser::Vec3f)hitpoint - (parser::Vec3f)p1  );
    
    if( parser::dot( C1 , normal) < -1e-3 )
    {
        return false;
    }
    parser::Vec3f C2 = parser::cross( (parser::Vec3f)p3 - (parser::Vec3f)p2 , (parser::Vec3f)hitpoint - (parser::Vec3f)p2  );
    if( parser::dot( C2 , normal) < -1e-3 )
    {

        return false;
    }
    parser::Vec3f C3 = parser::cross( (parser::Vec3f)p1 - (parser::Vec3f)p3 , (parser::Vec3f)hitpoint - (parser::Vec3f)p3  );
    
    if( parser::dot( C3 , normal) < -1e-3 )
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
static bool ray_sphere_intersection(const Ray& ray , const parser::Sphere& sphere , const parser::Vec3f & center , parser::Vec3f &normal ,  parser::Vec3f & hitpoint  )
{
    // at^2 + bT + c = 0
    parser::Vec3f O = ray.origin; 
    parser::Vec3f k = ray.direction; 

    /*float a  = 3.0f;
    float b = 2 * ( O.x * k.x + O.y * k.y + O.z * k.z ) -2*(k.x * center.x + k.y * center.y + k.z * center.z ); 
    float c = parser::dot(O , O ) + parser::dot(center,center) + -2*(O.x*center.x + O.y*center.y + O.z*center.z) -3 * sphere.radius * sphere.radius;
    */

    float a = parser::dot(k , k );
    float b = 2 * parser::dot(k , O-center );
    float c = parser::dot(O-center , O-center) - (sphere.radius * sphere.radius );
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
    return true; 

}


//a mesh intersection 
static bool calculate_intersection(parser::Scene& scene ,parser::Mesh& object ,const Ray & ray , parser::Vec3f & intersection_normal ,  parser::Vec3f & intersection_point  )
{
    std::vector<parser::Vec3f> hit_points; // there can be multiple hitpoints for an object 
    std::vector<parser::Vec3f> normals; 
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
        if( ray_triangle_intersection(ray , p1 , p2 , p3 , normal , temp_intersection_point) )
        {
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
    if( ray_triangle_intersection(ray , p1 , p2 , p3 , normal , temp_intersection_point) )
    {
        intersection_normal = normal;
        intersection_point = temp_intersection_point; 
        return true;
    }
    //false otherwise
    return false; 
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
    int object_id = 0; 

    while( true )
    {
        // if all done s
        
        if( mesh_count == scene.meshes.size() && sphere_count == scene.spheres.size() && triangle_count == scene.triangles.size() )
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
            std::cout << " cannot found second hit_point in mesh" << std::endl;
            return false;
        }
        return true;  
    }
    else if( if_sphere )
    {

        bool is_intersection = calculate_intersection(scene, sphere , refracted_ray , second_normal , second_hit_point);
        if(!is_intersection)
        {
            std::cout << " cannot found second hit_point in sphere" << std::endl;
            return false;

        }
        return true;  
    }
    else if( if_triangle )
    {
        bool is_intersection = calculate_intersection(scene, triangle , refracted_ray , second_normal , second_hit_point);
        if(!is_intersection)
        {
            std::cout << " cannot found second hit_point in triangle" << std::endl;
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
