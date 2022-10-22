#include "Renderer.h"
// intersections and stuff 
static Plane generate_plane(const glm::vec3 & p1 , const glm::vec3 & normal )
{
    //calculate normal
    
}
static bool ray_triangle_intersection( const Ray &ray , const parser::Scene& scene ,  glm::vec3 intersection_point  )
{
     
}
static bool ray_plane_intersection(const Ray& ray , const Plane& plane , glm::vec3 & hitpoint  )
{

}
static bool ray_triangle_intersection(const Ray& ray , const glm::vec3 & p1 ,const  glm::vec3 & p2 , const glm::vec3 &  p3 , const glm::vec3  & normal , glm::vec3 & hitpoint )
{
  // I will do the inefficent method since I am used to it
  
  // 1 - generate plane from 3 points
  Plane p = generate_plane(p1 , normal );
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


//mesh intersection 
static bool calculate_intersection(parser::Scene& scene ,parser::Mesh& object ,  glm::vec3 & intersection_point )
{
    
    for (size_t i = 0; i < object.faces.size(); i++)
    {
        //get points
        glm::vec3 p1 = glm::vec3( scene.vertex_data[ (object.faces[i].v0_id-1) / 3  ].x , scene.vertex_data[ (object.faces[i].v0_id-1) / 3  ].y , scene.vertex_data[ (object.faces[i].v0_id-1) / 3  ].z); 
        glm::vec3 p2 = glm::vec3( scene.vertex_data[ (object.faces[i].v1_id-1) / 3  ].x , scene.vertex_data[ (object.faces[i].v1_id-1) / 3  ].y , scene.vertex_data[ (object.faces[i].v1_id-1) / 3  ].z);
        glm::vec3 p3 = glm::vec3( scene.vertex_data[ (object.faces[i].v2_id-1) / 3  ].x , scene.vertex_data[ (object.faces[i].v2_id-1) / 3  ].y , scene.vertex_data[ (object.faces[i].v2_id-1) / 3  ].z);
        //calculate normal
        glm::vec3 normal = glm::normalize(glm::cross((p2-p1) , (p3-p1)) );
        ray_triangle_intersection();

    }
    
    return true; 
}
static bool calculate_intersection(parser::Scene& scene ,parser::Sphere& object ,  glm::vec3 & intersection_point )
{
    return true; 
}