#include "Renderer.h"
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