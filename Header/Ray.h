#pragma once
#include "../Dependencies/glm_headers/glm.hpp"
#include "../System_Files/support_files/parser.h"
struct Object
{
    glm::vec3 color = glm::vec3(1.0f , 1.0f , 0.0f ); 
};
struct Plane : Object
{
    Plane(const float & A , const float & B , const float & C , const float & D )
    {
        this->A = A;
        this->B = B;
        this->C = C;
        this->D = D;
    }
    float A;
    float B;
    float C;
    float D;
    // Ax + By + Cz + D = 0
};

struct Sphere : Object
{

    Sphere(const glm::vec3 & center , const float & radius )
    {
        this->radius = radius;
        this->center = center; 
    }
    glm::vec3 center; 
    float radius; 
};
class Ray
{
    public:
    Ray(glm::vec3 origin , glm::vec3 point2 ); //direction comes out normalized
    //Ray( glm::vec3 origin , glm::vec3  p2  );
    glm::vec3 origin;
    glm::vec3 direction;
};

