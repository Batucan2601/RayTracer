#pragma once
#include "../Dependencies/glm_headers/glm.hpp"
#include "../System_Files/support_files/parser.h"
struct Object
{
    parser::Vec3f color = parser::Vec3f(1.0f , 1.0f , 0.0f ); 
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

    Sphere(const parser::Vec3f & center , const float & radius )
    {
        this->radius = radius;
        this->center = center; 
    }
    parser::Vec3f center; 
    float radius; 
};
class Ray
{
    public:
    Ray(parser::Vec3f origin , parser::Vec3f point2 ); //direction comes out normalized
    //Ray( parser::Vec3f origin , parser::Vec3f  p2  );
    parser::Vec3f origin;
    parser::Vec3f direction;
};

