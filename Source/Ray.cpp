#pragma once 
#include "../Header/Ray.h"
#include "../System_Files/support_files/parser.h"
#include "../Header/Prototypes.h"
#include <iostream>
Ray::Ray(parser::Vec3f origin , parser::Vec3f point2 )
{
    
    this->origin = origin;
    this->direction = parser::normalize(point2 - origin);

     
}



