#pragma once 
#include "../Header/Ray.h"
#include "../System_Files/support_files/parser.h"
#include "../Header/Prototypes.h"
#include <math.h>       /* isnan, std */
#include <random>
#include <iostream>

Ray::Ray(parser::Vec3f origin , parser::Vec3f point2 )
{
    
    this->origin = origin;
    this->direction = parser::normalize(point2 - origin);
    
     
}

void Ray::generate_time()
{
    // uniform dist for motion blur
    std::random_device rd_mot_blur;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen_mot_blur(rd_mot_blur()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis_mot_blur(0 , 1);
    this->time_parameter = dis_mot_blur(gen_mot_blur);
}


