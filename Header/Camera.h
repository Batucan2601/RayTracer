#pragma once
#include "../Dependencies/glm_headers/glm.hpp"
#include "Image.h" 
class Camera
{
    public:
    Camera(glm::vec3 position);
    ~Camera();
    
    //variables
    glm::vec3 position; 
    float vertical_angle = 45.0f;
    float horizontal_angle = 45.0f;

    glm::vec3 up = glm::vec3( 0.0f , 1.0f , 0.0f );
    glm::vec3 front;
    Image* image; 
};