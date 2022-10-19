#include "../Header/Camera.h"
Camera::Camera(glm::vec3 position )
{
    this->position = position; 
}

Camera::~Camera()
{
    delete this->image;
}