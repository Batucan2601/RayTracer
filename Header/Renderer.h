#pragma once
#include "../Dependencies/glm_headers/glm.hpp"
#include "Camera.h" 
#include "Ray.h"
#include "Scene.h"

//assume the direction is always +y 
class Renderer
{
    public:
    Renderer();
    void set_h_angle(const float&  angle);
    void set_w_angle(const float&  angle);
    void generate_camera(const glm::vec3& origin);
    void generate_image(const glm::vec3 & mid_point , const  int &  width ,const int&   height, const int &  channel_no ,const char* filename  );
    void ray_tracer();

    Camera* camera;
    Image* image; 


};
