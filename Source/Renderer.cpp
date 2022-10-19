#include "../Header/Renderer.h"
#include "../Header/Ray.h"
Renderer::Renderer()
{

}
void Renderer::generate_camera(const glm::vec3& origin)
{
    this->camera = new Camera(origin);
}

void Renderer::generate_image(const glm::vec3 & mid_point , const  int &  width ,const int&   height, const int &  channel_no ,const char* filename  )
{
    this->image = new Image( mid_point ,  width ,height,  channel_no , filename );
}
void Renderer::set_h_angle(const float&  angle)
{
    this->camera->horizontal_angle = angle; 
}
void Renderer::set_w_angle(const float&  angle)
{
    this->camera->vertical_angle = angle; 
}
void Renderer::ray_tracer()
{
    //get 0 point
    glm::vec3 image_mid_point = this->camera->position + this->image->image_center;  // world oordinate of the mid point
    //now go to the 0'th point (0 , 0 )
    float left = glm::tan(glm::radians(this->camera->horizontal_angle/2.0f)) *  glm::length(this->image->image_center); // tan(alpha) * a 
    float up = glm::tan(glm::radians(this->camera->vertical_angle/2.0f)) *  glm::length(this->image->image_center); // tan(alpha) * a 
    glm::vec3 up_left_point = glm::vec3(image_mid_point.x - left , image_mid_point.y , image_mid_point.z +  up ); // found 0 , 0 
    for (size_t x = 0; x < this->image->width; x++)
    {
        for (size_t y = 0; y < this->image->height; y++)
        {
            // calculate the ray
            //first calculate the pixel world position
           
            
        }
        
    }

}