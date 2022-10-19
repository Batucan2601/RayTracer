#pragma once
#include "../Dependencies/glm_headers/glm.hpp"
class Image
{
    public:
    Image( const glm::vec3 & image_center ,   const int & height , const int &  width ,const int & channel_no , const char* filename   );
    ~Image();

    //variable 
    int height; 
    int width;
    int channel_no;
    const char * filename; 
    unsigned char* img;  // height * width * channel 3 ,
    glm::vec3 image_center;  // image center is relative to camera's space 
};