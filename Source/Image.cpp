#include "../Header/Image.h"

Image::Image( const glm::vec3 & image_center ,   const int & height , const int &  width ,const int & channel_no , const char* filename   )
{
    this->height = height;
    this->width = width;
    this->channel_no = channel_no;
    this->filename = filename; 
    this->img = new unsigned char[height * width * channel_no ];
    this->image_center = image_center;
}
Image::~Image()
{
    delete[] this->img; 
}