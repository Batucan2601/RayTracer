#pragma once
#include "../System_Files/support_files/parser.h"
#include <random>


#define TINYEXR_IMPLEMENTATION
#include "../System_Files/support_files/tinyexr.h"

std::random_device rd_env;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen_rd(rd_env()); // Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis_env(0 ,1);


float* hdr_out; // width * height * RGBA
int hdr_width;
int hdr_height;
const char* hdr_err = NULL; // or nullptr in C++11  

void sphere_env_open_hdr(std::string hdr_img_name)
{
    
    int ret = LoadEXR(&hdr_out, &hdr_width, &hdr_height, hdr_img_name.c_str(), &hdr_err);
    if (ret != TINYEXR_SUCCESS) {
    if (hdr_err) {
       fprintf(stderr, "ERR : %s\n", hdr_err);
       FreeEXRErrorMessage(hdr_err); // release memory of error message.
    }
  }
}
void release_sphere_env()
{
    free(hdr_out);
}
parser::Vec3f get_random_vector_rejection_sampling(parser::Vec3f &normal )
{
    while( true)
    {
        // 1- generate random vector
        float x = dis_env(gen_rd);
        float y = dis_env(gen_rd);
        float z = dis_env(gen_rd);
        // put them to -1 0 
        x = 2* x - 1; 
        y = 2* y - 1; 
        z = 2* z - 1; 
        
        parser::Vec3f rand_vec( x ,y ,z );
        if(  parser::length(rand_vec) <= 1 && parser::dot( rand_vec , normal) > 0 )
        {
            return parser::normalize( rand_vec );
        }
    }   

}
parser::Vec3f get_hdr_image_color(parser::Vec3f &l   )
{
   float u = (  (1 + std::atan2(l.x , -1 * l.z)/M_PI ) /  2);
   float v = (  std::acos(l.y) / M_PI);

    float i = u * hdr_width;
    float j = v * hdr_height;
    int p = std::floor(i);
    int q = std::floor(j);
    float dx = i - p;
    float dy = j - q; 

    parser::Vec3f color_left_bot;
    color_left_bot.x =  hdr_out[ 4 * ( q * hdr_width + p)  ];
    color_left_bot.y =  hdr_out[ 4 * ( q * hdr_width + p) + 1  ];
    color_left_bot.z =  hdr_out[ 4 * ( q * hdr_width + p) + 2  ];

    color_left_bot = color_left_bot * (float)( 2.0f * M_PI ) ; 

    return  color_left_bot ;
}

parser::Vec3f spherical_background(parser::Vec3f & ray_direction)
{
    parser::Vec3f ray_dir_normalized = parser::normalize(ray_direction);
    return get_hdr_image_color(ray_dir_normalized);
}