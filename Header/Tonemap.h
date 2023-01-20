#pragma once
#include "../System_Files/support_files/parser.h"
#include <bits/stdc++.h>
double average_world_luminance(float * tonemap_image , int image_width , int image_height , std::vector<float>& luminances /* 100 luminances */)
{
    double average_world_lum = 0.0f;
    std::vector<float> lum_vec;
    int index = 0;   
    for (size_t y = 0; y < image_height; y++)
    {
        for (size_t x = 0; x < image_width; x++)
        {
            float L_w = tonemap_image[ index++] * 0.2126 + tonemap_image[ index++ ] * 0.7152 + tonemap_image[ index++ ] * 0.0722;  
            
            bool unique_lw = true;
            // for (size_t i = 0; i < lum_vec.size(); i++)
            // {
            //     if( std::abs(lum_vec[i] - L_w )  < 1e-5)
            //     {
            //         unique_lw = false; 
            //         break;
            //     }
            // }
            if( unique_lw)
            {
                lum_vec.push_back(L_w);
            }
            average_world_lum += std::log(L_w + 1e-5); 
        }
    }
    average_world_lum = std::exp( average_world_lum    / (image_width * image_height) ) ;

    std::sort(lum_vec.begin(), lum_vec.end());
 
    
    //fill first 100 luuminances 
    if( lum_vec.size() > 100 )
    {
        for (size_t i = 0; i < lum_vec.size(); i += lum_vec.size()/100 )
        {
            luminances.push_back(lum_vec[i]); 
        }
    }
    
    lum_vec.erase(lum_vec.begin());
    luminances.push_back(lum_vec[lum_vec.size() - 1]);
    return average_world_lum; //equation 1 from reinhard 


}
parser::Vec3f apply_tonemap( parser::Vec3f color ,parser::Camera  & camera  , double average_world_lum  , std::vector<float>& luminances )
{
    parser::ToneMap tonemap = camera.tonemap;
    // 1- compute luminance
    float luminance  = 0.2126 *color.x + 0.7152*color.y + 0.0722*color.z;
    float s = tonemap.saturation;
    // 2 - tonemapping algorithm
    float key = tonemap.tmoOptions[0];
    float burn_percent = tonemap.tmoOptions[1];
    float L_white = 0;
    double scaled_luminance = (key / average_world_lum) * luminance  ; //equation 2 from reinhard 
    
    if( std::abs(luminance) < 1e-6 )
    {
        luminance += 1e-6;
    }

    if( std::abs(burn_percent - 0.0f ) < 1e-5 ) //if burn is 0
    {
        scaled_luminance = scaled_luminance  / (scaled_luminance + 1);
    }
    else
    {
        L_white = luminances[100 - int(burn_percent) - 1];
        scaled_luminance = scaled_luminance * ( 1 + (scaled_luminance/std::pow(L_white, 2)) ) / (1 + scaled_luminance ); 
    }

    parser::Vec3f new_color( scaled_luminance * std::pow( color.x / luminance, s) , scaled_luminance * std::pow(color.y / luminance , s) , scaled_luminance * std::pow( color.z / luminance ,s )); 
    
    //clamp to 0-1 range 
    if( new_color.x > 1.0f )
    {
        new_color.x = 1.0f; 
    }
    else if( new_color.x < 0.0f )
    {
        new_color.x =0.0f;
    }
    if( new_color.y > 1.0f )
    {
        new_color.y = 1.0f; 
    }
    else if( new_color.y < 0.0f )
    {
        new_color.y =0.0f;
    }
    if( new_color.z > 1.0f )
    {
        new_color.z = 1.0f; 
    }
    else if( new_color.z < 0.0f )
    {
        new_color.z =0.0f;
    }
    new_color.x = 255 * std::pow(new_color.x , 1.0f/tonemap.gama);
    new_color.y = 255 * std::pow(new_color.y , 1.0f/tonemap.gama);
    new_color.z = 255 * std::pow(new_color.z , 1.0f/tonemap.gama);



    return new_color;
}