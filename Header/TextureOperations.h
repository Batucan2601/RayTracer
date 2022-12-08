#include "Intersections.h"
#include <math.h>       /* isnan, std */
#include "../System_Files/support_files/parser.h"

bool is_texture_present( parser::Scene & scene ,  int object_id , parser::Vec3f & intersection_point  , parser::Vec3f & intersection_normal , parser::Face &intersection_face , parser::Vec3f & texture_color ,  parser::Material &new_material)
{
    if( object_id < scene.meshes.size())
    {
        //object is mesh 
        if( scene.meshes[object_id].textures.size() > 0 )
        {
            get_texture_color_from_mesh( scene ,  object_id ,  intersection_point  ,  intersection_normal , intersection_face ,  texture_color , new_material);
            return true; 
        }
        return false; 
    }
    else if(object_id < scene.meshes.size() + scene.spheres.size() )
    {
        //object is sphere 
        if( scene.spheres[object_id - scene.meshes.size()].textures.size() > 0 )
        {
            get_texture_color_from_sphere(  scene ,   object_id ,  intersection_point  ,  intersection_normal  ,  texture_color , new_material ) ;
            return true; 
        }
        return false; 
    }
    else if( object_id < scene.meshes.size() + scene.spheres.size() + scene.triangles.size() )
    {
        //object is triangle 
    }
    else if(object_id < scene.meshes.size() + scene.spheres.size() + scene.triangles.size() + scene.mesh_instances.size() )
    {
        //object is mesh instance 
    }
}
parser::Vec3f get_texture_color_from_mesh( parser::Scene & scene ,  int object_id , parser::Vec3f & intersection_point  , parser::Vec3f & intersection_normal , parser::Face &intersection_face , parser::Vec3f &  texture_color, parser::Material & new_material  )
{
    // here goes the texture getting stuff
    parser::Mesh *mesh =  &scene.meshes[object_id];
    new_material = scene.materials[mesh->material_id - 1 ]; // material is equal to old material for now 
    std::vector<int> textures = mesh->textures; 
    //intersection_point * mesh_transformation_inverse
    parser::Vec4f intersection_4f;
    intersection_4f.x = intersection_point.x;
    intersection_4f.y = intersection_point.y;
    intersection_4f.z = intersection_point.z;
    intersection_4f.w = 1.0f;
    intersection_4f = mesh->transformation_inverse * intersection_4f;

    parser::Vec3f intersection_object_space;
    intersection_object_space.x = intersection_4f.x / intersection_4f.w;
    intersection_object_space.y = intersection_4f.y / intersection_4f.w;
    intersection_object_space.z = intersection_4f.z / intersection_4f.w;

    //intersection_normal * mesh_transformation_inverse
    parser::Vec4f normal_4f;
    normal_4f.x = intersection_normal.x;
    normal_4f.y = intersection_normal.y;
    normal_4f.z = intersection_normal.z;
    normal_4f.w = 1.0f;
    normal_4f = mesh->transformation_inverse * normal_4f;

    parser::Vec3f normal_object_space;
    normal_object_space.x = normal_4f.x / normal_4f.w;
    normal_object_space.y = normal_4f.y / normal_4f.w;
    normal_object_space.z = normal_4f.z / normal_4f.w;



    for (size_t i = 0; i < textures.size(); i++)
    {
        if( scene.textureMaps[textures[i] - 1 ].is_image ) //image 
        {
            //parse interpolation
            parser::Vec3f v1 = scene.vertex_data[intersection_face.v0_id-1]; //A 
            parser::Vec3f v2 = scene.vertex_data[intersection_face.v1_id-1]; // B 
            parser::Vec3f v3 = scene.vertex_data[intersection_face.v2_id-1];// C
            
            int v1_texcoord[2]; // a
            v1_texcoord[0] = scene.texcoords[intersection_face.v0_id-1][0]; 
            v1_texcoord[1] = scene.texcoords[intersection_face.v0_id-1][1];  

            int v2_texcoord[2]; // b 
            v2_texcoord[0] = scene.texcoords[intersection_face.v1_id-1][0];
            v2_texcoord[1] = scene.texcoords[intersection_face.v1_id-1][1];

            int v3_texcoord[2]; // c 
            v3_texcoord[0] = scene.texcoords[intersection_face.v2_id-1][0];
            v3_texcoord[1] = scene.texcoords[intersection_face.v2_id-1][1];

            //barycentric calculation 
            // float den = 1 /( (v2.y - v3.y) * (v1.x - v3.x) + (v3.x - v2.x) * (v1.y - v3.y));
            // float alpha = ( (v2.y - v3.y) * (intersection_point.x - v3.x) + (v3.x - v2.x) * (intersection_point.y - v3.y) ) * den;
            // float beta = ( (v3.y - v1.y) * (intersection_point.x - v3.x) + (v1.x - v3.x) * (intersection_point.y - v3.y) ) * den;
            // float gama = 1 - alpha - beta; 


            float alpha =  parser::length(parser::cross_non_normalize( v3 - v2 ,intersection_object_space - v2  )) / parser::length(parser::cross_non_normalize( v2 - v1  , v3 - v1  ) )  ;
            
            float beta = parser::length(parser::cross_non_normalize( v3 - v1 ,intersection_object_space - v1  ) ) / parser::length(parser::cross_non_normalize( v2 - v1  , v3 - v1  ) )  ;
            
            float gama = parser::length(parser::cross_non_normalize( v2 - v1 ,intersection_object_space - v1  ) ) / parser::length(parser::cross_non_normalize( v2 - v1  , v3 - v1  ) )  ; 
            // if( std::abs(gama - 1) < 1e-5   )
            // {
            //     gama = 1.0f; 
            // }
            // else if(std::abs(gama) < 1e-5 )
            // {
            //     gama = 0.0f; 
            // }

            //now we have found the barycentirc coordinates of hitpoint

            //now find u v coord  
            float hit_texcoord[2];
            hit_texcoord[0] = alpha * v1_texcoord[0] + beta * v2_texcoord[0] + gama * v3_texcoord[0];
            hit_texcoord[1] = alpha * v1_texcoord[1] + beta * v2_texcoord[1] + gama * v3_texcoord[1]; 
            if(scene.textureMaps[textures[i] - 1 ].interpolation == "bilinear" )
            {
                texture_color = get_bilinear_coord_color(scene.textureMaps[textures[i] - 1 ] , hit_texcoord[0] ,hit_texcoord[1]);
                
            }
            else if(scene.textureMaps[textures[i] - 1 ].interpolation == "nearest")
            {
                texture_color = get_nearest_coord_color(scene.textureMaps[textures[i] - 1 ] , hit_texcoord[0] ,hit_texcoord[1]);
            }

            //parse Decal
            if(scene.textureMaps[textures[i] - 1 ].decalMode == "replace_kd" )
            {
                new_material.diffuse = texture_color /255; 
            }
            else if(scene.textureMaps[textures[i] - 1 ].decalMode == "blend_kd")
            {
                //new_material.diffuse =  ( texture_color/255 + new_material.diffuse  )  /2 ; 
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "replace_ks" )
            {
               // new_material.specular =  texture_color/255  ; 
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "replace_normal" )
            {
                
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "bump_normal" )
            {
                
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "replace_all" )
            {
                
            }
            

            

        }
        else //perlin
        {

        }

    }
    

}
parser::Vec3f get_texture_color_from_sphere( parser::Scene & scene ,  int object_id , parser::Vec3f & intersection_point  , parser::Vec3f & intersection_normal  , parser::Vec3f &  texture_color ,  parser::Material & new_material )
{
    // here goes the texture getting stuff
    parser::Sphere *sphere =  &scene.spheres[object_id];
    new_material = scene.materials[sphere->material_id - 1 ]; // material is equal to old material for now 
    std::vector<int> textures = sphere->textures; 
    //intersection_point * mesh_transformation_inverse
    parser::Vec4f intersection_4f;
    intersection_4f.x = intersection_point.x;
    intersection_4f.y = intersection_point.y;
    intersection_4f.z = intersection_point.z;
    intersection_4f.w = 1.0f;
    intersection_4f = sphere->transformation_inverse * intersection_4f;

    parser::Vec3f intersection_object_space;
    intersection_object_space.x = intersection_4f.x / intersection_4f.w;
    intersection_object_space.y = intersection_4f.y / intersection_4f.w;
    intersection_object_space.z = intersection_4f.z / intersection_4f.w;

    //intersection_normal * mesh_transformation_inverse
    parser::Vec4f normal_4f;
    normal_4f.x = intersection_normal.x;
    normal_4f.y = intersection_normal.y;
    normal_4f.z = intersection_normal.z;
    normal_4f.w = 1.0f;
    normal_4f = sphere->transformation_inverse * normal_4f;

    parser::Vec3f normal_object_space;
    normal_object_space.x = normal_4f.x / normal_4f.w;
    normal_object_space.y = normal_4f.y / normal_4f.w;
    normal_object_space.z = normal_4f.z / normal_4f.w;



    for (size_t i = 0; i < textures.size(); i++)
    {
        if( scene.textureMaps[textures[i] - 1 ].is_image ) //image 
        {
            float omega = std::acos(intersection_object_space.y /  sphere->radius); 
            float phi = std::atan2(intersection_object_space.z ,  intersection_object_space.x); 

            float hit_texcoord[2];
            hit_texcoord[0] = (-phi + M_PI) / (2 * M_PI);
            hit_texcoord[1] = (-phi + M_PI) / (2 * M_PI);

            if(scene.textureMaps[textures[i] - 1 ].interpolation == "bilinear" )
            {
                texture_color = get_bilinear_coord_color(scene.textureMaps[textures[i] - 1 ] , hit_texcoord[0] ,hit_texcoord[1]);
                
            }
            else if(scene.textureMaps[textures[i] - 1 ].interpolation == "nearest")
            {
                texture_color = get_nearest_coord_color(scene.textureMaps[textures[i] - 1 ] , hit_texcoord[0] ,hit_texcoord[1]);
            }

            //parse Decal
            if(scene.textureMaps[textures[i] - 1 ].decalMode == "replace_kd" )
            {
                new_material.diffuse = texture_color /255; 
            }
            else if(scene.textureMaps[textures[i] - 1 ].decalMode == "blend_kd")
            {
                //new_material.diffuse =  ( texture_color/255 + new_material.diffuse  )  /2 ; 
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "replace_ks" )
            {
               // new_material.specular =  texture_color/255  ; 
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "replace_normal" )
            {
                
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "bump_normal" )
            {
                
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "replace_all" )
            {
                
            }
            

            

        }
        else //perlin
        {

        }

    }
}

parser::Vec3f get_bilinear_coord_color(  parser::TextureMap & texturemap , float u , float v )
{
    parser::Image * img =  texturemap.image; 
    float i = u * img->w;
    float j = v * img->h;
    int p = std::floor(i);
    int q = std::floor(j);
    float dx = i - p;
    float dy = j - q; 

    parser::Vec3f color_left_bot;
    color_left_bot.x = img->image[ img->comp * ( q * img->w + p)  ];
    color_left_bot.y = img->image[ img->comp * ( q * img->w + p) + 1  ];
    color_left_bot.z = img->image[ img->comp * ( q * img->w + p) + 2  ];

    parser::Vec3f color_left_top;
    color_left_bot.x = img->image[ img->comp * ( (q+1) * img->w + p)  ];
    color_left_bot.y = img->image[ img->comp * ( (q+1) * img->w + p) + 1  ];
    color_left_bot.z = img->image[ img->comp * ( (q+1) * img->w + p) + 2  ];

    parser::Vec3f color_right_top;
    color_right_top.x = img->image[ img->comp * ( (q+1) * img->w + p+1)  ];
    color_right_top.y = img->image[ img->comp * ( (q+1) * img->w + p+1) + 1  ];
    color_right_top.z = img->image[ img->comp * ( (q+1) * img->w + p+1) + 2  ];

    parser::Vec3f color_right_bot;
    color_right_bot.x = img->image[ img->comp * ( (q) * img->w + p+1)  ];
    color_right_bot.y = img->image[ img->comp * ( (q) * img->w + p+1) + 1  ];
    color_right_bot.z = img->image[ img->comp * ( (q) * img->w + p+1) + 2  ];

    return  color_left_bot * (1 - dx )* (1 - dy ) + color_right_bot * (dx) * (1-dy) + color_left_top *(1 - dx) * (dy) + color_right_top * dx * dy; 
}
parser::Vec3f get_nearest_coord_color(  parser::TextureMap & texturemap , float u , float v )
{
    // 1 - get image
    parser::Image * img =  texturemap.image; 
    int pix_x = std::floor(u * img->w);
    int pix_y = std::floor(v * img->h);

    parser::Vec3f color;
    color.x = img->image[ img->comp * ( pix_y * img->w + pix_x)  ];
    color.y = img->image[ img->comp * ( pix_y * img->w + pix_x) + 1  ];
    color.z = img->image[ img->comp * ( pix_y * img->w + pix_x) + 2  ];

    return color;


}