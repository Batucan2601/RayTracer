#include "Intersections.h"
#include <math.h>       /* isnan, std */
#include "../System_Files/support_files/parser.h"
#include <algorithm>
static int permutation_table[256]; // ap ermutation table for perlin noise 

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
            else 
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
                 new_material.diffuse =  ( texture_color/255 + new_material.diffuse  )  /2 ; 
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "replace_ks" )
            {
                new_material.specular =  texture_color/255  ; 
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "replace_normal" )
            {
                intersection_normal = replace_normal_mesh(v1_texcoord ,v2_texcoord ,v3_texcoord , v1 ,v2,v3, texture_color );
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "bump_normal" )
            {
                intersection_normal = replace_normal_mesh( v1_texcoord , v2_texcoord , v3_texcoord , v1 ,v2,v3, texture_color  , scene.textureMaps[textures[i] - 1 ] ,  hit_texcoord[0], hit_texcoord[1]);
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
    texture_color.x = 0; 
    texture_color.y = 0; 
    texture_color.z = 0; 

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

    parser::Vec3f center = scene.vertex_data[sphere->center_vertex_id-1];
    intersection_object_space = intersection_object_space - center; 
    for (size_t i = 0; i < textures.size(); i++)
    {
        if( scene.textureMaps[textures[i] - 1 ].is_image ) //image 
        {
            float omega = std::acos(intersection_object_space.y /  sphere->radius); 
            float phi = std::atan2(intersection_object_space.z ,  intersection_object_space.x); 

            float hit_texcoord[2];
            hit_texcoord[0] = (-phi + M_PI) / (2 * M_PI);
            hit_texcoord[1] = (omega) / ( M_PI);

            if(scene.textureMaps[textures[i] - 1 ].interpolation == "bilinear" )
            {
                texture_color = texture_color + get_bilinear_coord_color(scene.textureMaps[textures[i] - 1 ] , hit_texcoord[0] ,hit_texcoord[1]);
                
            }
            else if(scene.textureMaps[textures[i] - 1 ].interpolation == "nearest")
            {
                texture_color = texture_color + get_nearest_coord_color(scene.textureMaps[textures[i] - 1 ] , hit_texcoord[0] ,hit_texcoord[1]);
            }
            else 
            {
                texture_color = texture_color + get_nearest_coord_color(scene.textureMaps[textures[i] - 1 ] , hit_texcoord[0] ,hit_texcoord[1]);
            }
            //parse Decal
            if(scene.textureMaps[textures[i] - 1 ].decalMode == "replace_kd" )
            {
                new_material.diffuse = texture_color /255; 
            }
            else if(scene.textureMaps[textures[i] - 1 ].decalMode == "blend_kd")
            {
                new_material.diffuse =  ( texture_color/255 + new_material.diffuse  )  /2 ; 
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "replace_ks" )
            {
                new_material.specular =  texture_color/255  ; 
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "replace_normal" )
            {
                // intersction point saibeli !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               intersection_normal = replace_normal_sphere( sphere , intersection_point ,   texture_color  ,  omega ,  phi );
               int a = 1;
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "bump_normal" )
            {
                intersection_normal = replace_normal_sphere_bump(sphere , intersection_point ,   texture_color  ,  omega ,  phi ,scene.textureMaps[textures[i] - 1 ] , hit_texcoord[0]  , hit_texcoord[1] );
            }
            else if( scene.textureMaps[textures[i] - 1 ].decalMode == "replace_all" )
            {
                
            }
            

            

        }
        else //perlin
        {
            //get the box
            parser::Vec3f center = scene.vertex_data[sphere->center_vertex_id-1];

            float left_x = center.x - sphere->radius;
            float right_x = center.x + sphere->radius;
            float bot_y = center.y - sphere->radius;
            float up_y = center.y + sphere->radius;
            float near_z = center.z - sphere->radius;
            float far_z = center.z + sphere->radius;



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
    color_left_bot.x = img->image[ 3 * ( q * img->w + p)  ];
    color_left_bot.y = img->image[ 3 * ( q * img->w + p) + 1  ];
    color_left_bot.z = img->image[ 3 * ( q * img->w + p) + 2  ];

    parser::Vec3f color_left_top;
    color_left_top.x = img->image[ 3 * ( (q-1) * img->w + p)  ];
    color_left_top.y = img->image[ 3 * ( (q-1) * img->w + p) + 1  ];
    color_left_top.z = img->image[ 3 * ( (q-1) * img->w + p) + 2  ];

    parser::Vec3f color_right_top;
    color_right_top.x = img->image[ 3 * ( (q-1) * img->w + p+1)  ];
    color_right_top.y = img->image[ 3 * ( (q-1) * img->w + p+1) + 1  ];
    color_right_top.z = img->image[ 3 * ( (q-1) * img->w + p+1) + 2  ];

    parser::Vec3f color_right_bot;
    color_right_bot.x = img->image[ 3 * ( (q) * img->w + p+1)  ];
    color_right_bot.y = img->image[ 3 * ( (q) * img->w + p+1) + 1  ];
    color_right_bot.z = img->image[ 3 * ( (q) * img->w + p+1) + 2  ];

    return  color_left_bot * (1 - dx )* (1 - dy ) + color_right_bot * (dx) * (1-dy) + color_left_top *(1 - dx) * (dy) + color_right_top * dx * dy; 
}
parser::Vec3f get_nearest_coord_color(  parser::TextureMap & texturemap , float u , float v )
{
    // 1 - get image
    parser::Image * img =  texturemap.image; 
    int pix_x = std::floor(u * img->w);
    int pix_y = std::floor(v * img->h);

    parser::Vec3f color;
    color.x = img->image[ 3 * ( pix_y * img->w + pix_x)  ];
    color.y = img->image[ 3 * ( pix_y * img->w + pix_x) + 1  ];
    color.z = img->image[ 3 * ( pix_y * img->w + pix_x) + 2  ];

    return color;


}


parser::Vec3f replace_normal_sphere(parser::Sphere * sphere , parser::Vec3f & hitpoint ,  parser::Vec3f & textureColor  , float omega , float phi )
{
    //get x y z again
    parser::Vec3f new_normal;
    
    parser::Vec3f T;
    T.x = 2 * M_PI * hitpoint.z;
    T.y = 0;
    T.z = -2 * M_PI * hitpoint.x;
    T = parser::normalize(T);

    parser::Vec3f B;
    B.x = M_PI * std::cos(  phi * M_PI/180) * hitpoint.y;
    B.y = -sphere->radius * M_PI * std::sin(omega * M_PI / 180);
    B.z = M_PI * hitpoint.y * std::sin(phi * M_PI/180);
    B = parser::normalize(B);

    parser::Vec3f texture_normal_from_color;
    parser::Vec3f ones;
    ones.x = 1;
    ones.y = 1; 
    ones.z = 1; 

    texture_normal_from_color = parser::normalize(textureColor / 127.5f - ones ); 
    parser::Vec3f N = parser::cross(B ,T);
    N = parser::normalize(N );
    new_normal.x = T.x * texture_normal_from_color.x + B.x * texture_normal_from_color.y + N.x * texture_normal_from_color.z;
    new_normal.y = T.y * texture_normal_from_color.x + B.y * texture_normal_from_color.y + N.y * texture_normal_from_color.z;
    new_normal.z = T.z * texture_normal_from_color.x + B.z * texture_normal_from_color.y + N.z * texture_normal_from_color.z;

    return parser::normalize(new_normal);
}

parser::Vec3f replace_normal_mesh(int* v1_texcoord ,int* v2_texcoord ,int* v3_texcoord , parser::Vec3f &v1 ,parser::Vec3f &v2,parser::Vec3f &v3, parser::Vec3f & textureColor )
{
    float delta1[2];
    float delta2[2];

    //texture space 
    delta1[0] = v2_texcoord[0] - v1_texcoord[0];
    delta1[1] = v2_texcoord[1] - v1_texcoord[1];

    delta2[1] = v3_texcoord[1] - v2_texcoord[1]; 
    delta2[0] = v3_texcoord[0] - v2_texcoord[0];

    //world space 
    parser::Vec3f E1 = v2 - v1; 
    parser::Vec3f E2 = v3 - v2; 

    float inv_delta[4]; // row major a =0 b  = 1 c = 2 d = 3 
    float den = 1 / (delta1[0] * delta2[1] - delta1[1] * delta2[0]);
    inv_delta[0] = delta2[1] * den; 
    inv_delta[1] = - delta1[1] * den;
    inv_delta[2] = - delta2[0] * den;
    inv_delta[3] =  delta1[0] * den;


    // now matrix mul 
    parser::Vec3f T;
    T.x = inv_delta[0] * E1.x + inv_delta[1] * E2.x;
    T.y = inv_delta[0] * E1.y + inv_delta[1] * E2.y;
    T.z = inv_delta[0] * E1.z + inv_delta[1] * E2.z;

    parser::Vec3f B;
    B.x = inv_delta[2] * E1.x + inv_delta[3] * E2.x;
    B.y = inv_delta[2] * E1.y + inv_delta[3] * E2.y;
    B.z = inv_delta[2] * E1.z + inv_delta[3] * E2.z; 

    T = parser::normalize(T);
    B = parser::normalize(B);


    parser::Vec3f new_normal; 
    parser::Vec3f texture_normal_from_color;
    parser::Vec3f ones;
    ones.x = 1;
    ones.y = 1; 
    ones.z = 1; 

    texture_normal_from_color = parser::normalize(textureColor / 127.5f - ones ); 
    parser::Vec3f N = parser::cross(B ,T);
    N = parser::normalize(N );
    new_normal.x = T.x * texture_normal_from_color.x + B.x * texture_normal_from_color.y + N.x * texture_normal_from_color.z;
    new_normal.y = T.y * texture_normal_from_color.x + B.y * texture_normal_from_color.y + N.y * texture_normal_from_color.z;
    new_normal.z = T.z * texture_normal_from_color.x + B.z * texture_normal_from_color.y + N.z * texture_normal_from_color.z;

    return parser::normalize(new_normal);


}

parser::Vec3f replace_normal_sphere_bump(parser::Sphere * sphere , parser::Vec3f & hitpoint ,  parser::Vec3f & textureColor  , float omega , float phi , parser::TextureMap & texturemap ,  float u , float v  )
{
    //get x y z again
    parser::Vec3f T;
    T.x = 2 * M_PI * hitpoint.z;
    T.y = 0;
    T.z = -2 * M_PI * hitpoint.x;
    T = parser::normalize(T);

    parser::Vec3f B;
    B.x = M_PI * std::cos(  phi * M_PI/180) * hitpoint.y;
    B.y = -sphere->radius * M_PI * std::sin(omega * M_PI / 180);
    B.z = M_PI * hitpoint.y * std::sin(phi * M_PI/180);
    B = parser::normalize(B);

 
    parser::Vec3f N = parser::cross(B ,T);
    N = parser::normalize(N );

    parser::Image * img =  texturemap.image; 
    int pix_x = std::floor(u * img->w);
    int pix_y = std::floor(v * img->h);
    parser::Vec3f color1;
    color1.x = img->image[ 3 * ( pix_y * img->w + pix_x)  ];
    color1.y = img->image[ 3 * ( pix_y * img->w + pix_x) + 1  ];
    color1.z = img->image[ 3 * ( pix_y * img->w + pix_x) + 2  ];
    // dif u 
    parser::Vec3f color2;
    color2.x = img->image[ 3 * ( pix_y * img->w + pix_x + 1)  ];
    color2.y = img->image[ 3 * ( pix_y * img->w + pix_x + 1) + 1  ];
    color2.z = img->image[ 3 * ( pix_y * img->w + pix_x +1) + 2  ];
    // diff v 
    parser::Vec3f color3;
    color3.x = img->image[ 3 * ( (pix_y + 1) * img->w + pix_x )  ];
    color3.y = img->image[ 3 * ( (pix_y + 1) * img->w + pix_x ) + 1  ];
    color3.z = img->image[ 3 * ( (pix_y + 1) * img->w + pix_x ) + 2  ];

    // W and H does not give good results !!!!!!!!!!!!!!!!!!!!!!!! 
    parser::Vec3f dif_u = (color2 - color1) / 100;
    parser::Vec3f dif_v = (color3 - color1) / 100;

    float average_dif_u = (dif_u.x + dif_u.y + dif_u.z) / 3;
    float average_dif_v = (dif_v.x + dif_v.y + dif_v.z) / 3;
    if( average_dif_v != 0 )
    {
        int a =1; 
    }
    parser::Vec3f q_u = T + N * average_dif_u  ;
    q_u = parser::normalize(q_u);
    parser::Vec3f q_v = B + N * average_dif_v ;
    q_v = parser::normalize(q_v);

     return parser::cross(q_v , q_u );

}
parser::Vec3f replace_normal_mesh(int* v1_texcoord ,int* v2_texcoord ,int* v3_texcoord , parser::Vec3f &v1 ,parser::Vec3f &v2,parser::Vec3f &v3, parser::Vec3f & textureColor  , parser::TextureMap & texturemap ,  float u , float v)
{
    float delta1[2];
    float delta2[2];
    //texture space 
    delta1[0] = v2_texcoord[0] - v1_texcoord[0];
    delta1[1] = v2_texcoord[1] - v1_texcoord[1];

    delta2[1] = v3_texcoord[1] - v2_texcoord[1]; 
    delta2[0] = v3_texcoord[0] - v2_texcoord[0];

    //world space 
    parser::Vec3f E1 = v2 - v1; 
    parser::Vec3f E2 = v3 - v2; 

    float inv_delta[4]; // row major a =0 b  = 1 c = 2 d = 3 
    float den = 1 / (delta1[0] * delta2[1] - delta1[1] * delta2[0]);
    inv_delta[0] = delta2[1] * den; 
    inv_delta[1] = - delta1[1] * den;
    inv_delta[2] = - delta2[0] * den;
    inv_delta[3] =  delta1[0] * den;


    // now matrix mul 
    parser::Vec3f T;
    T.x = inv_delta[0] * E1.x + inv_delta[1] * E2.x;
    T.y = inv_delta[0] * E1.y + inv_delta[1] * E2.y;
    T.z = inv_delta[0] * E1.z + inv_delta[1] * E2.z;

    parser::Vec3f B;
    B.x = inv_delta[2] * E1.x + inv_delta[3] * E2.x;
    B.y = inv_delta[2] * E1.y + inv_delta[3] * E2.y;
    B.z = inv_delta[2] * E1.z + inv_delta[3] * E2.z; 

    T = parser::normalize(T);
    B = parser::normalize(B);

    parser::Vec3f N = parser::cross(B ,T);
    N = parser::normalize(N );

 parser::Image * img =  texturemap.image; 
    int pix_x = std::floor(u * img->w);
    int pix_y = std::floor(v * img->h);
    parser::Vec3f color1;
    color1.x = img->image[ 3 * ( pix_y * img->w + pix_x)  ];
    color1.y = img->image[ 3 * ( pix_y * img->w + pix_x) + 1  ];
    color1.z = img->image[ 3 * ( pix_y * img->w + pix_x) + 2  ];
    // dif u 
    parser::Vec3f color2;
    color2.x = img->image[ 3 * ( pix_y * img->w + pix_x + 1)  ];
    color2.y = img->image[ 3 * ( pix_y * img->w + pix_x + 1) + 1  ];
    color2.z = img->image[ 3 * ( pix_y * img->w + pix_x +1) + 2  ];
    // diff v 
    parser::Vec3f color3;
    color3.x = img->image[ 3 * ( (pix_y + 1) * img->w + pix_x )  ];
    color3.y = img->image[ 3 * ( (pix_y + 1) * img->w + pix_x ) + 1  ];
    color3.z = img->image[ 3 * ( (pix_y + 1) * img->w + pix_x ) + 2  ];

    // W and H does not give good results !!!!!!!!!!!!!!!!!!!!!!!! 
    parser::Vec3f dif_u = (color2 - color1) / 100;
    parser::Vec3f dif_v = (color3 - color1) / 100;

    float average_dif_u = (dif_u.x + dif_u.y + dif_u.z) / 3;
    float average_dif_v = (dif_v.x + dif_v.y + dif_v.z) / 3;
    if( average_dif_v != 0 )
    {
        int a =1; 
    }
    parser::Vec3f q_u = T + N * average_dif_u  ;
    q_u = parser::normalize(q_u);
    parser::Vec3f q_v = B + N * average_dif_v ;
    q_v = parser::normalize(q_v);

    return parser::cross(q_v , q_u );

}
void init_perlin_noise()
{
    for (size_t i = 0; i < 256; i++)
    {
        permutation_table[i] = i; 
    }
    // now shuffle 
    std::random_shuffle(std::begin(permutation_table) , std::end(permutation_table));
}
parser::Vec3f hash_function(int index)
{
    parser::Vec3f vec; 
    int index_mod = index % 8;
    if( index_mod == 0 )
    {
        vec.x = 1.0f;
        vec.y = 0.0f;
        vec.z = 0.0f;

    }
    else if( index_mod == 1 )
    {
        vec.x = -1.0f;
        vec.y = 0.0f;
        vec.z = 0.0f;

    }
    else if( index_mod == 2 )
    {
        vec.x = 0.0f;
        vec.y = 1.0f;
        vec.z = 0.0f;
        

    }
    else if( index_mod == 3 )
    {
        vec.x = 0.0f;
        vec.y = -1.0f;
        vec.z = 0.0f;

    }
    else if( index_mod == 4 )
    {
        vec.x = 0.0f;
        vec.y = 0.0f;
        vec.z = 1.0f;

    }
    else if( index_mod == 5 )
    {
        vec.x = 0.0f;
        vec.y = 0.0f;
        vec.z = -1.0f;

    }
    else if( index_mod == 6 )
    {
        vec.x = 1.0f;
        vec.y = 0.0f;
        vec.z = 1.0f;
        vec = parser::normalize(vec);

    }
    else if( index_mod == 7 )
    {
        vec.x = 1.0f;
        vec.y = 1.0f;
        vec.z = 1.0f;
        vec = parser::normalize(vec);
    }
    return vec;

}