#pragma once
#include "../System_Files/support_files/ppm.h"
#include "../System_Files/support_files/parser.h"
#include "../Header/Ray.h"
#include "../Header/BVH.h"
static void generate_image(parser::Scene & scene);
static void calculate_image_plane(const parser::Camera &current_camera , parser::Vec3f & starting_point,  parser::Vec3f & interwal_row, parser::Vec3f & interwal_col);
static bool ray_plane_intersection(const Ray& ray , const parser::Vec3f & p1 , const parser::Vec3f & normal , parser::Vec3f & hitpoint  );
static bool ray_sphere_intersection(const Ray& ray ,  parser::Sphere& sphere , const parser::Vec3f & center , parser::Vec3f &normal ,  parser::Vec3f & hitpoint  );
static parser::Vec3f color_pixel(parser::Scene& scene , Ray & ray );
static bool calculate_intersection(parser::Scene& scene ,parser::Mesh& object ,const Ray & ray , parser::Vec3f & intersection_normal ,  parser::Vec3f & intersection_point ,parser::Face &hit_face,  bool is_shadow_rays_active );
static bool calculate_intersection(parser::Scene& scene ,parser::MeshInstance& object , const Ray & ray ,parser::Vec3f & intersection_normal ,   parser::Vec3f & intersection_point , parser::Face &hit_face , bool is_shadow_rays_active  );
static bool calculate_intersection(parser::Scene& scene ,parser::Sphere& object , const Ray & ray , parser::Vec3f & intersection_normal ,   parser::Vec3f & intersection_point , bool is_shadow_rays_active  );
static bool ray_triangle_intersection(const Ray& ray , const parser::Vec3f & p1 ,const  parser::Vec3f & p2 , const parser::Vec3f &  p3 , parser::Vec3f & hitpoint );
static bool is_point_in_triangle(const parser::Vec3f & p1 ,const parser::Vec3f & p2 , const parser::Vec3f & p3 , const parser::Vec3f & hitpoint);
static bool ray_object_intersection( const Ray & ray , parser::Scene & scene , parser::Vec3f &hitpoint , parser::Vec3f & normal , parser::Material &material ,  int &  prev_object_id  , parser::Face & hit_face , bool is_shadow_rays_active  );
static bool calculate_intersection(parser::Scene& scene ,parser::Triangle& object ,const Ray & ray , parser::Vec3f & intersection_normal ,  parser::Vec3f & intersection_point  );
static bool calculate_second_hitpoint_in_same_object( parser::Scene & scene , const Ray & refracted_ray ,  parser::Vec3f & hit_point , parser::Vec3f & normal  , int & object_id , parser::Vec3f & second_hit_point  , parser::Vec3f & second_normal ,parser::Face & hit_face, bool is_shadow_rays_active  );
//static bool ray_object_intersection( const Ray & ray , parser::Scene & scene , parser::Vec3f &hitpoint , parser::Vec3f & normal , parser::Material &material ,   int &  prev_object_id  , bool is_shadow_rays_active, int & object_id  );

//texture operations
bool is_texture_present( parser::Scene & scene ,  int object_id , parser::Vec3f & intersection_point  , parser::Vec3f & intersection_normal , parser::Face &intersection_face , parser::Vec3f & texture_color ,  parser::Material &modified_material );
parser::Vec3f get_texture_color_from_mesh( parser::Scene & scene ,  int object_id , parser::Vec3f & intersection_point  , parser::Vec3f & intersection_normal , parser::Face &intersection_face , parser::Vec3f &  texture_color , parser::Material & modified_material );
parser::Vec3f get_texture_color_from_sphere( parser::Scene & scene ,  int object_id , parser::Vec3f & intersection_point  , parser::Vec3f & intersection_normal  , parser::Vec3f &  texture_color ,  parser::Material & new_material );
parser::Vec3f get_bilinear_coord_color(  parser::TextureMap & texturemap , float u , float v );
parser::Vec3f get_nearest_coord_color(  parser::TextureMap & texturemap , float u , float v );
parser::Vec3f replace_normal_sphere(parser::Sphere * sphere , parser::Vec3f & hitpoint ,  parser::Vec3f & textureColor  , float omega , float phi );
parser::Vec3f replace_normal_mesh(int* v1_texcoord ,int* v2_texcoord ,int* v3_texcoord , parser::Vec3f &v1 ,parser::Vec3f &v2,parser::Vec3f &v3, parser::Vec3f & textureColor );
parser::Vec3f replace_normal_sphere_bump(parser::Sphere * sphere , parser::Vec3f & hitpoint ,  parser::Vec3f & textureColor  , float omega , float phi , parser::TextureMap & texturemap ,  float u , float v  );
parser::Vec3f replace_normal_mesh(int* v1_texcoord ,int* v2_texcoord ,int* v3_texcoord , parser::Vec3f &v1 ,parser::Vec3f &v2,parser::Vec3f &v3, parser::Vec3f & textureColor  , parser::TextureMap & texturemap ,  float u , float v);
void init_perlin_noise();
parser::Vec2f hash_function(int index);
parser::Vec3f get_perlin_color_sphere(parser::Sphere * sphere , parser::TextureMap & texturemap ,  parser::Vec3f hit_point   );
parser::Vec3f replace_normal_sphere_bump_perlin(parser::Sphere * sphere , parser::TextureMap & texturemap ,  parser::Vec3f hit_point , parser::Vec3f  hit_normal);
parser::Vec3f parse_background(  const parser::Camera &current_camera , Ray  ray,  parser::Vec3f  interval_row, parser::Vec3f  interval_col  );
parser::Vec3f get_perlin_color_mesh( parser::TextureMap & texturemap ,  parser::Vec3f hit_point   );