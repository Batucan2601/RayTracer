#pragma once
#include "../System_Files/support_files/ppm.h"
#include "../System_Files/support_files/parser.h"
#include "../Dependencies/glm_headers/glm.hpp"
static void generate_image(parser::Scene & scene);
static void calculate_image_plane(const parser::Camera &current_camera , parser::Vec3f & starting_point,  parser::Vec3f & interwal_row, parser::Vec3f & interwal_col);
static bool ray_plane_intersection(const Ray& ray , const glm::vec3 & p1 , const glm::vec3 & normal , glm::vec3 & hitpoint  );
static bool ray_sphere_intersection(const Ray& ray , const parser::Sphere& sphere , const glm::vec3 & center , glm::vec3 &normal ,  glm::vec3 & hitpoint  );
static glm::vec3 color_pixel(parser::Scene& scene , Ray & ray );
static bool calculate_intersection(parser::Scene& scene ,parser::Mesh& object , const Ray & ray ,glm::vec3 & intersection_normal ,   glm::vec3 & intersection_point );
static bool calculate_intersection(parser::Scene& scene ,parser::Sphere& object , const Ray & ray , glm::vec3 & intersection_normal ,   glm::vec3 & intersection_point );
static bool ray_triangle_intersection(const Ray& ray , const glm::vec3 & p1 ,const  glm::vec3 & p2 , const glm::vec3 &  p3 , glm::vec3 & hitpoint );
static bool is_point_in_triangle(const glm::vec3 & p1 ,const glm::vec3 & p2 , const glm::vec3 & p3 , const glm::vec3 & hitpoint);
static bool ray_object_intersection( const Ray & ray , parser::Scene & scene , glm::vec3 &hitpoint , glm::vec3 & normal , parser::Material &material ,  int &  prev_object_id  , bool is_shadow_rays_active  );
static bool calculate_intersection(parser::Scene& scene ,parser::Triangle& object ,const Ray & ray , glm::vec3 & intersection_normal ,  glm::vec3 & intersection_point  );


