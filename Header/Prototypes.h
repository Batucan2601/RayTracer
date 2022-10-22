#pragma once
#include "../System_Files/support_files/ppm.h"
#include "../System_Files/support_files/parser.h"
static void generate_image(parser::Scene & scene);
static void calculate_image_plane(const parser::Camera &current_camera , parser::Vec3f & starting_point,  parser::Vec3f & interwal_row, parser::Vec3f & interwal_col);
static bool ray_plane_intersection(const Ray& ray , const Plane& plane , glm::vec3 & hitpoint  );
static bool ray_sphere_intersection(const Ray& ray , const Sphere& sphere , glm::vec3 & hitpoint  );
static float color_pixel(parser::Scene& scene , Ray & ray );
static bool calculate_intersection(parser::Scene& scene ,parser::Mesh& object ,  glm::vec3 & intersection_point );
static bool calculate_intersection(parser::Scene& scene ,parser::Sphere& object ,  glm::vec3 & intersection_point );
