#include "Intersections.h"
#include <math.h>       /* isnan, std */
#include "../System_Files/support_files/parser.h"

bool is_texture_present( parser::Scene & scene ,  int object_id , parser::Vec3f & intersection_point  , parser::Vec3f & intersection_normal , parser::Vec3f &intersection_face , parser::Vec3f & texture_color )
{
    if( object_id < scene.meshes.size())
    {
        //object is mesh 
        if( scene.meshes[object_id].textures.size() > 0 )
        {
            get_texture_color_from_mesh( scene ,  object_id ,  intersection_point  ,  intersection_normal , intersection_face ,  texture_color);
            return true; 
        }
        return false; 
    }
    else if(object_id < scene.meshes.size() + scene.spheres.size() )
    {
        //object is sphere 
        if( scene.spheres[object_id - scene.meshes.size()].textures.size() > 0 )
        {
            get_texture_color_from_sphere(  scene ,   object_id ,  intersection_point  ,  intersection_normal  ,  texture_color );
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
parser::Vec3f get_texture_color_from_mesh( parser::Scene & scene ,  int object_id , parser::Vec3f & intersection_point  , parser::Vec3f & intersection_normal , parser::Vec3f &intersection_face , parser::Vec3f &  texture_color )
{
    // here goes the texture getting stuff
    parser:Mesh *mesh =  &scene.meshes[object_id];
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
        

    }
    

}
parser::Vec3f get_texture_color_from_sphere( parser::Scene & scene ,  int object_id , parser::Vec3f & intersection_point  , parser::Vec3f & intersection_normal  , parser::Vec3f &  texture_color )
{

}