#pragma once 
#include "parser.h"
#include "tinyxml2.h"
#include <sstream>
#include <stdexcept>
#include <iostream>
void parser::Scene::loadFromXml(const std::string &filepath)
{
    tinyxml2::XMLDocument file;
    std::stringstream stream;

    auto res = file.LoadFile(filepath.c_str());
    if (res)
    {
        throw std::runtime_error("Error: The xml file cannot be loaded.");
    }

    auto root = file.FirstChild();
    if (!root)
    {
        throw std::runtime_error("Error: Root is not found.");
    }

    //Get BackgroundColor
    auto element = root->FirstChildElement("BackgroundColor");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0 0 0" << std::endl;
    }
    stream >> background_color.x >> background_color.y >> background_color.z;

    //Get ShadowRayEpsilon
    element = root->FirstChildElement("ShadowRayEpsilon");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0.001" << std::endl;
    }
    stream >> shadow_ray_epsilon;

    //Get MaxRecursionDepth
    element = root->FirstChildElement("MaxRecursionDepth");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0" << std::endl;
    }
    stream >> max_recursion_depth;

    //Get Cameras
    element = root->FirstChildElement("Cameras");
    element = element->FirstChildElement("Camera");
    Camera camera;
    std::cout << "1 " << std::endl;  
    while (element)
    {
        auto child = element->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Gaze");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Up");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("NearPlane");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("NearDistance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("ImageResolution");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("ImageName");
        stream << child->GetText() << std::endl;

        stream >> camera.position.x >> camera.position.y >> camera.position.z;
        stream >> camera.gaze.x >> camera.gaze.y >> camera.gaze.z;
        stream >> camera.up.x >> camera.up.y >> camera.up.z;
        stream >> camera.near_plane.x >> camera.near_plane.y >> camera.near_plane.z >> camera.near_plane.w;
        stream >> camera.near_distance;
        stream >> camera.image_width >> camera.image_height;
        stream >> camera.image_name;

        cameras.push_back(camera);
        element = element->NextSiblingElement("Camera");
    }

    //Get Lights
    element = root->FirstChildElement("Lights");
    auto child = element->FirstChildElement("AmbientLight");
    stream << child->GetText() << std::endl;
    stream >> ambient_light.x >> ambient_light.y >> ambient_light.z;
    element = element->FirstChildElement("PointLight");
    PointLight point_light;

    std::cout << "2 " << std::endl;  

    while (element)
    {
        child = element->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Intensity");
        stream << child->GetText() << std::endl;

        stream >> point_light.position.x >> point_light.position.y >> point_light.position.z;
        stream >> point_light.intensity.x >> point_light.intensity.y >> point_light.intensity.z;

        point_lights.push_back(point_light);
        element = element->NextSiblingElement("PointLight");
    }

    //Get Materials
    element = root->FirstChildElement("Materials");
    element = element->FirstChildElement("Material");
    Material material;


    std::cout << "3 " << std::endl;  

    while (element)
    {
        material.is_mirror = (element->Attribute("type", "mirror") != NULL);
        material.is_dielectric = (element->Attribute("type", "dielectric") != NULL);
        material.is_conductor = (element->Attribute("type", "conductor") != NULL);
        child = element->FirstChildElement("AmbientReflectance");
        if( child != NULL)
        {
            std::cout << " ambient ok " << std::endl; 
            stream << child->GetText() << std::endl;
            stream >> material.ambient.x >> material.ambient.y >> material.ambient.z;
        }
        else
        {
            std::cout << " ambient not ok " << std::endl; 

            material.ambient.x = 0.0f; 
            material.ambient.y = 0.0f; 
            material.ambient.z = 0.0f; 

        }
        child = element->FirstChildElement("DiffuseReflectance");
        if( child != NULL)
        {
            stream << child->GetText() << std::endl;
            stream >> material.diffuse.x >> material.diffuse.y >> material.diffuse.z;
        }
        else
        {
            material.diffuse.x = 0.0f; 
            material.diffuse.y = 0.0f; 
            material.diffuse.z = 0.0f; 

        }
        child = element->FirstChildElement("SpecularReflectance");
        if( child != NULL)
        {
        stream << child->GetText() << std::endl;
        stream >> material.specular.x >> material.specular.y >> material.specular.z;

        }
        else
        {
            material.specular.x = 0.0f;
            material.specular.y = 0.0f;
            material.specular.z = 0.0f;

        }
        child = element->FirstChildElement("MirrorReflectance");
        if( child != NULL)
        {
            stream << child->GetText() << std::endl;
            stream >> material.mirror.x >> material.mirror.y >> material.mirror.z;
        }
        else
        {
            material.mirror.x =0.0f;
            material.mirror.y =0.0f;
            material.mirror.z =0.0f;

        }
        child = element->FirstChildElement("PhongExponent");
        if( child != NULL)
        {
            stream << child->GetText() << std::endl;
            stream >> material.phong_exponent;
        }
        else
        {
            material.phong_exponent = 0.0f; 
        }
        child = element->FirstChildElement("AbsorptionCoefficient");
        if( child != NULL)
        {
            stream << child->GetText() << std::endl;
            stream >> material.absorptionCoefficient.x >> material.absorptionCoefficient.y >> material.absorptionCoefficient.z ;
        }
        else
        {
            material.absorptionCoefficient.x = 0.0f;
            material.absorptionCoefficient.y = 0.0f;
            material.absorptionCoefficient.z = 0.0f;

        }
        child = element->FirstChildElement("RefractionIndex");
        if( child != NULL)
        {
            stream << child->GetText() << std::endl;
            stream >> material.refraction_index;
        }
        else
        {
            material.refraction_index =  0.0f; 
        }
        materials.push_back(material);
        element = element->NextSiblingElement("Material");
    }
    std::cout << "4 " << std::endl;

    //Get VertexData
    element = root->FirstChildElement("VertexData");
    stream << element->GetText() << std::endl;
    std::cout << element->GetText();
    Vec3f vertex;
    while (!(stream >> vertex.x).eof() )
    {
         
        stream >> vertex.y >> vertex.z;
        vertex_data.push_back(vertex);
        std::cout << vertex.x << " " << vertex.y << " " << vertex.z  << std::endl;  

    }

    stream.clear();

    //Get Meshes
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Mesh");
    Mesh mesh;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> mesh.material_id;

        child = element->FirstChildElement("Faces");
        stream << child->GetText() << std::endl;
        Face face;
        while (!(stream >> face.v0_id).eof())
        {
            stream >> face.v1_id >> face.v2_id;
            mesh.faces.push_back(face);
        }
        stream.clear();

        meshes.push_back(mesh);
        mesh.faces.clear();
        element = element->NextSiblingElement("Mesh");
    }
    std::cout << "6 " << std::endl;  

    stream.clear();

    //Get Triangles
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Triangle");
    Triangle triangle;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> triangle.material_id;

        child = element->FirstChildElement("Indices");
        stream << child->GetText() << std::endl;
        stream >> triangle.indices.v0_id >> triangle.indices.v1_id >> triangle.indices.v2_id;

        triangles.push_back(triangle);
        element = element->NextSiblingElement("Triangle");
    }
    //Get Spheres
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Sphere");
    Sphere sphere;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> sphere.material_id;

        child = element->FirstChildElement("Center");
        stream << child->GetText() << std::endl;
        stream >> sphere.center_vertex_id;

        child = element->FirstChildElement("Radius");
        stream << child->GetText() << std::endl;
        stream >> sphere.radius;

        spheres.push_back(sphere);
        element = element->NextSiblingElement("Sphere");
    }
}
// geometry functions
float parser::dot( const Vec3f& vec1 , const Vec3f& vec2   )
{
    return vec1.x * vec2.x +  vec1.y * vec2.y + vec1.z * vec2.z;
}

parser::Vec3f parser::cross( const parser::Vec3f& vec1 , const parser::Vec3f & vec2 )
{
    parser::Vec3f cross; 
    cross.x = (vec1.y * vec2.z) - (vec1.z * vec2.y);
    cross.y = -1 * (vec1.x * vec2.z) - (vec1.z * vec2.x);
    cross.z = (vec1.x * vec2.y) - (vec1.y * vec2.x);

    return cross;
}
float parser::length(const parser::Vec3f & vec1 )
{
    return std::sqrt( vec1.x*vec1.x + vec1.y*vec1.y + vec1.z*vec1.z );
}
parser::Vec3f parser::normalize( const  parser::Vec3f &  vec1 )
{
    parser::Vec3f new_vec; 
    float len = length(vec1);
    new_vec.x =  vec1.x / len;
    new_vec.y =  vec1.y / len;
    new_vec.z =  vec1.z / len;
    return new_vec;

}
parser::Vec3f::Vec3f(float x  , float y , float z )
{
    this->x = x;
    this->y = y;
    this->z = z;

}
parser::Vec3f::Vec3f()
{
    this->x = 0;
    this->y = 0;
    this->z = 0;

}