#ifndef __HW1__PARSER__
#define __HW1__PARSER__
#include <string>
#include <vector>
#include <cmath> 
#include "tinyxml2.h"
#include <sstream>
#include <stdexcept>
#include <iostream>
#include "../support_files/happly.h"
#include <unistd.h>
namespace parser
{
    //Notice that all the structures are as simple as possible
    //so that you are not enforced to adopt any style or design.
    class Vec3f
    {
        public:
        Vec3f();
        Vec3f( float x , float y , float z );
        Vec3f operator +(const Vec3f & vec1   );
        Vec3f operator -(const Vec3f & vec1   );
        Vec3f operator *(const Vec3f & vec1   );
        Vec3f operator /(const Vec3f & vec1   );
        void operator =(const Vec3f & vec1 );
        Vec3f operator *(const int & num   );
        Vec3f operator +(const int & num   );
        Vec3f operator -(const int & num   );
        Vec3f operator /(const int & num   );
        Vec3f operator *(const float & num   );
        Vec3f operator +(const float & num   );
        Vec3f operator -(const float & num   );
        Vec3f operator /(const float & num   );



        //members
        float x, y, z;
    };

    struct Vec3i
    {
        int x, y, z;
    };

    struct Vec4f
    {
        float x, y, z, w;
    };
    
    class Matrix
    {
        public:
        Matrix( int row , int col );
        Matrix();

        float at(int row_index ,int col_index );
        void set(int row , int col , float num);
        void Identity();
        void Transpose();
        Matrix inverse();
        //operator
        Vec4f operator *(const Vec4f & vec1   ); // 4x4 matrix vs 4x1 vector 
        Matrix operator*( Matrix & matrix   ); // 4x4 matrix vs 4x1 vector 

        std::vector<float> elements;
        int col_no; 
        int row_no; 
    }; 
    struct Camera
    {
        Vec3f position;
        Vec3f gaze;
        Vec3f up;
        Vec4f near_plane;
        float near_distance;
        int image_width, image_height;
        std::string image_name;
        int number_of_samples; 
        float focus_distance;
        float aperture_size;

    };

    struct PointLight
    {
        Vec3f position;
        Vec3f intensity;
    };
    struct OrthonormalBasis
    {
        Vec3f normal; 
        Vec3f u;
        Vec3f v;

    };
    struct AreaLight
    {
        Vec3f position;
        Vec3f normal;
        Vec3f radiance;
        int size; 
        //orthonormal bases 
        OrthonormalBasis ortho_basis; 
    };
    struct Image
    {
        std::string path;
        unsigned char* image; 
        int w, h, comp;
    };
    struct TextureMap
    {
        Image* image;
        std::string decalMode; 
        std::string interpolation;  
        bool is_image; // no if perlin yes if texture
        std::string noiseConversion; // no if perlin yes if texture 
        std::string noiseScale; // no if perlin yes if texture 

    };
    struct Material
    {
        bool is_mirror;
        bool is_dielectric;
        bool is_conductor;
        Vec3f ambient;
        Vec3f diffuse;
        Vec3f specular;
        Vec3f mirror;
        Vec3f absorptionCoefficient; 
        float phong_exponent;
        float refraction_index;
        float absorption_index; 
        float roughness;
    };

    struct Face
    {
        int v0_id;
        int v1_id;
        int v2_id;
    };

    struct Transformations
    {
        std::vector<Matrix> scaling;
        std::vector<Matrix> translation;
        std::vector<Matrix> rotation; 
    };
    struct Mesh
    {
        int material_id;
        std::vector<Face> faces;
        std::vector<Matrix> transformations;
        Vec3f motion_blur; 
        Vec3f current_motion_blur;  // for shadow rays 
        Matrix transformation; 
        Matrix transformation_inverse; 
        std::vector<int> textures; 

    };
    struct MeshInstance
    {
        Mesh* mesh_ptr; 
        bool reset_transform; 
        std::vector<Matrix> transformations;
        int material_id;
        Vec3f motion_blur;
        Vec3f current_motion_blur; 

        Matrix transformation; 
        Matrix transformation_inverse; 
    };
    struct Triangle
    {
        int material_id;
        Face indices;
        std::vector<Matrix> transformations;

        Matrix transformation; 
        Matrix transformation_inverse; 
    };

    struct Sphere
    {
        int material_id;
        int center_vertex_id;
        float radius;
        Vec3f motion_blur; 
        Vec3f current_motion_blur; 

        std::vector<Matrix> transformations;
        std::vector<int> scale_indices; 

        Matrix transformation; 
        Matrix transformation_inverse; 
        std::vector<int> textures; 

    };

    struct Scene
    {
        //Data
        Vec3i background_color;
        float shadow_ray_epsilon;
        int max_recursion_depth;
        std::vector<Camera> cameras;
        Vec3f ambient_light;
        std::vector<PointLight> point_lights;
        std::vector<AreaLight> area_lights;
        std::vector<Image> images;
        std::vector<TextureMap> textureMaps;
        std::vector<std::array<int , 2U > > texcoords; 


        std::vector<Material> materials;
        std::vector<Vec3f> vertex_data;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;
        std::vector<MeshInstance> mesh_instances;
        Transformations transformations;
        //Functions
        void loadFromXml(const std::string &filepath);

        
    };


    float dot( const Vec3f& vec1 , const Vec3f& vec2   );
    Vec3f cross( const Vec3f& vec1 , const Vec3f & vec2 );
    float length(const Vec3f & vec1 );
    Vec3f normalize( const  Vec3f &  vec1 );
    float distance(const  Vec3f &  vec1 ,  const  Vec3f & vec2   );

    
    void print_matrix(const Matrix & mat );
}

#endif
