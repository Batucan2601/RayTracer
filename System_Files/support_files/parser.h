#ifndef __HW1__PARSER__
#define __HW1__PARSER__
#include <string>
#include <vector>
#include <cmath> 
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
        float at(int row_index ,int col_index );
        void set(int row , int col , float num);
        void Identity();

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
    };

    struct PointLight
    {
        Vec3f position;
        Vec3f intensity;
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
    };
    struct MeshInstance
    {
        Mesh* mesh_ptr; 
        bool reset_transform; 
        std::vector<Matrix> transformations;
        int material_id;
    };
    struct Triangle
    {
        int material_id;
        Face indices;
        std::vector<Matrix> transformations;

    };

    struct Sphere
    {
        int material_id;
        int center_vertex_id;
        float radius;
        std::vector<Matrix> transformations;
        int scale_index; 
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
