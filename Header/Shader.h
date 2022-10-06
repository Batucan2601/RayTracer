#pragma once 
#include <string>
#include "../Dependencies/glm_headers/glm.hpp"
struct ShaderProgramSource{
            std::string Vertex; 
            std::string Fragment; 
        };
class Shader{
    private:
        
        std::string filepath_location; 
        unsigned int id;
    public:
        Shader(  std::string filepath );
        ~Shader();
        
        void bind();
        void unbind();
        bool compileShader();
        std::string getFilePath();
        unsigned int getid();
        //uniforms
        void setUniform4f(const std::string &name , float v0 , float v1, float f2 , float f3 );
        void setUniform1I( const std::string &name , float value );
        void setUniformMat4f(const std::string & name , const glm::mat4& matrix );
        void setUniformVec4(const std::string&name, const glm::vec4& vector );
        void setUniformVec3(const std::string&name, const glm::vec3& vector );
    private:
        unsigned int CreateShaders(const std::string& vertexShader , const std::string& fragmentShader );
        unsigned int GetUniformLocation(const std::string&name );
        unsigned int CompileShader( const unsigned int type  , const std::string source  );
        ShaderProgramSource  parseShader(const std::string & filepath );
};