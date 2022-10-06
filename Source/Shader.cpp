#include "../Header/Shader.h"
#include <GL/glew.h>
#include <GL/glut.h>
#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

Shader::Shader(std::string filepath ){
    filepath_location = filepath;
    id = 0;
    ShaderProgramSource source = parseShader(filepath); // shader.shader to fix 
    
    id = CreateShaders(source.Vertex , source.Fragment);
    //std::cout << source.Vertex << "\n" << source.Fragment << std::endl;
}
Shader::~Shader(){
    glDeleteProgram(id);

}
unsigned int Shader::CreateShaders(const std::string& vertexShader , const std::string& fragmentShader ){
   unsigned int program =  glCreateProgram();
   unsigned int vs = CompileShader(GL_VERTEX_SHADER , vertexShader);
   unsigned int fs = CompileShader(GL_FRAGMENT_SHADER , fragmentShader);
   
   glAttachShader(program , vs );
   glAttachShader(program , fs );

   glLinkProgram(program);
   glValidateProgram(program);
   glDeleteShader(vs);
   glDeleteShader(fs);

    return program; 
}
 unsigned int Shader::CompileShader( const unsigned int type  , const std::string source  ){
    unsigned int id = glCreateShader(type);
    const char* src = source.c_str();
    glShaderSource(id,1,&src , nullptr);
    glCompileShader(id);

    int result; 
    glGetShaderiv(id , GL_COMPILE_STATUS , &result);
    if( result == GL_FALSE ){
        int length;
        glGetShaderiv(id , GL_INFO_LOG_LENGTH, &length );
        char* message = (char*)alloca(length *sizeof(char)  );
        glGetShaderInfoLog(id,length , &length  , message);
    }
    return id;
}
ShaderProgramSource Shader::parseShader(const std::string & filepath ){
    std::ifstream stream(filepath);
    std::cout << filepath << std::endl;
    enum class ShaderType{
        NONE = -1, VERTEX = 0, FRAGMENT = 1
    };
    std::string line; 
    std::stringstream ss[2];
    ShaderType type = ShaderType::NONE;
    while(getline(stream,line) ){
        if( line.find("#shader") != std::string::npos ){
            if( line.find("vertex") != std::string::npos ){
                type = ShaderType::VERTEX;
            }
            else if( line.find("fragment") != std::string::npos ){
                type = ShaderType::FRAGMENT;
            }
        }
        else{
            ss[(int)type] << line << '\n';
        }
    }
    return { ss[0].str() , ss[1].str()};
}
void Shader::bind(){
    glUseProgram(id);
}

void Shader::unbind(){
    glUseProgram(0);
}
unsigned int Shader::GetUniformLocation(const std::string&name ){
    unsigned int location = glGetUniformLocation(id , name.c_str());
    return location;
}
void Shader::setUniform4f( const std::string &name , float v0 , float v1, float v2 , float v3){
     glUniform4f(GetUniformLocation(name) , v0 , v1 , v2 , v3 );
}

std::string  Shader::getFilePath(){
    return filepath_location;
}

void Shader::setUniform1I(const std::string& name , float value ){
glUniform1f( GetUniformLocation( name ) , value );
}

void Shader::setUniformMat4f(const std::string& name , const glm::mat4& matrix  ){
    glUniformMatrix4fv(GetUniformLocation(name) , 1, GL_FALSE , &matrix[0][0] );
    
}

void Shader::setUniformVec3(const std::string&name, const glm::vec3& vector ){
    glUniform3fv(GetUniformLocation(name ) ,  1,  &vector[0]);
}

void Shader::setUniformVec4(const std::string&name, const glm::vec4& vector ){
    glUniform4fv(GetUniformLocation(name ) ,  1,  &vector[0]);
}

unsigned int Shader::getid(){
    return id;
}