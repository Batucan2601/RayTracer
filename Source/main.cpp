#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#define DEBUG
#include "../Header/Ray.h"
#include "../System_Files/support_files/ppm.h"
#include "../System_Files/support_files/parser.h"
parser::Scene scene;
parser::Camera current_camera;
#include "../Header/Intersections.h"




int main(int argc, char* argv[])
{
    
    std::cout << argv[1] << std::endl; 
    scene.loadFromXml(argv[1]); // load the scene
    std::cout << "loaded succesful y " << std::endl;
    generate_image(scene);
    Ray ray( parser::Vec3f( scene.cameras[0].position.x  , scene.cameras[0].position.y ,  scene.cameras[0].position.z)  , parser::Vec3f(-0.12f, 0.2f ,  0.5f) );
    parser::Vec3f p1(-0.0303071, 0.0339954,  0.0147064);
    parser::Vec3f p2(-0.0348346, 0.0347214, 0.0345782);
    parser::Vec3f p3(-0.0352833, 0.0341798, 0.0227449);
    parser::Vec3f hit; 
    parser::Vec3f normal; 
    bool is_hit = calculate_intersection(scene , scene.meshes[0] , ray , hit , normal  );
    std::cout <<  " is hit " << is_hit << std::endl;
    std::cout << "hitpoint "  << std::endl; 
    
    //glm::vec3 p1(0.0f , 0.0f , 0.0f );
    
    //std::cout << calculate_intersection( scene , o , ray , p1   ) << std::endl; 
    //std::cout << p1.x << " " << p1.y << " " << p1.z << std::endl; 
    return 1; 
}