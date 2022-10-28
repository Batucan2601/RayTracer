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
    //Ray ray( glm::vec3( scene.cameras[0].position.x  , scene.cameras[0].position.y ,  scene.cameras[0].position.z)  , glm::vec3( scene.cameras[0].gaze.x , scene.cameras[0].gaze.y , scene.cameras[0].gaze.z ) );
    Ray r_test( glm::vec3( 0 , 0 , 0 ) , glm::vec3( 0 , 0 , -1 ) );
    glm::vec3 norm; 
    glm::vec3 hit_point ; 

    bool sphere = ray_sphere_intersection(r_test , scene.spheres[0] , glm::vec3(0 ,  0 ,  -2) , norm , hit_point   );
    


    //glm::vec3 p1(0.0f , 0.0f , 0.0f );
    
    //std::cout << calculate_intersection( scene , o , ray , p1   ) << std::endl; 
    //std::cout << p1.x << " " << p1.y << " " << p1.z << std::endl; 
    return 1; 
}