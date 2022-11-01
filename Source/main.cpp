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
    generate_image(scene);
    std::cout << "loaded succesful y " << std::endl;
    parser::Matrix m1(4 ,4 );
    parser::print_matrix(m1);
    std::cout << " //////////////////" << std::endl; 
    m1.Identity();
    m1.set(1, 3 , 10 );
    parser::print_matrix(m1);
    std::cout << m1.at(0 , 0 ) <<  " " << m1.at(1,3);
    //std::cout << calculate_intersection( scene , o , ray , p1   ) << std::endl; 
    //std::cout << p1.x << " " << p1.y << " " << p1.z << std::endl; 
    return 1; 
}