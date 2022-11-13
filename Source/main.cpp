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
    std::cout << " LOADING  DONE " << std::endl;
    generate_image(scene);
    std::cout << " IMAGE GENERATION DONE " << std::endl;
    
   

    //std::vector< std::pair<BVH::BoundingBox, int > > vec = BVH::generate_bounding_boxes(scene);
    //std::cout << vec.size()<<  " intersection  " << std::endl;;
    //std::cout << p1.x << " " << p1.y << " " << p1.z << std::endl; 
    return 1; 
}