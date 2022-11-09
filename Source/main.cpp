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
#include "../Header/BVH.h"



int main(int argc, char* argv[])
{
    
    std::cout << argv[1] << std::endl; 
    scene.loadFromXml(argv[1]); // load the scene
    std::cout << " LOADING  DONE " << std::endl;
    generate_image(scene);
    std::cout << " IMAGE GENERATION DONE " << std::endl;
    BVH::BoundingBox b(-10.0f , 10.0f , 11.0f , -10.0f , 12.0f , 20.0f  );
    parser::Vec3f O(0.0f,0.0f , 0.0f );
    parser::Vec3f d(0.1f,0.0f , 0.9f ); 

    Ray r(  O, d );
    //std::vector< std::pair<BVH::BoundingBox, int > > vec = BVH::generate_bounding_boxes(scene);
    //std::cout << vec.size()<<  " intersection  " << std::endl;;
    //std::cout << p1.x << " " << p1.y << " " << p1.z << std::endl; 
    return 1; 
}