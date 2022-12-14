#pragma once 
#include "parser.h"
#define STB_IMAGE_IMPLEMENTATION
#include "/home/batu/Desktop/RayTracer/RayTracer/Header/stb_image.h"
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
    while (element)
    {

        if( element->Attribute("type") )
        {

            auto child = element->FirstChildElement("Position");
            stream << child->GetText() << std::endl;
            child = element->FirstChildElement("GazePoint");
            stream << child->GetText() << std::endl;
            child = element->FirstChildElement("Up");
            stream << child->GetText() << std::endl;
            child = element->FirstChildElement("FovY");
            stream << child->GetText() << std::endl;
            child = element->FirstChildElement("NearDistance");
            stream << child->GetText() << std::endl;
            child = element->FirstChildElement("ImageResolution");
            stream << child->GetText() << std::endl;
            child = element->FirstChildElement("ImageName");
            stream << child->GetText() << std::endl;
            child = element->FirstChildElement("NumSamples"); //assume square 
            stream << child->GetText() << std::endl;
            
            Vec3f gazePoint;  
            float fovy; 
            stream >> camera.position.x >> camera.position.y >> camera.position.z;
            stream >> gazePoint.x >>  gazePoint.y >>  gazePoint.z;
            stream >> camera.up.x >> camera.up.y >> camera.up.z;
            stream >> fovy;
            stream >> camera.near_distance;
            stream >> camera.image_width >> camera.image_height;
            stream >> camera.image_name;
            stream >> camera.number_of_samples;
            float t = camera.near_distance *  std::tan( fovy * M_PI/180  / 2 );
            float bottom = -t;
            float aspect = camera.image_width /  camera.image_height;
            float b = -t; 
            float r = t * aspect; 
            float l = -r; 
            camera.near_plane.x = l; 
            camera.near_plane.y = r;
            camera.near_plane.z = b; 
            camera.near_plane.w = t; 

            camera.gaze = parser::normalize( gazePoint - camera.position);

            camera.focus_distance = -1;
            camera.aperture_size = -1 ;
        }   
        else
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
            child = element->FirstChildElement("NumSamples"); //assume square 
            if (child != NULL )
            {
            stream << child->GetText() << std::endl;
            }
            else
            {
                stream << "1" << std::endl;
            }
            child = element->FirstChildElement("FocusDistance"); //assume square 
            if (child != NULL )
            {
            stream << child->GetText() << std::endl;
            }
            else
            {
                stream << "-1" << std::endl;
            }
            child = element->FirstChildElement("ApertureSize"); //assume square 
            if (child != NULL )
            {
            stream << child->GetText() << std::endl;
            }
            else
            {
                stream << "-1" << std::endl;
            }
            stream >> camera.position.x >> camera.position.y >> camera.position.z;
            stream >> camera.gaze.x >> camera.gaze.y >> camera.gaze.z;
            stream >> camera.up.x >> camera.up.y >> camera.up.z;
            stream >> camera.near_plane.x >> camera.near_plane.y >> camera.near_plane.z >> camera.near_plane.w;
            stream >> camera.near_distance;
            stream >> camera.image_width >> camera.image_height;
            stream >> camera.image_name;
            stream >> camera.number_of_samples;
            stream >> camera.focus_distance;
            stream >> camera.aperture_size;


        }
        

        cameras.push_back(camera);
        element = element->NextSiblingElement("Camera");
    }
    std::cout << "2" << std::endl;  

    //Get Lights
    element = root->FirstChildElement("Lights");
    auto child = element->FirstChildElement("AmbientLight");
    stream << child->GetText() << std::endl;
    stream >> ambient_light.x >> ambient_light.y >> ambient_light.z;

    //point light 
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
    // area light
    element = root->FirstChildElement("Lights");
    element = element->FirstChildElement("AreaLight");
    AreaLight area_light;
    while( element )
    {
        child = element->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Normal");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Size");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Radiance");
        stream << child->GetText() << std::endl;
        stream >> area_light.position.x >> area_light.position.y >> area_light.position.z;
        stream >> area_light.normal.x >> area_light.normal.y >> area_light.normal.z;
        stream >> area_light.size;
        stream >> area_light.radiance.x >> area_light.radiance.y >> area_light.radiance.z;


        //calculate ortho basis
        area_light.ortho_basis.normal = area_light.normal; 
        // n'
        if( std::abs(area_light.normal.x) < std::abs(area_light.normal.y) && std::abs(area_light.normal.x) < std::abs(area_light.normal.z))
        {
            // x is smallest
            area_light.ortho_basis.normal.x =  1.0f; 
        }
        else if( std::abs(area_light.normal.y) < std::abs(area_light.normal.x) && std::abs(area_light.normal.y) < std::abs(area_light.normal.z))
        {
            // y is smallest
            area_light.ortho_basis.normal.y =  1.0f; 
        }
        else
        {
            // z is smallest
            area_light.ortho_basis.normal.z =  1.0f; 
        }
        //compute u
        area_light.ortho_basis.u = parser::normalize(parser::cross(area_light.ortho_basis.normal , area_light.normal ));
        area_light.ortho_basis.v = parser::normalize(parser::cross(area_light.normal , area_light.ortho_basis.u ));
        area_lights.push_back(area_light);

        element = element->NextSiblingElement("AreaLight");

    }
    //get textures 
    element = root->FirstChildElement("Textures");
    element = element->FirstChildElement("Images");
    if( element != NULL )
    {
        element = element->FirstChildElement("Image");
        while( element )
        {
            Image image; 
            char cwd[500];
            std::string current_path;
            getcwd(cwd , sizeof(cwd));
            current_path = cwd;  
            current_path = current_path + "/hw4/inputs/";
            current_path = current_path + element->GetText();
            int w , h , comp;
            image.path = current_path;   
            image.image = stbi_load(current_path.c_str(), &w, &h, &comp, STBI_rgb);
            image.w = w;
            image.h = h;
            image.comp = comp;
            images.push_back(image);   
            element = element->NextSiblingElement("Image");
        }
    }
    
    // get texture maps 
    element = root->FirstChildElement("Textures");
    element = element->FirstChildElement("TextureMap");
    if( element != NULL )
    {
        while(element)
        {
            TextureMap textureMap; 
            if( element->Attribute("type" , "image") != NULL  )
            {
                textureMap.is_image = true;
                child =  element->FirstChildElement("ImageId");
                int image_id = std::stoi(child->GetText());
                textureMap.image = &images[image_id - 1 ];
                child =  element->FirstChildElement("DecalMode");
                std::string decalMode = child->GetText();
                textureMap.decalMode = decalMode;
                if( decalMode == "replace_background")
                {
                    //is_texture_background = true; 
                    child =  element->FirstChildElement("Interpolation");
                    textureMap.interpolation = "";
                    if( child != NULL )
                    {
                        std::string interpolation = child->GetText();
                        textureMap.interpolation = interpolation;
                    }

                    child =  element->FirstChildElement("BumpFactor");
                    if( child != NULL )
                    {
                        float bumpFactor = std::stoi(child->GetText());
                        textureMap.bumpFactor = bumpFactor;
                    }
                    //background = textureMap;
                    continue; 
                }
                    

                child =  element->FirstChildElement("Interpolation");
                textureMap.interpolation = "";
                if( child != NULL )
                {
                    std::string interpolation = child->GetText();
                    textureMap.interpolation = interpolation;
                }

                child =  element->FirstChildElement("BumpFactor");
                if( child != NULL )
                {
                    float bumpFactor = std::stoi(child->GetText());
                    textureMap.bumpFactor = bumpFactor;
                }
                
                

            } 
            if( element->Attribute("type" , "perlin") != NULL  )
            {
                textureMap.is_image = false;
                child =  element->FirstChildElement("NoiseConversion");
                std::string NoiseConversion = child->GetText();
                textureMap.noiseConversion = NoiseConversion;

                    

                child =  element->FirstChildElement("NoiseScale");
                float NoiseScale = std::stoi(child->GetText());
                textureMap.noiseScale = NoiseScale;

                child =  element->FirstChildElement("DecalMode");
                if( child != NULL )
                {
                    std::string decalMode = child->GetText();
                    textureMap.decalMode = decalMode;
                }
                

            } 
            
            textureMaps.push_back(textureMap);
            
            element = element->NextSiblingElement("TextureMap");
        }
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
            stream << child->GetText() << std::endl;
            stream >> material.ambient.x >> material.ambient.y >> material.ambient.z;
        }
        else
        {

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
            material.phong_exponent = 1.0f; 
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
        child = element->FirstChildElement("AbsorptionIndex");
        if( child != NULL)
        {
            stream << child->GetText() << std::endl;
            stream >> material.absorption_index;
        }
        else
        {
            material.absorption_index = 0.0f;
        }
        child = element->FirstChildElement("Roughness");
        if( child != NULL)
        {
            stream << child->GetText() << std::endl;
            stream >> material.roughness;
        }
        else
        {
            material.roughness = 0.0f;
        }
        materials.push_back(material);
        element = element->NextSiblingElement("Material");
    }
    //Get transformations
    element = root->FirstChildElement("Transformations");
    if ( element != NULL )
    {
   //Get Scaling
    element = element->FirstChildElement("Scaling");
    while( element)
    {
        

        stream << element->GetText() << std::endl;
        std::cout << "element get text " <<  element->GetText() << std::endl;
        Matrix s(4,4);
        Vec3f vertex;
        while (!(stream >> vertex.x).eof() )
        {
            stream >> vertex.y >> vertex.z;
        }
        std::cout << " vertex " << vertex.x << " " << vertex.y << " " << vertex.z    << std::endl;

        //scale matrix
        s.set(0, 0 , vertex.x );
        s.set(1, 1 , vertex.y );  
        s.set(2, 2 , vertex.z );  
        s.set(3, 3 , 1 );  
        transformations.scaling.push_back(s);

        element = element->NextSiblingElement("Scaling");
        stream.clear();

        if( element == 0 )
        {
            break;
        }
        
    }
    std::cout << " scaling done " << std::endl;
    stream.clear();
    //Get Translation
    element = root->FirstChildElement("Transformations");
    element = element->FirstChildElement("Translation");
    while( element)
    {
        stream << element->GetText() << std::endl;
        Matrix s(4,4);
        Vec3f vertex;
        while (!(stream >> vertex.x).eof() )
        {
            stream >> vertex.y >> vertex.z;
        }
        //translation matrix
        s.set(0, 0 , 1 );
        s.set(1, 1 , 1 );  
        s.set(2, 2 , 1 );  
        s.set(3, 3 , 1 );

        s.set(0, 3 , vertex.x );
        s.set(1, 3 , vertex.y );
        s.set(2, 3 , vertex.z );
        s.set(3, 3 , 1 );

        transformations.translation.push_back(s);
        element = element->NextSiblingElement("Translation");
        stream.clear();
        if( element == 0 )
        {
            break;
        }
        
    }
    std::cout << " translation  done " << std::endl;

    stream.clear();
    //Get Translation
    element = root->FirstChildElement("Transformations");
    element = element->FirstChildElement("Rotation");
    while( element)
    {
        stream << element->GetText() << std::endl;
        Matrix s(4,4);
        Vec4f vertex;
        while (!(stream >> vertex.x).eof() )
        {
            stream >> vertex.y >> vertex.z >> vertex.w;
        }

        float degree = vertex.x;
        float Rx = vertex.y;
        float Ry = vertex.z;
        float Rz = vertex.w;

        float cosine_theta = std::cos(M_PI/180  * degree );
        float sine_theta = std::sin(M_PI/180  * degree );

        //rotation matrix 
        //  x is degree 
        // y z w arbitrary axis
        s.set(0, 0 , cosine_theta + Rx * Rx * (1 - cosine_theta )  );
        s.set(0, 1 , Rx * Ry * (1 - cosine_theta ) - Rz * sine_theta  );  
        s.set(0, 2 , Rx * Rz * (1 -cosine_theta ) + Ry *sine_theta );  
        s.set(0, 3 , 0 );

        s.set(1, 0 , Ry*Rx*(1-cosine_theta) + Rz*sine_theta);
        s.set(1, 1 , cosine_theta + Ry*Ry*(1-cosine_theta) );
        s.set(1, 2 , Ry*Rz*(1-cosine_theta)-Rx*sine_theta );
        s.set(1, 3 , 0 );

        s.set(2, 0 , Ry*Rx*(1-cosine_theta) - Ry*sine_theta);
        s.set(2, 1 ,Rz*Ry*(1-cosine_theta) + Rx*sine_theta );
        s.set(2, 2 , cosine_theta + Rz*Rz*(1-cosine_theta) );
        s.set(2, 3 , 0 );

        s.set(3,0,0);
        s.set(3,1,0);
        s.set(3,2,0);
        s.set(3,3,1);


        transformations.rotation.push_back(s);
        element = element->NextSiblingElement("Rotation");
        stream.clear();
        if( element == 0 )
        {
            break;
        }
        
    }
    std::cout << " rotation  done " << std::endl;

    }
 
    stream.clear();
    //Get VertexData
    element = root->FirstChildElement("VertexData");
    stream << element->GetText() << std::endl;
    Vec3f vertex;
    while (!(stream >> vertex.x).eof() )
    {
         
        stream >> vertex.y >> vertex.z;
        vertex_data.push_back(vertex);

    }
    std::cout << " vertex data   done " << std::endl;
        std::cout << " here 0 " << std::endl;

    stream.clear();
    //Get VertexData
    element = root->FirstChildElement("TexCoordData");
    if( element != NULL )
    {
        stream << element->GetText() << std::endl;
        std::array<int, 2> tex;
        while (!(stream >> tex[0]).eof() )
        {
            stream >> tex[1];
            texcoords.push_back(tex);
        }
    }
    
    std::cout << " tex data   done " << std::endl;
    std::cout << " here 0 " << std::endl;
    //Get Meshes
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Mesh");
    while (element)
    {
        Mesh mesh;
        mesh.transformation.Identity(); // set I 
        mesh.transformation_inverse.Identity(); // set I 

        child = element->FirstChildElement("Material");
        mesh.material_id = std::stoi(child->GetText());

        child = element->FirstChildElement("Textures");
        stream << child->GetText() << std::endl;
        int texture_id; 
        while( stream >> texture_id  )
        {
            mesh.textures.push_back(texture_id);
        }
        stream.clear();


        child = element->FirstChildElement("Faces");

        

        if(!child->Attribute("plyFile")) // no ply
        {
            stream << child->GetText() << std::endl;
            
            Face face;
            while (!(stream >> face.v0_id).eof())
            {
                stream >> face.v1_id >> face.v2_id;
                mesh.faces.push_back(face);
            }
        }        
        else //ply file 
        {
            //chekc for ply
            const char* base_mesh = child->Attribute("plyFile");
            std::string base_mesh_str = base_mesh;
            happly::PLYData plyIn(base_mesh_str);
            std::vector<std::array<double , 3U > > vertex_pos = plyIn.getVertexPositions();
		    std::vector<std::vector<size_t>> face_indices=  plyIn.getFaceIndices();
            int vertex_data_initial_size = vertex_data.size();
            std::cout << " enter read " << std::endl; 
            for (size_t i = 0; i < vertex_pos.size(); i++)
            {
                //add them to vertex_data
                Vec3f point;
                point.x = (float)vertex_pos[i][0];
                point.y = (float)vertex_pos[i][1];
                point.z = (float)vertex_pos[i][2];

                vertex_data.push_back( point );
            }

            //now increase all of the face indicies by vertex_data_initial_size
            for (size_t i = 0; i < face_indices.size(); i++)
            {
                if(  face_indices[i].size() == 3 )
                {
                    Face face1;
                    face1.v0_id = face_indices[i][0] +  vertex_data_initial_size + 1;
                    face1.v1_id = face_indices[i][1] + vertex_data_initial_size  + 1 ;
                    face1.v2_id = face_indices[i][2] + vertex_data_initial_size + 1;
                    mesh.faces.push_back(face1 );
                }
                else if(face_indices[i].size() == 4  )
                {
                    Face face1;
                    face1.v0_id = face_indices[i][0] +  vertex_data_initial_size + 1;
                    face1.v1_id = face_indices[i][1] + vertex_data_initial_size  + 1 ;
                    face1.v2_id = face_indices[i][2] + vertex_data_initial_size + 1;
                    mesh.faces.push_back(face1 );
                    Face face2;
                    face2.v0_id = face_indices[i][0] +  vertex_data_initial_size + 1;
                    face2.v1_id = face_indices[i][2] + vertex_data_initial_size  + 1 ;
                    face2.v2_id = face_indices[i][3] + vertex_data_initial_size + 1;
                    mesh.faces.push_back(face2 );
                }
                


            }
            
            
        }
        child = element->FirstChildElement("Transformations");
        if( child != NULL )
        {
            stream << child->GetText() << std::endl;
            std::istringstream iss(child->GetText());
            std::string trans = child->GetText();
            while(std::getline(iss , trans , ' '))
            {
                char first = trans.at(0);
                std::cout << "substr ==>  " << trans << " " << trans.substr(0,trans.length()) << std::endl; 
                int no =  std::stoi(trans.substr(1, trans.length()) ) ;
                if( first == 't') //translate 
                {
                    Matrix m = transformations.translation[no-1];
                    mesh.transformations.push_back(m);

                    mesh.transformation = m * mesh.transformation; 
                }
                else if(first == 's') //sclae 
                {
                    Matrix m = transformations.scaling[no-1];
                    mesh.transformations.push_back(m);

                    mesh.transformation = m * mesh.transformation; 

                }
                else if( first == 'r') //rotate 
                {
                    Matrix m = transformations.rotation[no-1];
                    mesh.transformations.push_back(m);

                    mesh.transformation = m * mesh.transformation; 

                }
            }

            //get the inverse
            mesh.transformation_inverse = mesh.transformation.inverse();

            stream.clear();
        }
        child = element->FirstChildElement("MotionBlur");
        if( child != NULL )
        {
            std::istringstream iss(child->GetText());
            std::string trans = child->GetText();
            Vec3f motion_blur ;
            int count = 0; 
            while(std::getline(iss , trans , ' '))
            {
                if(count == 0 )
                {
                    motion_blur.x = std::stof(trans);
                }
                else if(count == 1 )
                {
                    motion_blur.y = std::stof(trans);
                }
                else if(count == 2 )
                {
                    motion_blur.z = std::stof(trans);
                }
                count++;

                mesh.motion_blur = motion_blur; 
            }

        }
        else
        {
            mesh.motion_blur.x = 0;
            mesh.motion_blur.y = 0;
            mesh.motion_blur.z = 0;
        }
        mesh.current_motion_blur.x = 0;
        mesh.current_motion_blur.y = 0;
        mesh.current_motion_blur.z = 0;


        meshes.push_back(mesh);
        mesh.faces.clear();
        element = element->NextSiblingElement("Mesh");
        stream.clear();
    }
    std::cout << "reading meshes done  " << std::endl;  
    stream.clear();


    //Get Mesh Instances


    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("MeshInstance");
    while (element)
    {
        MeshInstance mesh_instance;
        mesh_instance.transformation.Identity();
        mesh_instance.transformation_inverse.Identity();

        std::string base_mesh = element->Attribute("baseMeshId");
        mesh_instance.mesh_ptr = &this->meshes[std::stoi(base_mesh) - 1 ];
        std::cout << "base mesh id " << base_mesh << std::endl; 
        std::string reset_transform = element->Attribute("resetTransform");
        std::cout << "reset transfomr  " << reset_transform << std::endl; 
        
        if( reset_transform.at(0) == 'f')
        {
            //means false
            mesh_instance.reset_transform = false; 
            mesh_instance.transformation = this->meshes[std::stoi(base_mesh) - 1 ].transformation;

        }
        else
        {
            mesh_instance.reset_transform = true; 
            
        }

        child = element->FirstChildElement("Material");
        stream.clear();
        stream << child->GetText() << std::endl;
        std::cout << child->GetText() << std::endl; 
        mesh_instance.material_id = std::stoi(child->GetText());


        child = element->FirstChildElement("Transformations");
        if( child != NULL )
        {
            stream << child->GetText() << std::endl;
            std::istringstream iss(child->GetText());
            std::string trans = child->GetText();
            while(std::getline(iss , trans , ' '))
            {
                char first = trans.at(0);
                int no =  std::stoi(trans.substr(1, trans.length()) ) ;
                if( first == 't') //translate 
                {
                    Matrix m = transformations.translation[no-1];
                    mesh_instance.transformations.push_back(m);
                    mesh_instance.transformation = m *  mesh_instance.transformation;   
                }
                else if(first == 's') //sclae 
                {
                    Matrix m = transformations.scaling[no-1];
                    mesh_instance.transformations.push_back(m);
                    mesh_instance.transformation = m * mesh_instance.transformation ;   

                }
                else if( first == 'r') //rotate 
                {
                    Matrix m = transformations.rotation[no-1];
                    mesh_instance.transformations.push_back(m);
                    mesh_instance.transformation = m * mesh_instance.transformation;   

                }
            }
            mesh_instance.transformation_inverse =     mesh_instance.transformation.inverse();
            stream.clear();
        }
        child = element->FirstChildElement("MotionBlur");
        if( child != NULL )
        {
            std::istringstream iss(child->GetText());
            std::string trans = child->GetText();
            Vec3f motion_blur ;
            int count = 0; 
            while(std::getline(iss , trans , ' '))
            {
                if(count == 0 )
                {
                    motion_blur.x = std::stof(trans);
                }
                else if(count == 1 )
                {
                    motion_blur.y = std::stof(trans);
                }
                else if(count == 2 )
                {
                    motion_blur.z = std::stof(trans);
                }
                count++;

                mesh_instance.motion_blur = motion_blur; 
            }

        }
        else
        {
            mesh_instance.motion_blur.x = 0;
            mesh_instance.motion_blur.y = 0;
            mesh_instance.motion_blur.z = 0;
        }
        mesh_instance.current_motion_blur.x = 0;
        mesh_instance.current_motion_blur.y = 0;
        mesh_instance.current_motion_blur.z = 0;
        mesh_instances.push_back(mesh_instance);
        element = element->NextSiblingElement("MeshInstance");
    } 
    std::cout << " mesh instances " << std::endl;

    //Get Triangles
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Triangle");
    while (element)
    {
        Triangle triangle;
        triangle.transformation.Identity();
        triangle.transformation_inverse.Identity();

        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> triangle.material_id;

        child = element->FirstChildElement("Indices");
        stream << child->GetText() << std::endl;
        stream >> triangle.indices.v0_id >> triangle.indices.v1_id >> triangle.indices.v2_id;
        child = element->FirstChildElement("Transformations");
        if( child != NULL )
        {
            stream << child->GetText() << std::endl;
            std::istringstream iss(child->GetText());
            std::string trans = child->GetText();
            while(std::getline(iss , trans , ' '))
            {
                char first = trans.at(0);
                int no =  std::stoi(trans.substr(1, trans.length()) ) ;
                if( first == 't') //translate 
                {
                    Matrix m = transformations.translation[no-1];
                    triangle.transformations.push_back(m);

                    triangle.transformation = m * triangle.transformation;
                }
                else if(first == 's') //sclae 
                {
                    Matrix m = transformations.scaling[no-1];
                    triangle.transformations.push_back(m);
                    triangle.transformation = m * triangle.transformation;

                }
                else if( first == 'r') //rotate 
                {
                    Matrix m = transformations.rotation[no-1];
                    triangle.transformations.push_back(m);
                    triangle.transformation = m * triangle.transformation;

                }

                triangle.transformation_inverse = triangle.transformation.inverse();
            }
        }
        stream.clear();
        triangles.push_back(triangle);
        element = element->NextSiblingElement("Triangle");
    }
    std::cout << "reading triangles done  " << std::endl;  

    //Get Spheres
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Sphere");
    while (element)
    {
        Sphere sphere;
        sphere.transformation.Identity();
        sphere.transformation_inverse.Identity();

        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        sphere.material_id = std::stoi(child->GetText());
        std::cout << sphere.material_id << std::endl; 

        child = element->FirstChildElement("Center");
        stream << child->GetText() << std::endl;
        sphere.center_vertex_id = std::stoi(child->GetText());

        child = element->FirstChildElement("Radius");
        stream << child->GetText() << std::endl;
        stream >> sphere.radius;
        sphere.radius = std::stof(child->GetText());

        stream.clear();

        child = element->FirstChildElement("Textures");
        if( child != NULL )
        {
            std::istringstream iss(child->GetText());
            std::string temp;
            while( std::getline(iss ,temp , ' ' ))
            {
                int texture_id; 
                sphere.textures.push_back(std::stoi(temp));
            }
        }

        
        stream.clear();
        


        child = element->FirstChildElement("Transformations");
        if( child != NULL )
        {
            stream << child->GetText() << std::endl;
            std::istringstream iss(child->GetText());
            std::string trans = child->GetText();
            while(std::getline(iss , trans , ' '))
            {
                char first = trans.at(0);
                int no =  std::stoi(trans.substr(1, trans.length()) ) ;
                if( first == 't') //translate 
                {
                    Matrix m = transformations.translation[no-1];
                    sphere.transformations.push_back(m);
                    sphere.transformation = m * sphere.transformation;
                }
                else if(first == 's') //scale 
                {
                    Matrix m = transformations.scaling[no-1];
                    sphere.transformations.push_back(m);
                    sphere.scale_indices.push_back( sphere.transformations.size() -1) ;
                    sphere.transformation = m * sphere.transformation ;

                }
                else if( first == 'r') //rotate 
                {
                    Matrix m = transformations.rotation[no-1];
                    sphere.transformations.push_back(m);
                    sphere.transformation = m * sphere.transformation ;
                }
            }
            sphere.transformation_inverse = sphere.transformation.inverse();
        }
        stream.clear();
        child = element->FirstChildElement("MotionBlur");
        if( child != NULL )
        {
            std::istringstream iss(child->GetText());
            std::string trans = child->GetText();
            Vec3f motion_blur ;
            int count = 0; 
            while(std::getline(iss , trans , ' '))
            {
                if(count == 0 )
                {
                    motion_blur.x = std::stof(trans);
                }
                else if(count == 1 )
                {
                    motion_blur.y = std::stof(trans);
                }
                else if(count == 2 )
                {
                    motion_blur.z = std::stof(trans);
                }
                count++;

                sphere.motion_blur = motion_blur; 
            }

        }
        else
        {
            sphere.motion_blur.x = 0;
            sphere.motion_blur.y = 0;
            sphere.motion_blur.z = 0;
        }
        sphere.current_motion_blur.x = 0;
        sphere.current_motion_blur.y = 0;
        sphere.current_motion_blur.z = 0;
        spheres.push_back(sphere);
        element = element->NextSiblingElement("Sphere");
    }
    std::cout << "reading spheres done  " << std::endl;  

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
    cross.y = -1 * ( (vec1.x * vec2.z) - (vec1.z * vec2.x) ) ;
    cross.z = (vec1.x * vec2.y) - (vec1.y * vec2.x);
    return parser::normalize(cross);
}
parser::Vec3f parser::cross_non_normalize( const parser::Vec3f& vec1 , const parser::Vec3f & vec2 )
{
    parser::Vec3f cross; 
    cross.x = (vec1.y * vec2.z) - (vec1.z * vec2.y);
    cross.y = -1 * ( (vec1.x * vec2.z) - (vec1.z * vec2.x) ) ;
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
    if( std::abs( len - 0.0f ) <  1e-4  )
    {
        new_vec = vec1;
        return new_vec; 
    }
    new_vec.x =  vec1.x / len;
    new_vec.y =  vec1.y / len;
    new_vec.z =  vec1.z / len;
    return new_vec;

}
float parser::distance(const  parser::Vec3f &  vec1 ,  const  parser::Vec3f & vec2   )
{
    return std::sqrt(  (vec1.x - vec2.x ) * (vec1.x - vec2.x )  + (vec1.y - vec2.y ) * (vec1.y - vec2.y ) + (vec1.z - vec2.z ) * (vec1.z - vec2.z )   ); 
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
parser::Vec2f::Vec2f(float x  , float y )
{
    this->x = x;
    this->y = y;

}
parser::Vec2f::Vec2f()
{
    this->x = 0;
    this->y = 0;

}
parser::Vec2f parser::normalize_vec2f( const  parser::Vec2f &  vec1    )
{
    parser::Vec2f new_vec;
    float len =  std::sqrt( vec1.x*vec1.x + vec1.y*vec1.y );
    new_vec.x = vec1.x / len;
    new_vec.y = vec1.y / len; 

    return new_vec;
}

float parser::dot_vec2f( const parser::Vec2f & v1 , const parser::Vec2f & v2  )
{
    return v1.x * v2.x + v1.y * v2.y;
} 
parser::Vec3f parser::Vec3f::operator +(const parser::Vec3f & vec1   )
{
    parser::Vec3f vec;
    vec.x =  this->x +  vec1.x ; 
    vec.y =  this->y +  vec1.y ; 
    vec.z =  this->z +  vec1.z ; 
    return vec;
}
parser::Vec3f parser::Vec3f::operator -(const parser::Vec3f & vec1   )
{
    parser::Vec3f vec;
    vec.x =  this->x - vec1.x ; 
    vec.y =  this->y - vec1.y ; 
    vec.z =  this->z - vec1.z ;
    return vec;

}
parser::Vec3f parser::Vec3f::operator *(const parser::Vec3f & vec1   )
{
    parser::Vec3f vec;
    vec.x = this->x * vec1.x ; 
    vec.y = this->y * vec1.y ; 
    vec.z = this->z * vec1.z ; 
    return vec;

}
parser::Vec3f parser::Vec3f::operator /(const parser::Vec3f & vec1   )
{
    parser::Vec3f vec;
    vec.x = this->x / vec1.x  ; 
    vec.y = this->y / vec1.y  ; 
    vec.z = this->z / vec1.z  ; 
    return vec;

}
// int 
parser::Vec3f parser::Vec3f::operator *(const int & num   )
{
    parser::Vec3f vec;
    vec.x = this->x * num ; 
    vec.y = this->y * num ; 
    vec.z = this->z * num ; 
    return vec;

}
parser::Vec3f parser::Vec3f::operator /(const int & num   )
{
    parser::Vec3f vec;
    vec.x = this->x / num; 
    vec.y = this->y / num; 
    vec.z = this->z / num; 
    return vec;

}
parser::Vec3f parser::Vec3f::operator +(const int & num   )
{
    parser::Vec3f vec;
    vec.x =  this->x + num; 
    vec.y =  this->y + num; 
    vec.z =  this->z + num; 
    return vec;

}
parser::Vec3f parser::Vec3f::operator -(const int & num   )
{
    parser::Vec3f vec;
    vec.x =  this->x - num ; 
    vec.y =  this->y - num ; 
    vec.z =  this->z - num ;
    return vec;
}

parser::Vec3f parser::Vec3f::operator *(const float & num   )
{
    parser::Vec3f vec;
    vec.x =  this->x * num ; 
    vec.y =  this->y * num ; 
    vec.z =  this->z * num ; 
    return vec;

}
parser::Vec3f parser::Vec3f::operator /(const float & num   )
{
    parser::Vec3f vec;
    vec.x = this->x / num; 
    vec.y = this->y / num; 
    vec.z = this->z / num; 
    return vec;

}
parser::Vec3f parser::Vec3f::operator +(const float & num   )
{
    parser::Vec3f vec;
    vec.x = this->x  + num; 
    vec.y = this->y  + num; 
    vec.z = this->z  + num; 
    return vec;

}
parser::Vec3f parser::Vec3f::operator -(const float & num   )
{
    parser::Vec3f vec;
    vec.x =  this->x - num; 
    vec.y =  this->y - num; 
    vec.z =  this->z - num; 
    return vec;
}
void parser::Vec3f::operator =(const parser::Vec3f & vec1 )
{
    this->x = vec1.x;
    this->y = vec1.y;
    this->z = vec1.z;
}

//matrix
parser::Matrix::Matrix(int row_no , int col_no )
{
    this->row_no = row_no;
    this->col_no = col_no;

    for (size_t i = 0; i < row_no * col_no; i++)
    {
        elements.push_back(0.0f);
    }
    
}
parser::Matrix::Matrix()
{
    this->row_no = 4;
    this->col_no = 4;

    for (size_t i = 0; i < row_no * col_no; i++)
    {
        elements.push_back(0.0f);
    }
    
}
float parser::Matrix::at(int row_index , int col_index)
{
    return elements[ row_index * this->col_no  + col_index];
}
void parser::Matrix::Identity()
{
    parser::Matrix mat(this->row_no , this->col_no);
    this->elements = mat.elements; 
    for (size_t i = 0; i < row_no; i++)
    {
        this->elements[i * col_no + i ] = 1.0f; 
    }
}
void parser::Matrix::set(int row , int col , float num )
{
    elements[ row * this->col_no  + col] = num; 
}

void parser::print_matrix(const parser::Matrix & mat )
{
    for (size_t i = 0; i < mat.row_no; i++)
    {
        for (size_t j = 0; j < mat.col_no; j++)
        {
           std::cout << mat.elements[i * mat.col_no + j] << " ";
        }
        std::cout<<"\n";
        
    }
    
}
parser::Vec4f parser::Matrix::operator *(const parser::Vec4f & vec1   ) // 4x4 matrix vs 4x1 vector 
{
    parser::Vec4f new_vec;
    new_vec.x = this->at(0,0) * vec1.x + this->at(0,1) * vec1.y + this->at(0,2) * vec1.z +this->at(0,3) * vec1.w; 
    new_vec.y = this->at(1,0) * vec1.x + this->at(1,1) * vec1.y + this->at(1,2) * vec1.z +this->at(1,3) * vec1.w; 
    new_vec.z = this->at(2,0) * vec1.x + this->at(2,1) * vec1.y + this->at(2,2) * vec1.z +this->at(2,3) * vec1.w;
    new_vec.w = this->at(3,0) * vec1.x + this->at(3,1) * vec1.y + this->at(3,2) * vec1.z +this->at(3,3) * vec1.w ;

    return new_vec;
}

parser::Matrix parser::Matrix::operator*( parser::Matrix & matrix   ) // 4x4 matrix vs 4x4 vector 
{
    parser::Matrix mat(4,4 );

    mat.set(0 , 0 , this->at( 0,0 ) * matrix.at(0,0) + this->at( 0,1 ) * matrix.at(1,0) + this->at( 0,2 ) * matrix.at(2,0) + this->at( 0,3 ) * matrix.at(3,0)  );
    mat.set(0 , 1 , this->at( 0,0 ) * matrix.at(0,1) + this->at( 0,1 ) * matrix.at(1,1) + this->at( 0,2 ) * matrix.at(2,1) + this->at( 0,3 ) * matrix.at(3,1)  );
    mat.set(0 , 2 , this->at( 0,0 ) * matrix.at(0,2) + this->at( 0,1 ) * matrix.at(1,2) + this->at( 0,2 ) * matrix.at(2,2) + this->at( 0,3 ) * matrix.at(3,2)  );
    mat.set(0 , 3 , this->at( 0,0 ) * matrix.at(0,3) + this->at( 0,1 ) * matrix.at(1,3) + this->at( 0,2 ) * matrix.at(2,3) + this->at( 0,3 ) * matrix.at(3,3)  );

    mat.set(1 , 0 , this->at( 1,0 ) * matrix.at(0,0) + this->at( 1,1 ) * matrix.at(1,0) + this->at( 1,2 ) * matrix.at(2,0) + this->at( 1,3 ) * matrix.at(3,0)  );
    mat.set(1 , 1 , this->at( 1,0 ) * matrix.at(0,1) + this->at( 1,1 ) * matrix.at(1,1) + this->at( 1,2 ) * matrix.at(2,1) + this->at( 1,3 ) * matrix.at(3,1)  );
    mat.set(1 , 2 , this->at( 1,0 ) * matrix.at(0,2) + this->at( 1,1 ) * matrix.at(1,2) + this->at( 1,2 ) * matrix.at(2,2) + this->at( 1,3 ) * matrix.at(3,2)  );
    mat.set(1 , 3 , this->at( 1,0 ) * matrix.at(0,3) + this->at( 1,1 ) * matrix.at(1,3) + this->at( 1,2 ) * matrix.at(2,3) + this->at( 1,3 ) * matrix.at(3,3)  );

    mat.set(2 , 0 , this->at( 2,0 ) * matrix.at(0,0) + this->at( 2,1 ) * matrix.at(1,0) + this->at( 2,2 ) * matrix.at(2,0) + this->at( 2,3 ) * matrix.at(3,0)  );
    mat.set(2 , 1 , this->at( 2,0 ) * matrix.at(0,1) + this->at( 2,1 ) * matrix.at(1,1) + this->at( 2,2 ) * matrix.at(2,1) + this->at( 2,3 ) * matrix.at(3,1)  );
    mat.set(2 , 2 , this->at( 2,0 ) * matrix.at(0,2) + this->at( 2,1 ) * matrix.at(1,2) + this->at( 2,2 ) * matrix.at(2,2) + this->at( 2,3 ) * matrix.at(3,2)  );
    mat.set(2 , 3 , this->at( 2,0 ) * matrix.at(0,3) + this->at( 2,1 ) * matrix.at(1,3) + this->at( 2,2 ) * matrix.at(2,3) + this->at( 2,3 ) * matrix.at(3,3)  );

    mat.set(3 , 0 , this->at( 3,0 ) * matrix.at(0,0) + this->at( 3,1 ) * matrix.at(1,0) + this->at( 3,2 ) * matrix.at(2,0) + this->at( 3,3 ) * matrix.at(3,0)  );
    mat.set(3 , 1 , this->at( 3,0 ) * matrix.at(0,1) + this->at( 3,1 ) * matrix.at(1,1) + this->at( 3,2 ) * matrix.at(2,1) + this->at( 3,3 ) * matrix.at(3,1)  );
    mat.set(3 , 2 , this->at( 3,0 ) * matrix.at(0,2) + this->at( 3,1 ) * matrix.at(1,2) + this->at( 3,2 ) * matrix.at(2,2) + this->at( 3,3 ) * matrix.at(3,2)  );
    mat.set(3 , 3 , this->at( 3,0 ) * matrix.at(0,3) + this->at( 3,1 ) * matrix.at(1,3) + this->at( 3,2 ) * matrix.at(2,3) + this->at( 3,3 ) * matrix.at(3,3)  );

    return mat; 
    
}
parser::Matrix parser::Matrix::inverse()
{
    Matrix inv(4,4);
    double det;
    int i;
    
    inv.elements[0] = this->elements[5]  * this->elements[10] * this->elements[15] - 
             this->elements[5]  * this->elements[11] * this->elements[14] - 
             this->elements[9]  * this->elements[6]  * this->elements[15] + 
             this->elements[9]  * this->elements[7]  * this->elements[14] +
             this->elements[13] * this->elements[6]  * this->elements[11] - 
             this->elements[13] * this->elements[7]  * this->elements[10];

    inv.elements[4] = -this->elements[4]  * this->elements[10] * this->elements[15] + 
              this->elements[4]  * this->elements[11] * this->elements[14] + 
              this->elements[8]  * this->elements[6]  * this->elements[15] - 
              this->elements[8]  * this->elements[7]  * this->elements[14] - 
              this->elements[12] * this->elements[6]  * this->elements[11] + 
              this->elements[12] * this->elements[7]  * this->elements[10];

    inv.elements[8] = this->elements[4]  * this->elements[9] * this->elements[15] - 
             this->elements[4]  * this->elements[11] * this->elements[13] - 
             this->elements[8]  * this->elements[5] * this->elements[15] + 
             this->elements[8]  * this->elements[7] * this->elements[13] + 
             this->elements[12] * this->elements[5] * this->elements[11] - 
             this->elements[12] * this->elements[7] * this->elements[9];

    inv.elements[12] = -this->elements[4]  * this->elements[9] * this->elements[14] + 
               this->elements[4]  * this->elements[10] * this->elements[13] +
               this->elements[8]  * this->elements[5] * this->elements[14] - 
               this->elements[8]  * this->elements[6] * this->elements[13] - 
               this->elements[12] * this->elements[5] * this->elements[10] + 
               this->elements[12] * this->elements[6] * this->elements[9];

    inv.elements[1] = -this->elements[1]  * this->elements[10] * this->elements[15] + 
              this->elements[1]  * this->elements[11] * this->elements[14] + 
              this->elements[9]  * this->elements[2] * this->elements[15] - 
              this->elements[9]  * this->elements[3] * this->elements[14] - 
              this->elements[13] * this->elements[2] * this->elements[11] + 
              this->elements[13] * this->elements[3] * this->elements[10];

    inv.elements[5] = this->elements[0]  * this->elements[10] * this->elements[15] - 
             this->elements[0]  * this->elements[11] * this->elements[14] - 
             this->elements[8]  * this->elements[2] * this->elements[15] + 
             this->elements[8]  * this->elements[3] * this->elements[14] + 
             this->elements[12] * this->elements[2] * this->elements[11] - 
             this->elements[12] * this->elements[3] * this->elements[10];

    inv.elements[9] = -this->elements[0]  * this->elements[9] * this->elements[15] + 
              this->elements[0]  * this->elements[11] * this->elements[13] + 
              this->elements[8]  * this->elements[1] * this->elements[15] - 
              this->elements[8]  * this->elements[3] * this->elements[13] - 
              this->elements[12] * this->elements[1] * this->elements[11] + 
              this->elements[12] * this->elements[3] * this->elements[9];

    inv.elements[13] = this->elements[0]  * this->elements[9] * this->elements[14] - 
              this->elements[0]  * this->elements[10] * this->elements[13] - 
              this->elements[8]  * this->elements[1] * this->elements[14] + 
              this->elements[8]  * this->elements[2] * this->elements[13] + 
              this->elements[12] * this->elements[1] * this->elements[10] - 
              this->elements[12] * this->elements[2] * this->elements[9];

    inv.elements[2] = this->elements[1]  * this->elements[6] * this->elements[15] - 
             this->elements[1]  * this->elements[7] * this->elements[14] - 
             this->elements[5]  * this->elements[2] * this->elements[15] + 
             this->elements[5]  * this->elements[3] * this->elements[14] + 
             this->elements[13] * this->elements[2] * this->elements[7] - 
             this->elements[13] * this->elements[3] * this->elements[6];

    inv.elements[6] = -this->elements[0]  * this->elements[6] * this->elements[15] + 
              this->elements[0]  * this->elements[7] * this->elements[14] + 
              this->elements[4]  * this->elements[2] * this->elements[15] - 
              this->elements[4]  * this->elements[3] * this->elements[14] - 
              this->elements[12] * this->elements[2] * this->elements[7] + 
              this->elements[12] * this->elements[3] * this->elements[6];

    inv.elements[10] = this->elements[0]  * this->elements[5] * this->elements[15] - 
              this->elements[0]  * this->elements[7] * this->elements[13] - 
              this->elements[4]  * this->elements[1] * this->elements[15] + 
              this->elements[4]  * this->elements[3] * this->elements[13] + 
              this->elements[12] * this->elements[1] * this->elements[7] - 
              this->elements[12] * this->elements[3] * this->elements[5];

    inv.elements[14] = -this->elements[0]  * this->elements[5] * this->elements[14] + 
               this->elements[0]  * this->elements[6] * this->elements[13] + 
               this->elements[4]  * this->elements[1] * this->elements[14] - 
               this->elements[4]  * this->elements[2] * this->elements[13] - 
               this->elements[12] * this->elements[1] * this->elements[6] + 
               this->elements[12] * this->elements[2] * this->elements[5];

    inv.elements[3] = -this->elements[1] * this->elements[6] * this->elements[11] + 
              this->elements[1] * this->elements[7] * this->elements[10] + 
              this->elements[5] * this->elements[2] * this->elements[11] - 
              this->elements[5] * this->elements[3] * this->elements[10] - 
              this->elements[9] * this->elements[2] * this->elements[7] + 
              this->elements[9] * this->elements[3] * this->elements[6];

    inv.elements[7] = this->elements[0] * this->elements[6] * this->elements[11] - 
             this->elements[0] * this->elements[7] * this->elements[10] - 
             this->elements[4] * this->elements[2] * this->elements[11] + 
             this->elements[4] * this->elements[3] * this->elements[10] + 
             this->elements[8] * this->elements[2] * this->elements[7] - 
             this->elements[8] * this->elements[3] * this->elements[6];

    inv.elements[11] = -this->elements[0] * this->elements[5] * this->elements[11] + 
               this->elements[0] * this->elements[7] * this->elements[9] + 
               this->elements[4] * this->elements[1] * this->elements[11] - 
               this->elements[4] * this->elements[3] * this->elements[9] - 
               this->elements[8] * this->elements[1] * this->elements[7] + 
               this->elements[8] * this->elements[3] * this->elements[5];

    inv.elements[15] = this->elements[0] * this->elements[5] * this->elements[10] - 
              this->elements[0] * this->elements[6] * this->elements[9] - 
              this->elements[4] * this->elements[1] * this->elements[10] + 
              this->elements[4] * this->elements[2] * this->elements[9] + 
              this->elements[8] * this->elements[1] * this->elements[6] - 
              this->elements[8] * this->elements[2] * this->elements[5];

    det = this->elements[0] * inv.elements[0] + this->elements[1] * inv.elements[4] + this->elements[2] * inv.elements[8] + this->elements[3] * inv.elements[12];

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        inv.elements[i] = inv.elements[i] * det;

    return inv;
}
void parser::Matrix::Transpose()
{
    std::vector<float> temp_element  =this->elements; 
    this->elements[0] = temp_element[0];
    this->elements[1] = temp_element[4];
    this->elements[2] = temp_element[8];
    this->elements[3] = temp_element[12]; 

    this->elements[4] = temp_element[1]; 
    this->elements[5] = temp_element[5];
    this->elements[6] = temp_element[9];
    this->elements[7] = temp_element[13]; 

    this->elements[8] = temp_element[2]; 
    this->elements[9] = temp_element[6];
    this->elements[10] = temp_element[10];
    this->elements[11] = temp_element[14]; 

    this->elements[12] = temp_element[3]; 
    this->elements[13] = temp_element[7];
    this->elements[14] = temp_element[11];
    this->elements[15] = temp_element[15]; 
}
