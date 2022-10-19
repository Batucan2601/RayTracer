#shader vertex
#version 330 core           
layout( location = 0 ) in vec3 position;
void main()   
 {
  
  gl_Position = position;   
 };



#shader fragment
#version 330 core      
out vec4 FragColor;
void main()   
  { 
    
    FragColor = vec4( 1.0f,1.0f,1.0f , 1.0f);
  };