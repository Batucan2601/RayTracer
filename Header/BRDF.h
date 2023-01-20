#include "../System_Files/support_files/parser.h"

parser::Vec3f get_BRDF(const parser::Scene & scene , const parser::Material &material , const float& alpha_angle , const float& theta_angle , const float & cos_half_vec , const parser::Camera & current_camera  , const parser::Vec3f & wi , const parser::Vec3f & wo , const parser::Vec3f & wh ,const parser::Vec3f & normal     )
{
    
    parser::BRDF brdf = material.brdf;
    if( brdf.type == "OriginalPhong") // 2 
    {
        if( theta_angle < 0 )
        {
            return parser::Vec3f(0.f , 0.0f , 0.0f ); 
        }
        parser::Vec3f temp_specular;

        temp_specular.x = material.specular.x  *  std::pow( alpha_angle  ,  brdf.exponent) / theta_angle ; // wtff ? 
        temp_specular.y = material.specular.y  *  std::pow( alpha_angle  ,  brdf.exponent) / theta_angle ;
        temp_specular.z = material.specular.z  *  std::pow( alpha_angle  ,  brdf.exponent) / theta_angle ;
        return  ( (parser::Vec3f)material.diffuse +  temp_specular );
    }
    else if(  brdf.type == "ModifiedPhong")  
    {
        parser::Vec3f temp_specular;
        if( theta_angle < 0 )
        {
            return parser::Vec3f(0.f , 0.0f , 0.0f ); 
        }
        if( brdf.is_normalized) // 3
        {
            temp_specular.x = material.specular.x  *  std::pow( alpha_angle  ,  brdf.exponent ) * ( brdf.exponent + 2 )  ;
            temp_specular.y = material.specular.y  *  std::pow( alpha_angle  ,  brdf.exponent ) * ( brdf.exponent + 2 )  ;
            temp_specular.z = material.specular.z  *  std::pow( alpha_angle  ,  brdf.exponent ) * ( brdf.exponent + 2 )  ;
            return  ( ((parser::Vec3f)material.diffuse / (float)M_PI)  +  (temp_specular  / (2.0f* (float)M_PI)) );
        }
        else // 4 
        {
            temp_specular.x = material.specular.x  *  std::pow( alpha_angle  ,  brdf.exponent)  ;
            temp_specular.y = material.specular.y  *  std::pow( alpha_angle  ,  brdf.exponent)  ;
            temp_specular.z = material.specular.z  *  std::pow( alpha_angle  ,  brdf.exponent)  ;
            return  ( (parser::Vec3f)material.diffuse +  temp_specular );

        }
        

    }
    else if( brdf.type == "OriginalBlinnPhong")  // 5 
    {
        if( cos_half_vec < 0  )
        {
            return parser::Vec3f(0.f , 0.0f , 0.0f ); 
        }

        parser::Vec3f temp_specular;
        temp_specular.x = material.specular.x  *  std::pow(  alpha_angle  ,  brdf.exponent) ;
        temp_specular.y = material.specular.y  *  std::pow(  alpha_angle  ,  brdf.exponent) ;
        temp_specular.z = material.specular.z  *  std::pow(  alpha_angle  ,  brdf.exponent) ;

        return  ( (parser::Vec3f)material.diffuse +  (temp_specular/cos_half_vec) );
    }
    else if( brdf.type == "ModifiedBlinnPhong")
    {
        if( theta_angle < 0 )
        {
            return parser::Vec3f(0.f , 0.0f , 0.0f ); 
        }
        
        parser::Vec3f temp_specular;
        if( brdf.is_normalized ) // 7 
        {
            temp_specular.x = material.specular.x  * (brdf.exponent + 8) *  std::pow(  alpha_angle  , brdf.exponent) ;
            temp_specular.y = material.specular.y  * (brdf.exponent + 8) *  std::pow(  alpha_angle  , brdf.exponent) ;
            temp_specular.z = material.specular.z  * (brdf.exponent + 8) *  std::pow(  alpha_angle  , brdf.exponent) ;
            

            return  (  ( (parser::Vec3f)material.diffuse /  (float)M_PI )   +   ( temp_specular /   (float) (8 * M_PI)  )   );
        }
        else // 6 
        {
            temp_specular.x = material.specular.x  *  std::pow(  alpha_angle  ,  brdf.exponent) ;
            temp_specular.y = material.specular.y  *  std::pow(  alpha_angle  ,  brdf.exponent) ;
            temp_specular.z = material.specular.z  *  std::pow(  alpha_angle  ,  brdf.exponent) ;
            

            return  ( (parser::Vec3f)material.diffuse +  temp_specular );
        }
        
    }
    else // torrance sparrow 8 
    {
        if( theta_angle < 0 )
        {
            return parser::Vec3f(0.f , 0.0f , 0.0f ); 
        }
        
        parser::Vec3f temp_specular;
        // step 1 - we already have wh 

        // step 2 
        float alpha = parser::dot( wh , normal) / ( parser::length(wh) * parser::length(normal));

        //step 3
        float d_alpha = (brdf.exponent + 2) / (2 * M_PI)  * std::pow( alpha, brdf.exponent);

        //step 4

        float g = std::min( 1.0f , std::min( 2 * (parser::dot(normal , wh ) / (parser::length(normal) * parser::length(wh) ) ) * (parser::dot(normal , wo ) / (parser::length(normal) * parser::length(wo ) ) ) /  (parser::dot(wo , wh ) / (parser::length(wo) * parser::length(wh) ) ) , 2 * (parser::dot(normal , wh ) / (parser::length(normal) * parser::length(wh) ) ) * (parser::dot(normal , wi ) / (parser::length(normal) * parser::length(wi ) ) ) /  (parser::dot(wo , wh ) / (parser::length(wi) * parser::length(wh) ) )  ) ); 

        //step 5 
        float ro =  (material.refraction_index - 1) * (material.refraction_index - 1) / ( (material.refraction_index + 1) * (material.refraction_index + 1) ) ;
        float f =  ro + ( 1 - ro ) * std::pow( 1 -  ( parser::dot(wo , wh )  / (parser::length(wo) * parser::length(wh) )  )  , 5  );


        temp_specular = (parser::Vec3f )material.specular * ( d_alpha * f * g )  /   ( 4 * theta_angle  ) ;
    
        return ( (parser::Vec3f)material.diffuse /  (float)M_PI )  + temp_specular;
    }
    
    
    
    
}