#include "../System_Files/support_files/parser.h"

parser::Vec3f get_BRDF(const parser::Scene & scene , const parser::Material &material , const float& alpha_angle , const float& theta_angle    )
{
    if( theta_angle < 0 )
    {
        return parser::Vec3f(0.f , 0.0f , 0.0f ); 
    }
    parser::BRDF brdf = material.brdf;
    if( brdf.type == "OriginalPhong")
    {
        parser::Vec3f temp_specular;
        temp_specular.x = material.specular.x  *  std::pow( alpha_angle  ,  brdf.exponent ) /   theta_angle  ;
        temp_specular.y = material.specular.y  *  std::pow( alpha_angle  ,  brdf.exponent ) /   theta_angle  ;
        temp_specular.z = material.specular.z  *  std::pow( alpha_angle  ,  brdf.exponent ) /   theta_angle  ;
        

        return  ( (parser::Vec3f)material.diffuse +  temp_specular );
    }
    if( brdf.type == "ModifiedBlinnPhong")
    {
         parser::Vec3f temp_specular;
        
        if( brdf.is_normalized )
        {
            temp_specular.x = material.specular.x  * (brdf.exponent + 2) *  std::pow(  alpha_angle  ,  scene.cameras[0].number_of_samples) ;
            temp_specular.y = material.specular.y  * (brdf.exponent + 2) *  std::pow(  alpha_angle  , scene.cameras[0].number_of_samples) ;
            temp_specular.z = material.specular.z  * (brdf.exponent + 2) *  std::pow(  alpha_angle  ,  scene.cameras[0].number_of_samples) ;
            

            return  ( (parser::Vec3f)material.diffuse /  (float)M_PI  +  temp_specular /   (float) (2 * M_PI)    );
        }
        else
        {
            temp_specular.x = material.specular.x  *  std::pow(  alpha_angle  ,  brdf.exponent) ;
            temp_specular.y = material.specular.y  *  std::pow(  alpha_angle  ,  brdf.exponent) ;
            temp_specular.z = material.specular.z  *  std::pow(  alpha_angle  ,  brdf.exponent) ;
            

            return  ( (parser::Vec3f)material.diffuse +  temp_specular );
        }
        
    }
    
}