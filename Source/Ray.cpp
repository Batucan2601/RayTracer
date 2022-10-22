#include "../Header/Ray.h"
#include "../System_Files/support_files/parser.h"
#include "../Header/Prototypes.h"
Ray::Ray(glm::vec3 origin , glm::vec3 point2 )
{
    this->origin = origin;
    this->direction = glm::normalize(origin - point2);
}



