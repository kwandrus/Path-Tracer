
#include "AreaLight.h"

AreaLight::AreaLight(const Point& position, const Color& color)
  : center(position), color(color)
{
}

AreaLight::~AreaLight()
{
}

void AreaLight::preprocess()
{
}

double AreaLight::getLight(Color& light_color, Vector& light_direction,
                            const RenderContext&, const Point& hitpos) const
{
    // generate sample point on the light
    Point sample_pt = 

    light_color = color;
    Vector dir = center - hitpos;
    double len = dir.normalize();
    light_direction = dir;
    return len;
}
