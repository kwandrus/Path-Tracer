
#include "AreaLight.h"

AreaLight::AreaLight(const Point& position, const Color& color, const Vector& a, const Vector& b)
  : cornerPos(position), color(color), a(a), b(b)
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
    double rand1 = double(rand()) / double(RAND_MAX);
    double rand2 = double(rand()) / double(RAND_MAX);
    Point samplePt = cornerPos + rand1 * a + rand2 * b;

    light_color = color;
    Vector dir = samplePt - hitpos;
    double len = dir.normalize();
    light_direction = dir;
    return len;
}

// used in direct rendering of area light (not its direct illumination on other surfaces)
Color AreaLight::areaShade(const Ray& ray, Vector normal) const{
    // source: Ray Tracing from the Ground Up, pg 332

    if (Dot(-normal, ray.direction()) > 0.0)
        return color * ls;
    else
        return Color(0, 0, 0);
}
