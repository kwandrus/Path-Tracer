
#include "LambertianMaterial.h"
#include "HitRecord.h"
#include "Light.h"
#include "Point.h"
#include "Primitive.h"
#include "Ray.h"
#include "RenderContext.h"
#include "Scene.h"
#include "Vector.h"
#include "Math.h"
using namespace std;

LambertianMaterial::LambertianMaterial(const Color& color, float Kd, float Ka, float Ks, int n, bool isReflective)
  :color(color), Kd(Kd), Ka(Ka), Ks(Ks), n(n), isReflective(isReflective)
{
}

LambertianMaterial::~LambertianMaterial()
{
}

bool LambertianMaterial::getReflective() const {
    return isReflective;
}

Color LambertianMaterial::getColor() const {
    return color;
}

float LambertianMaterial::getKs() const {
    return Ks;
}

Color LambertianMaterial::sample_f(const Vector& normal, const Vector& wo, Vector& wi, float& pdf) const {
    // source: Ray Tracing from the Ground Up, pg 548

    Vector v = Cross(Vector(0.0034, 1.0, 0.0071), normal);
    v.normalize();
    Vector u = Cross(v, normal);

    double rand1 = double(rand()) / double(RAND_MAX);
    double rand2 = double(rand()) / double(RAND_MAX);
    double rand3 = double(rand()) / double(RAND_MAX);

    Point sp = Point(rand1, rand2, rand3);
    wi = sp.x() * u + sp.y() * v + sp.z() * normal;
    wi.normalize();
    pdf = Dot(normal, wi) * (1.0 / M_PI);

    return color * Kd * (1.0 / M_PI);
}

//Color LambertianMaterial::shade(const RenderContext& context,
//                        const Ray& ray, const HitRecord& hit, const Color&, int depth) const
//{
//    const Scene* scene = context.getScene();
//    const vector<Light*>& lights = scene->getLights();
//    Point hitpos = ray.origin() + ray.direction() * hit.minT();
//
//    // compute normal at hitpos
//    Vector normal;
//    hit.getPrimitive()->normal(normal, context, hitpos, ray, hit);
//    normal.normalize();
//    double costheta = Dot(normal, ray.direction());
//    if (costheta > 0)
//        normal = -normal;
//
//    const Object* world = scene->getObject();
//
//    //Color light = scene->getAmbient() * Ka;
//    Color light = Color(0, 0, 0);
//
//#if 0
//    for (vector<Light*>::const_iterator iter = lights.begin(); iter != lights.end(); iter++) {
//#else
//    Light* const* begin = &lights[0];
//    Light* const* end = &lights[0] + lights.size();
//    while (begin != end) {
//#endif
//        if (depth == 0)
//            light += (*begin)->areaShade(ray, normal);
//
//        // compute direction vector from light to hitpos
//        //Color light_color;
//        //Vector light_direction;
//        //double dist = (*begin++)->getLight(light_color, light_direction, context, hitpos);
//        //light_direction.normalize();
//        //double cosphi = Dot(normal, light_direction);
//
//        //Vector reflect_direction = 2 * Dot(normal, light_direction) * normal - light_direction;
//        //reflect_direction.normalize();
//
//        //// multiply ray direction by -1 to use direction from hitpos back to eye
//        //double cosbet = Dot(reflect_direction, ray.direction() * -1);
//
//        Vector wi, wo = -ray.direction();
//        float pdf;
//        Color f = this->sample_f(normal, wo, wi, pdf);
//        float ndotwi = Dot(normal, wi);
//        Ray reflected_ray(hitpos, wi);
//        light += 
//
//        if (cosphi > 0) {
//            // Cast shadow rays...
//            HitRecord shadowhit(dist);
//            Ray shadowray(hitpos, light_direction);
//            world->intersect(shadowhit, context, shadowray);
//
//            if (!shadowhit.getPrimitive()) {
//                // No shadows...
//                if (cosbet > 0) {
//                    light += light_color * (Kd * cosphi + Ks * pow(cosbet, n));
//                }
//                else {
//                    light += light_color * (Kd * cosphi);
//                }
//            }
//        }
//    }
//
//    return light * color;
//}

Color LambertianMaterial::shade(const RenderContext& context,
                               const Ray& ray, const HitRecord& hit, const Color&, int) const
{
    const Scene* scene = context.getScene();
    const vector<Light*>& lights = scene->getLights();
    Point hitpos = ray.origin() + ray.direction() * hit.minT();

    // compute normal at hitpos
    Vector normal;
    hit.getPrimitive()->normal(normal, context, hitpos, ray, hit);
    normal.normalize();
    double costheta = Dot(normal, ray.direction());
    if(costheta > 0)
        normal = -normal;

    const Object* world = scene->getObject();

    Color light = scene->getAmbient()*Ka;

#if 0
    for(vector<Light*>::const_iterator iter = lights.begin(); iter != lights.end(); iter++){
#else
    Light*const* begin = &lights[0];
    Light*const* end = &lights[0] + lights.size();
    while(begin != end){
#endif
        // compute direction vector from light to hitpos
        Color light_color;
        Vector light_direction;
        double dist = (*begin++)->getLight(light_color, light_direction, context, hitpos);
        light_direction.normalize();
        double cosphi = Dot(normal, light_direction);

        Vector reflect_direction = 2 * Dot(normal, light_direction) * normal - light_direction;
        reflect_direction.normalize();

        // multiply ray direction by -1 to use direction from hitpos back to eye
        double cosbet = Dot(reflect_direction, ray.direction() * -1);

        if(cosphi > 0){
            // Cast shadow rays...
            HitRecord shadowhit(dist);
            Ray shadowray(hitpos, light_direction);
            world->intersect(shadowhit, context, shadowray);

            if (!shadowhit.getPrimitive()) {
                // No shadows...
                if (cosbet > 0) {
                    light += light_color * (Kd * cosphi + Ks * pow(cosbet, n));
                }
                else {
                    light += light_color * (Kd * cosphi);
                }
            }
        }
    }

    return light * color;
}

// source for the following functions: https://raytracing.github.io/books/RayTracingTheRestOfYourLife.html#importancesamplingmaterials

Vector random_cosine_direction() {
    /*float r1 = random_double();
    float r2 = random_double();*/
    double r1 = double(rand()) / double(RAND_MAX + 1);
    double r2 = double(rand()) / double(RAND_MAX + 1);
    float z = sqrt(1 - r2);
    float phi = 2 * M_PI * r1;
    float x = cos(phi) * sqrt(r2);
    float y = sin(phi) * sqrt(r2);
    return Vector(x, y, z);
}

bool LambertianMaterial::scatter(Point hitpos, Vector normal, Ray& scattered, float& pdf) const {
    // generate 3 random numbers between -1 and 1
    /*double rand1 = double(rand()) / double(RAND_MAX);
    rand1 = -1.0 + 2.0 * rand1;
    double rand2 = double(rand()) / double(RAND_MAX);
    rand2 = -1.0 + 2.0 * rand2;
    double rand3 = double(rand()) / double(RAND_MAX);
    rand3 = -1.0 + 2.0 * rand3;
    Vector randVec = Vector(rand1, rand2, rand3);

    Vector point = Vector(hitpos);
    Vector target = normal + point + randVec;
    Vector dir = Vector(target - point);
    dir.normalize();
    scattered = Ray(hitpos, dir);
    pdf = Dot(normal, scattered.direction()) / M_PI;*/

    onb uvw;
    uvw.build_from_w(normal);
    Vector direction = uvw.local(random_cosine_direction());
    direction.normalize();
    scattered = Ray(hitpos, direction);
    pdf = Dot(uvw.w(), scattered.direction()) / M_PI;

    return true;
}

float LambertianMaterial::scattering_pdf(const Vector& normal, const Ray& scattered) const {
    Vector scatteredDir = scattered.direction();
    scatteredDir.normalize();
    float cosine = Dot(normal, scatteredDir);
    if (cosine < 0)
        return 0;
    return cosine / M_PI;
}


void onb::build_from_w(Vector& n) {
    n.normalize();
    axis[2] = n;
    Vector a;
    if (fabs(w().x()) > 0.9)
        a = Vector(0, 1, 0);
    else
        a = Vector(1, 0, 0);
    axis[1] = Cross(w(), a);
    axis[1].normalize();
    axis[0] = Cross(w(), v());
}
