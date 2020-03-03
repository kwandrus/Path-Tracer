 
#include "Scene.h"
#include "Background.h"
#include "Camera.h"
#include "HitRecord.h"
#include "Image.h"
#include "Light.h"
#include "Material.h"
#include "Object.h"
#include "Ray.h"
#include "RenderContext.h"
#include "Primitive.h"
#include <float.h>
#include <iostream>
#include <stdlib.h>

using namespace std;

#define ALBEDO 0.18
#define BIAS 0.0001

Vector uniformSampleHemisphere(const float& r1, const float& r2);
void createCoordinateSystem(const Vector& N, Vector& Nt, Vector& Nb);

Scene::Scene()
{
    object = 0;
    background = 0;
    camera = 0;
    ambient = Color(0, 0, 0);
    image = 0;
    minAttenuation = 0;
    numSamples = 0;
}

Scene::~Scene()
{
    delete object;
    delete background;
    delete camera;
    delete image;
    for(int i = 0; i < static_cast<int>(lights.size()); i++){
        Light* light = lights[i];
        delete light;
    }
}

void Scene::preprocess()
{
    background->preprocess();
    for(int i = 0; i < static_cast<int>(lights.size()); i++){
        Light* light = lights[i];
        light->preprocess();
    }

    double aspect_ratio = image->aspect_ratio();
    camera->preprocess(aspect_ratio);
    object->preprocess();
}

int max(int n1, int n2) 
{
    if (n1 > n2)
        return n1;
    return n2;
}

void Scene::render()
{
    if(!object || !background || !camera || !image){
    cerr << "Incomplete scene, cannot render!\n";
    exit(1);
    }
    int xres = image->getXresolution();
    int yres = image->getYresolution();
    RenderContext context(this);
    double dx = 2. / xres;
    double xmin = -1. + dx / 2.;
    double dy = 2. / yres;
    double ymin = -1. + dy / 2.;
    Color atten(1, 1, 1);

    for (int i = 0; i < yres; i++) {
        for (int j = 0; j < xres; j++) {
            Ray ray;
            Color result(0, 0, 0);

            double x = xmin + j * dx;
            double y = ymin + i * dy;
            camera->makeRay(ray, context, x, y);
            result += traceRay(context, ray, atten, 0);

            image->set(j, i, result);
        }
    }
}

Color Scene::traceRay(const RenderContext& context, const Ray& ray, const Color& atten, int depth) const
{
    // code based off of: https://www.scratchapixel.com/code.php?id=34&origin=/lessons/3d-basic-rendering/global-illumination-path-tracing
    
    Color result(0, 0, 0), indirect(0, 0, 0), direct(0, 0, 0);
    HitRecord hit(DBL_MAX);
    object->intersect(hit, context, ray);

    if (depth <= maxRayDepth && hit.getPrimitive()) {
        // Ray hit something...
        const Material* matl = hit.getMaterial();

        // compute DIRECT light
        direct = matl->shade(context, ray, hit, atten, depth);

        // compute INDIRECT light
        // intersection info
        Point hitpos = ray.origin() + ray.direction() * hit.minT();
        Vector normal;
        hit.getPrimitive()->normal(normal, context, hitpos, ray, hit);
        normal.normalize();

        Vector Nt, Nb;
        createCoordinateSystem(normal, Nt, Nb);
        float pdf = float(1) / (float(2) * M_PI);
        for (int i = 0; i < numSamples; i++) {
            // create random numbers between [0, 1)
            double randNum1 = double(rand()) / double(RAND_MAX + 1);
            double randNum2 = double(rand()) / double(RAND_MAX + 1);

            Vector sample = uniformSampleHemisphere(randNum1, randNum2);
            Vector sampleWorld(
                sample.x * Nb.x + sample.y * normal.x + sample.z * Nt.x,
                sample.x * Nb.y + sample.y * normal.y + sample.z * Nt.y,
                sample.x * Nb.z + sample.y * normal.z + sample.z * Nt.z);

            // create ray
            Ray nextRay(hitpos + sampleWorld * BIAS, sampleWorld);

            // don't forget to divide by PDF and multiply by cos(theta)
            // randNum1 == cos(theta) from uniformSampleHemisphere()
            indirect += (traceRay(context, nextRay, atten, depth + 1) * randNum1) / pdf;
        }
        // divide by N
        indirect /= float(numSamples);

        result = (direct / M_PI + indirect * 2.0) * ALBEDO;     
    }
    else {
        background->getBackgroundColor(result, context, ray);
    }

    return result;
}

void createCoordinateSystem(const Vector& N, Vector& Nt, Vector& Nb)
{
    // source: https://www.scratchapixel.com/code.php?id=34&origin=/lessons/3d-basic-rendering/global-illumination-path-tracing

    if (std::fabs(N.x) > std::fabs(N.y))
        Nt = Vector(N.z, 0, -N.x) / sqrtf(N.x * N.x + N.z * N.z);
    else
        Nt = Vector(0, -N.z, N.y) / sqrtf(N.y * N.y + N.z * N.z);
    Nb = Cross(N, Nt);
}

Vector uniformSampleHemisphere(const float& r1, const float& r2)
{
    // source: https://www.scratchapixel.com/code.php?id=34&origin=/lessons/3d-basic-rendering/global-illumination-path-tracing
    
    float sinTheta = sqrtf(1 - r1 * r1);
    float phi = 2 * M_PI * r2;
    float x = sinTheta * cosf(phi);
    float z = sinTheta * sinf(phi);
    return Vector(x, r1, z);
}
