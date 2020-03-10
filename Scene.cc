 
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
#include "Vector.h"
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <cmath>

using namespace std;

#define ALBEDO 0.18
#define BIAS 0.0001


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

            int dim = (int)sqrt((float)numSamples);
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    // get x and y values within grid cell to sample
                    double randNum = double(rand()) / double(RAND_MAX);
                    double sampleX = double(j) - 0.5 + (double(l) + randNum) / double(dim);
                    randNum = double(rand()) / double(RAND_MAX);
                    double sampleY = double(i) - 0.5 + (double(k) + randNum) / double(dim);
                    double x = xmin + sampleX * dx;
                    double y = ymin + sampleY * dy;

                    camera->makeRay(ray, context, x, y);
                    result += traceRay(context, ray, atten, 0);
                }
            }

            result /= float(numSamples);

            image->set(j, i, result);
        }
    }
}

Color Scene::traceRay(const RenderContext& context, const Ray& ray, const Color& atten, int depth) const
{
    // algorithm: https://www.scratchapixel.com/code.php?id=34&origin=/lessons/3d-basic-rendering/global-illumination-path-tracing

    Color result(0, 0, 0), indirect(0, 0, 0), direct(0, 0, 0);
    HitRecord hit(DBL_MAX);
    object->intersect(hit, context, ray);

    if (hit.getPrimitive()) {
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

        Ray scattered;
        float pdf;

        if (depth < maxRayDepth && matl->scatter(hitpos, normal, scattered, pdf)) {
            indirect = traceRay(context, scattered, atten, depth + 1) * ALBEDO
                * matl->scattering_pdf(normal, scattered) / pdf;
        }

        result = (direct + indirect);
    }
    else {
        background->getBackgroundColor(result, context, ray);
    }

    return result;
}
