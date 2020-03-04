
#ifndef LambertianMaterial_h
#define LambertianMaterial_h

#include "Material.h"
#include "Color.h"
#include "Point.h"
#include "Ray.h"

class LambertianMaterial : public Material {
 public:
  LambertianMaterial(const Color& color, float Kd, float Ka, float Ks, int n, bool isReflective);
  virtual ~LambertianMaterial();

  virtual Color shade(const RenderContext& context, const Ray& ray,
                     const HitRecord& hit, const Color& atten, int depth) const;
  virtual bool getReflective() const;
  virtual Color getColor() const;
  virtual float getKs() const;


  // following 2 functions source: https://raytracing.github.io/books/RayTracingTheRestOfYourLife.html#importancesamplingmaterials
  
  virtual bool scatter(Point hitpos, Vector normal, Ray& scattered, float& pdf) const {
      // generate 3 random numbers between -1 and 1
      double rand1 = double(rand()) / double(RAND_MAX);
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
      pdf = Dot(normal, scattered.direction()) / M_PI;
      return true;
  }

  virtual float scattering_pdf(const Vector& normal, const Ray& scattered) const {
      Vector scatteredDir = scattered.direction();
      scatteredDir.normalize();
      float cosine = Dot(normal, scatteredDir);
      if (cosine < 0)
          return 0;
      return cosine / M_PI;
  }


 private:
  LambertianMaterial(const LambertianMaterial&);
  LambertianMaterial& operator=(const LambertianMaterial&);

  Color color;
  float Kd;
  float Ka;
  float Ks;
  int n;
  bool isReflective;
};

#endif
