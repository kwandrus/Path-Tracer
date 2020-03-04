
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


  virtual float scattering_pdf(const Vector& normal, const Ray& scattered) const;
  virtual bool scatter(Point hitpos, Vector normal, Ray& scattered, float& pdf) const;
  Color sample_f(const Vector& normal, const Vector& wo, Vector& wi, float& pdf) const;


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

class onb
{
public:
    onb() {}
    inline Vector operator[](int i) const { return axis[i]; }
    Vector u() const { return axis[0]; }
    Vector v() const { return axis[1]; }
    Vector w() const { return axis[2]; }
    Vector local(float a, float b, float c) const { return a * u() + b * v() + c * w(); }
    Vector local(const Vector& a) const { return a.x() * u() + a.y() * v() + a.z() * w(); }
    void build_from_w(Vector&);
    Vector axis[3];
};

#endif
