
#ifndef Material_h
#define Material_h

#include "Point.h"

class Color;
class HitRecord;
class Ray;
class RenderContext;

class Material {
 public:
  Material();
  virtual ~Material();

  virtual void preprocess();
  virtual Color shade(const RenderContext& context, const Ray& ray,
                     const HitRecord& hit, const Color& atten, int depth) const = 0;
  virtual bool getReflective() const = 0;
  virtual Color getColor() const = 0;
  virtual float getKs() const = 0;

  virtual bool scatter(Point hitpos, Vector normal, Ray& scattered, float& pdf) const {
      return false;
  }
  virtual float scattering_pdf(const Vector& normal, const Ray& scattered) const {
      return 0;
  }
  virtual Vector emitted(float u, float v, const Vector& p) const {
      return Vector(0, 0, 0);
  }

 private:
  Material(const Material&);
  Material& operator=(const Material&);
};

#endif

