
#ifndef AreaLight_h
#define AreaLight_h

#include "Light.h"
#include "Point.h"
#include "Color.h"
#include "Ray.h"

class AreaLight : public Light {
 public:
  AreaLight(const Point& position, const Color& color, const Vector& a, const Vector& b);
  virtual ~AreaLight();

  virtual void preprocess();
  virtual double getLight(Color& light_color, Vector& light_direction,
                          const RenderContext& context, const Point& pos) const;
  Color areaShade(const Ray& ray, Vector normal) const;

 private:
  AreaLight(const AreaLight&);
  AreaLight& operator=(const AreaLight&);

  Point cornerPos;
  Color color;
  Vector a, b; // represent dimensions of rectangle  

  float ls; // radiance scaling factor
};

#endif

