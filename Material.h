
#ifndef Material_h
#define Material_h

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

  /*virtual bool scatter(const Ray& ray, const HitRecord& hit, Vector albedo, Ray& scattered, float& pdf) const {
      return false;
  }*/

 private:
  Material(const Material&);
  Material& operator=(const Material&);
};

#endif

