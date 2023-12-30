#ifndef SPARKIUM_HIT_RECORD_GLSL
#define SPARKIUM_HIT_RECORD_GLSL

struct HitRecord {
  int hit_entity_id;
  vec3 position;
  vec3 normal;
  vec3 geometry_normal;
  vec3 tangent;
  vec2 tex_coord;
  bool front_face;

  vec3 base_color;
  vec3 emission;
  float emission_strength;
  float alpha;
  uint material_type;
};

#endif // SPARKIUM_HIT_RECORD_GLSL
