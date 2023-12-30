#ifndef SPARKIUM_HIT_RECORD_GLSL
#define SPARKIUM_HIT_RECORD_GLSL

#extension GL_GOOGLE_include_directive : require

// clang-format off
#include "disney.glsl"
// clang-format on

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

  DisneyParams disney_params;
};

#endif // SPARKIUM_HIT_RECORD_GLSL
