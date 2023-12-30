#ifndef SPARKIUM_MATERIAL_GLSL
#define SPARKIUM_MATERIAL_GLSL

#extension GL_GOOGLE_include_directive : require

// clang-format off
#include "disney.glsl"
// clang-format on

struct Material {
  vec3 albedo_color;
  int albedo_texture_id;
  vec3 emission;
  float emission_strength;
  float alpha;
  uint material_type;

  DisneyParams disney_params;
};

#define MATERIAL_TYPE_LAMBERTIAN 0
#define MATERIAL_TYPE_SPECULAR 1
#define MATERIAL_TYPE_TRANSMISSIVE 2
#define MATERIAL_TYPE_PRINCIPLED 3
#define MATERIAL_TYPE_EMISSION 4

#endif // SPARKIUM_MATERIAL_GLSL
