#ifndef SPARKIUM_BRDF_GLSL
#define SPARKIUM_BRDF_GLSL

#extension GL_GOOGLE_include_directive : require

// clang-format off
#include "constants.glsl"
#include "disney_impl.glsl"
#include "hit_record.glsl"
#include "material.glsl"
#include "random.glsl"
// clang-format on

void OrthonormalBasis(vec3 n, out vec3 u, out vec3 v) {
  if (abs(n.z) > 0.58) {
    // u = normalize(cross(vec3(1, 0, 0), n));
    u = normalize(vec3(n.z, 0, -n.x));
  } else {
    // u = normalize(cross(vec3(0, 0, 1), n));
    u = normalize(vec3(-n.y, n.x, 0));
  }
  v = normalize(cross(n, u));
}

void SampleCosineHemisphere(vec3 normal, out vec3 dir, out float pdf) {
  vec3 u, v;
  OrthonormalBasis(normal, u, v);

  float r1 = RandomFloat() * 2.0 * PI;
  float r2 = RandomFloat();

  dir = vec3(cos(r1) * sqrt(r2), sin(r1) * sqrt(r2), sqrt(1.0 - r2));
  pdf = dir.z / PI;
  dir = dir.x * u + dir.y * v + dir.z * normal;
}

void SampleBRDFLambertian(vec3 normal, out vec3 L, out float eval, out float pdf) {
  SampleCosineHemisphere(normal, L, pdf);
  eval = pdf;
}

void SampleBRDFSpecular(vec3 normal, vec3 V, out vec3 L, out float eval, out float pdf) {
  float cosTheta = dot(normal, V);
  if (cosTheta <= 0.0) {
    eval = 0.0;
    pdf = 0.0;
  } else {
    L = 2.0 * cosTheta * normal - V;
    pdf = INF;
    eval = INF;
  }
}

void SampleBRDFPrincipled(vec3 normal, vec3 tangent, vec3 base_color, DisneyParams disney_params, vec3 V, out vec3 L, out vec3 eval, out float pdf) {
  // TODO: Use importance sampling.
  SampleCosineHemisphere(normal, L, pdf);
  vec3 cotangent = cross(normal, tangent);
  eval = BRDF_Disney(L, V, normal, tangent, cotangent, base_color, disney_params);
}

void SampleBRDF(HitRecord hit_record, vec3 V, out vec3 L, out vec3 eval, out float pdf) {
  float brdf_value;
  switch (hit_record.material_type) {
    case MATERIAL_TYPE_LAMBERTIAN:
      SampleBRDFLambertian(hit_record.normal, L, brdf_value, pdf);
      eval = brdf_value * hit_record.base_color;
      break;
    case MATERIAL_TYPE_SPECULAR:
      SampleBRDFSpecular(hit_record.normal, V, L, brdf_value, pdf);
      eval = brdf_value * hit_record.base_color;
      break;
    case MATERIAL_TYPE_PRINCIPLED:
      SampleBRDFPrincipled(hit_record.normal, hit_record.tangent, hit_record.base_color, hit_record.disney_params, V, L, eval, pdf);
      break;
    default:
      eval = vec3(0.0);
      pdf = 0.0;
      return;
  }
  if (dot(L, hit_record.geometry_normal) <= 0.0) {
    eval = vec3(0.0);
    pdf = 0.0;
  }
}

#endif // SPARKIUM_BRDF_GLSL
