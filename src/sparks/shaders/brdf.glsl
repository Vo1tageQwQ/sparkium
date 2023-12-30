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

vec3 PositionAlongNormal(vec3 normal, vec3 dir) {
  vec3 u, v;
  OrthonormalBasis(normal, u, v);
  return dir.x * u + dir.y * v + dir.z * normal;
}

void SampleCosineHemisphere(vec3 normal, out vec3 dir, out float pdf) {
  float r1 = RandomFloat() * 2.0 * PI;
  float r2 = RandomFloat();

  dir = vec3(cos(r1) * sqrt(r2), sin(r1) * sqrt(r2), sqrt(1.0 - r2));
  pdf = dir.z / PI;
  dir = PositionAlongNormal(normal, dir);
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

vec3 SampleGTR2Aniso(vec3 normal, vec3 tangent, vec3 cotangent, float ax, float ay, float r1, float r2) {
  r1 *= 2.0 * PI;
  return normalize(normal + sqrt(r2 / (1.0 - r2)) * (ax * cos(r1) * tangent + ay * sin(r1) * cotangent));
}

vec3 SampleGTR1(vec3 normal, vec3 tangent, vec3 cotangent, float a, float r1, float r2) {
  float aSqr = a * a;
  float cosThetaSqr = (1.0 - pow(aSqr, r1)) / (1.0 - aSqr);
  float cosTheta = sqrt(cosThetaSqr);
  float sinTheta = sqrt(1.0 - cosThetaSqr);
  float phi = 2.0 * PI * r2;
  return cosTheta * normal + sinTheta * (cos(phi) * tangent + sin(phi) * cotangent);
}

// See https://shihchinw.github.io/2015/07/implementing-disney-principled-brdf-in-arnold.html
void SampleBRDFPrincipled(vec3 normal, vec3 tangent, vec3 base_color, DisneyParams params, vec3 V, out vec3 L, out vec3 eval, out float pdf) {
  vec3 cotangent = cross(normal, tangent);
  float aspect = sqrt(1.0 - params.anisotropic * 0.9);
  float a = params.roughness * params.roughness;
  float ax = max(0.001, a / aspect);
  float ay = max(0.001, a * aspect);

  float r1 = RandomFloat(), r2 = RandomFloat();
  float gtr2Weight = 1.0 / (1.0 + params.clearcoat);
  if (r1 < gtr2Weight) {
    r1 /= gtr2Weight;
    L = SampleGTR2Aniso(normal, tangent, cotangent, ax, ay, r1, r2);
  } else {
    r1 = (r1 - gtr2Weight) / (1.0 - gtr2Weight);
    L = SampleGTR1(normal, tangent, cotangent, a, r1, r2);
  }

  vec3 m = normalize(L + V);
  float LdotM = dot(L, m);
  float NdotM = dot(normal, m);
  if (NdotM <= 0.0) {
    pdf = 0.0;
    return;
  }
  float NdotMSqr = NdotM * NdotM;
  float clearcoatWeight = params.clearcoat / (1.0 + params.clearcoat);
  float d = mix(GTR2_aniso(NdotM, dot(tangent, m), dot(cotangent, m), ax, ay), GTR1(NdotM, a), clearcoatWeight);
  pdf = d * abs(NdotM) / (4.0 * LdotM);

  eval = BRDF_Disney(L, V, normal, tangent, cotangent, base_color, params);
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
