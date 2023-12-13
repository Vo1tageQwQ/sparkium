#include "sparks/renderer/path_tracer.h"
#include <vector>
#include <tuple>
#include <random>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "sparks/util/util.h"
#include "sparks/assets/brdf.h"

namespace sparks {
PathTracer::PathTracer(const RendererSettings *render_settings,
                       const Scene *scene) {
  render_settings_ = render_settings;
  scene_ = scene;
  rd = std::mt19937(std::random_device()());
  uniform = std::uniform_real_distribution<float>(0, 1);
}

void PathTracer::SampleFromLight(glm::vec3 &lpos, glm::vec3 &lnorm, float &area) {
  const std::vector<Entity> &entities = scene_->GetEntities();
  std::vector<float> areaAccum;
  float totalArea = 0.0f;
  std::vector<std::tuple<glm::vec3, glm::vec3, glm::vec3>> triangles;
  for (int i = 0; i < entities.size(); i++) {
    if (entities[i].GetMaterial().material_type != 4) {
      continue;
    }
    const std::vector<Vertex> vertices = entities[i].GetModel()->GetVertices();
    const std::vector<uint32_t> indices = entities[i].GetModel()->GetIndices();
    glm::mat4 transform = entities[i].GetTransformMatrix();
    for (int j = 0; j < indices.size(); j+=3) {
      glm::vec4 v0(vertices[indices[j]].position, 1.0f);
      glm::vec4 v1(vertices[indices[j + 1]].position, 1.0f);
      glm::vec4 v2(vertices[indices[j + 2]].position, 1.0f);
      glm::vec3 u0 = glm::vec3(transform * v0);
      glm::vec3 u1 = glm::vec3(transform * v1);
      glm::vec3 u2 = glm::vec3(transform * v2);
      float deltaArea = glm::length(glm::cross(u2 - u0, u1 - u0)) / 2.0f;
      totalArea += deltaArea;
      areaAccum.emplace_back(totalArea);
      triangles.emplace_back(std::make_tuple(u0, u1, u2));
    }
  }
  if(totalArea == 0.0f) {
    lpos = glm::vec3{0.0f};
    lnorm = glm::vec3{0.0f};
    area = 0.0f;
    return;
  }
  float areaSample = totalArea * uniform(rd);
  int target = std::lower_bound(areaAccum.begin(), areaAccum.end(), areaSample) - areaAccum.begin();
  float r1 = uniform(rd), r2 = uniform(rd);
  glm::vec3 v[3];
  std::tie(v[0], v[1], v[2]) = triangles[target];
  lpos = (1.0f - glm::sqrt(r1)) * v[0] + (glm::sqrt(r1) * (1.0f - r2)) * v[1] + r2 * glm::sqrt(r1) * v[2];
  lnorm = glm::normalize(glm::cross(v[2] - v[0], v[1] - v[0]));\
  area = totalArea;
}

glm::vec3 PathTracer::HemisphereSample(float u, float v, glm::vec3 norm) {
  glm::vec3 dir;
  float r, phi;
  u = 2.0f * u - 1.0f;
  v = 2.0f * v - 1.0f;
  if (u > -v) {
    if (u > v) {
      r = u;
      phi = (PI / 4.0f) * (v / u);
    }
    else {
      r = v;
      phi = (PI / 4.0f) * (2.0f - (u / v));
    }
  }
  else {
    if (u < v) {
      r = -u;
      phi = (PI / 4.0f) * (4.0f + (v / u));
    }
    else {
      r = -v;
      if (v != 0.0f) {
        phi = (PI / 4.0f) * (6.0f - (u / v));
      }
      else {
        phi = 0.0f;
      }
    }
  }
  dir.x = r * glm::cos(phi);
  dir.y = r * glm::sin(phi);
  dir.z = glm::sqrt(1.0f - r * r);
  glm::vec3 axis = glm::normalize(glm::vec3(-norm.y, norm.x, 0.0f));
  float c = norm.z;
  float s = glm::sqrt(1.0f - c * c);
  glm::mat3 tm{c + axis.x * axis.x * (1.0f - c), axis.y * axis.x * (1.0f - c) + axis.z * s, axis.z * axis.x * (1.0f - c) - axis.y * s,
              axis.x * axis.y * (1.0f - c) - axis.z * s, c + axis.y * axis.y * (1.0f - c), axis.z * axis.y * (1.0f - c) + axis.x * s,
              axis.x * axis.z * (1.0f - c) + axis.y * s, axis.y * axis.z * (1.0f - c) - axis.x * s, c + axis.z * axis.z * (1.0f - c)};
  glm::vec3 samp = tm * dir;
  return samp;
}


glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample) {
  glm::vec3 throughput{1.0f};
  glm::vec3 radiance{0.5f};
  HitRecord hit_record;
  int max_bounce = render_settings_->num_bounces;
  float t0 = scene_->TraceRay(origin, direction, 1e-4f, 1e7f, &hit_record);
  if (t0 < 0.0f) {
    radiance += throughput * scene_->GetEnvmapMajorColor();
    return radiance;
  }
  for (int i = 0; i < max_bounce; i++) {
    glm::vec3 pos = hit_record.position, norm = hit_record.geometry_normal;
    Material mater = scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
    if (mater.material_type == 4) {
      radiance += throughput * mater.emission * mater.emission_strength;
      break;
    }
    glm::vec3 brdf = mater.albedo_color;
    //Constribution from the light source.
    glm::vec3 lpos, lnorm;
    float area;
    SampleFromLight(lpos, lnorm, area);
    if (area != 0.0f) {
      glm::vec3 ldir = glm::normalize(lpos - pos), lrad;
      HitRecord lhit;
      float lt = scene_->TraceRay(pos, ldir, 1e-4f, 1e7f, &lhit);
      if (glm::length(pos + lt * ldir - lpos) < 0.1f) {
        Material lmater = scene_->GetEntity(lhit.hit_entity_id).GetMaterial();
        lrad = lmater.emission * lmater.emission_strength * brdf * area / glm::dot(lpos - pos, lpos - pos) * glm::dot(norm, ldir);
      }
      radiance += lrad * throughput;
    }
    //Constribution from other reflectors.
    float RRTest = uniform(rd);
    if (RRTest >= RRProb) {
      break;
    }
    glm::vec3 ndir = HemisphereSample(uniform(rd), uniform(rd), norm);
    float nt = scene_->TraceRay(pos, ndir, 1e-4f, 1e7f, &hit_record);
    throughput *= brdf * 2.0f * PI * glm::dot(norm, ndir) / RRTest;
    if (nt < 0.0f) {
      radiance += throughput * scene_->GetEnvmapMajorColor();
      break;
    }
  }
  return radiance;
}
}  // namespace sparks
