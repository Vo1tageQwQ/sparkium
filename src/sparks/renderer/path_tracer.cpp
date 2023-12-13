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


glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample) {
  glm::vec3 throughput{1.0f};
  glm::vec3 radiance{0.0f};
  HitRecord hit_record;
  int max_bounce = render_settings_->num_bounces;
  float t0 = scene_->TraceRay(origin, direction, 1e-4f, 1e7f, &hit_record);
  if (t0 < 0.0f) {
    return radiance;
  }
  if (scene_->GetEntity(hit_record.hit_entity_id).GetMaterial().material_type == 4) {
    Material lmater = scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
    radiance = lmater.emission * lmater.emission_strength / glm::dot(hit_record.position - origin, hit_record.position - origin);
    return radiance;
  }
  for (int i = 0; i < max_bounce; i++) {
    //Constribution from the light source.
    glm::vec3 lpos, lnorm;
    float area;
    SampleFromLight(lpos, lnorm, area);
    if (area != 0.0f) {
      glm::vec3 ldir = glm::normalize(lpos - origin);
      HitRecord lhit;
      float lt = scene_->TraceRay(origin, ldir, 1e-4f, 1e7f, &lhit);
      if (glm::length(origin + lt * ldir - lpos) < 0.1f) {
        Material lmater = scene_->GetEntity(lhit.hit_entity_id).GetMaterial();
        glm::vec3 lrad = lmater.emission * lmater.emission_strength * area ;
      }
    }
    //Constribution from other reflectors.
    float RRTest = uniform(rd);
    if (RRTest >= RRProb) {
      break;
    }

  }
  return radiance;
}
}  // namespace sparks
