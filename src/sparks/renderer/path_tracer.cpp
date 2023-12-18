#include <random>
#include <tuple>
#include <vector>

#include <glm/glm.hpp>

#include "sparks/assets/brdf.h"
#include "sparks/renderer/path_tracer.h"
#include "sparks/util/util.h"

namespace sparks {
PathTracer::PathTracer(const RendererSettings *render_settings,
                       const Scene *scene)
    : render_settings_(render_settings),
      scene_(scene),
      rd(std::random_device()()) {
}

void PathTracer::SampleFromLight(glm::vec3 &lpos,
                                 glm::vec3 &lnorm,
                                 float &area) const {
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
    for (int j = 0; j < indices.size(); j += 3) {
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
  if (totalArea == 0.0f) {
    lpos = glm::vec3{0.0f};
    lnorm = glm::vec3{0.0f};
    area = 0.0f;
    return;
  }
  float areaSample = totalArea * random();
  int target =
      std::lower_bound(areaAccum.begin(), areaAccum.end(), areaSample) -
      areaAccum.begin();
  float r1 = random(), r2 = random();
  glm::vec3 v[3];
  std::tie(v[0], v[1], v[2]) = triangles[target];
  lpos = (1.0f - glm::sqrt(r1)) * v[0] + (glm::sqrt(r1) * (1.0f - r2)) * v[1] +
         r2 * glm::sqrt(r1) * v[2];
  lnorm = glm::normalize(glm::cross(v[2] - v[0], v[1] - v[0]));
  area = totalArea;
}

glm::vec3 PathTracer::HemisphereSample(float u, float v, glm::vec3 norm) const {
  glm::vec3 dir;
  float r, phi;
  u = 2.0f * u - 1.0f;
  v = 2.0f * v - 1.0f;
  if (u > -v) {
    if (u > v) {
      r = u;
      phi = (PI / 4.0f) * (v / u);
    } else {
      r = v;
      phi = (PI / 4.0f) * (2.0f - (u / v));
    }
  } else {
    if (u < v) {
      r = -u;
      phi = (PI / 4.0f) * (4.0f + (v / u));
    } else {
      r = -v;
      if (v != 0.0f) {
        phi = (PI / 4.0f) * (6.0f - (u / v));
      } else {
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
  glm::mat3 tm{c + axis.x * axis.x * (1.0f - c),
               axis.y * axis.x * (1.0f - c) + axis.z * s,
               axis.z * axis.x * (1.0f - c) - axis.y * s,
               axis.x * axis.y * (1.0f - c) - axis.z * s,
               c + axis.y * axis.y * (1.0f - c),
               axis.z * axis.y * (1.0f - c) + axis.x * s,
               axis.x * axis.z * (1.0f - c) + axis.y * s,
               axis.y * axis.z * (1.0f - c) - axis.x * s,
               c + axis.z * axis.z * (1.0f - c)};
  glm::vec3 samp = tm * dir;
  return samp;
}

glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample) const {
  glm::vec3 throughput{1.0f};
  glm::vec3 radiance{0.0f};
  glm::vec3 pos = origin, dir = direction, norm;
  float tt;
  HitRecord hit_record;
  int max_bounce = render_settings_->num_bounces;
  for (int i = 0; i <= max_bounce; i++) {
    tt = scene_->TraceRay(pos, dir, 1e-4f, 1e7f, &hit_record);
    if (tt < 0.0f) {
      radiance += throughput * scene_->GetEnvmapMajorColor();
      break;
    }
    Material mater = scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
    glm::vec3 brdf = mater.albedo_color;
    pos = hit_record.position, norm = hit_record.geometry_normal;
    pos += norm * 1e-3f;
    if (mater.material_type == MATERIAL_TYPE_EMISSION) {
      const Entity &lentity = scene_->GetEntity(hit_record.hit_entity_id);
      const std::vector<Vertex> vertices = lentity.GetModel()->GetVertices();
      const std::vector<uint32_t> indices = lentity.GetModel()->GetIndices();
      float lArea = 0.0f;
      glm::mat4 transform = lentity.GetTransformMatrix();
      for (int j = 0; j < indices.size(); j += 3) {
        glm::vec4 v0(vertices[indices[j]].position, 1.0f);
        glm::vec4 v1(vertices[indices[j + 1]].position, 1.0f);
        glm::vec4 v2(vertices[indices[j + 2]].position, 1.0f);
        glm::vec3 u0 = glm::vec3(transform * v0);
        glm::vec3 u1 = glm::vec3(transform * v1);
        glm::vec3 u2 = glm::vec3(transform * v2);
        float deltaArea = glm::length(glm::cross(u2 - u0, u1 - u0)) / 2.0f;
        lArea += deltaArea;
      }
      radiance += throughput * mater.emission * mater.emission_strength *
                  lArea / (tt * tt);
      break;
    }
    // Constribution from the light source.
    glm::vec3 lpos, lnorm;
    float area;
    SampleFromLight(lpos, lnorm, area);
    if (area != 0.0f) {
      glm::vec3 ldir = glm::normalize(lpos - pos), lrad;
      HitRecord lhit;
      float lt = scene_->TraceRay(pos, ldir, 1e-4f, 1e7f, &lhit);
      if (lt > 0.0f || glm::length(pos + lt * ldir - lpos) < 0.1f) {
        Material lmater = scene_->GetEntity(lhit.hit_entity_id).GetMaterial();
        lrad = lmater.emission * lmater.emission_strength * brdf * area *
               glm::dot(norm, ldir) * glm::dot(-lnorm, ldir) / (lt * lt);
      }
      radiance += lrad * throughput;
    }
    // Constribution from other reflectors.
    float RRTest = random();
    if (RRTest >= RRProb && i > 1) {
      break;
    }
    dir = HemisphereSample(random(), random(), norm);
    if (glm::dot(dir, norm) < 0.0f) {
      dir = -dir;
    }
    throughput *= brdf * 2.0f * PI * glm::dot(norm, dir) / RRTest;
  }
  radiance = glm::clamp(radiance, glm::vec3{0.0f}, glm::vec3{1.0f});
  return radiance;
}

}  // namespace sparks
