#pragma once
#include "random"
#include "sparks/assets/scene.h"
#include "sparks/renderer/renderer_settings.h"

namespace sparks {
class PathTracer {
 public:
  PathTracer(const RendererSettings *render_settings, const Scene *scene);
  [[nodiscard]] glm::vec3 SampleRay(glm::vec3 origin,
                                    glm::vec3 direction,
                                    int x,
                                    int y,
                                    int sample) const;
  void SampleFromLight(glm::vec3 &lpos, glm::vec3 &lnorm, float &area) const; 
  glm::vec3 HemisphereSample(float u, float v, glm::vec3 norm) const;

 private:
  const RendererSettings *render_settings_{};
  const Scene *scene_{};
  static constexpr float RRProb = 0.95f;
  mutable std::mt19937 rd;
  
  float random() const {
    return std::uniform_real_distribution<float>(0.0f, 1.0f)(rd);
  }
};

}  // namespace sparks
