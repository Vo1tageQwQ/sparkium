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
                                    int sample);
  void SampleFromLight(glm::vec3 &lpos, glm::vec3 &lnorm, float &area); 
  glm::vec3 HemisphereSample(float u, float v, glm::vec3 norm);

 private:
  const RendererSettings *render_settings_{};
  const Scene *scene_{};
  const float RRProb = 0.95f;
  const float PI = 3.1415926535f;
  std::uniform_real_distribution<float> uniform;
  std::mt19937 rd;
};
}  // namespace sparks
