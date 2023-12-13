#pragma once
#include "glm/glm.hpp"
#include "sparks/assets/hit_record.h"
#include "sparks/assets/material.h"
#include "sparks/util/util.h"

namespace sparks{

    inline float cos(float theta) {
        return glm::cos(theta);
    }
    
    inline float cos2(float theta) {
        return glm::cos(theta) * glm::cos(theta);
    }

    inline float sin(float theta) {
        return glm::sin(theta);
    }

    inline float sin2(float theta) {
        return glm::sin(theta) * glm::sin(theta);
    }

    inline float sqrt(float r) {
        return glm::sqrt(r);
    }

    inline float mix(float u, float v, float r) {
        return u * (1.0f - r) + v * r;
    }

    inline glm::vec3 tint(glm::vec3 C) {
        float lum = 0.2126f * C.x + 0.7152 * C.y + 0.0722 * C.z;
        return C / lum;
    }
    glm::vec3 brdf(glm::vec3 &idir, glm::vec3 &odir, Material material);
}