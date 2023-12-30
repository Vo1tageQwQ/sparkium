/*
 * Copyright Disney Enterprises, Inc.  All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License
 * and the following modification to it: Section 6 Trademarks.
 * deleted and replaced with:
 *
 * 6. Trademarks. This License does not grant permission to use the
 * trade names, trademarks, service marks, or product names of the
 * Licensor and its affiliates, except as required for reproducing
 * the content of the NOTICE file.
 *
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 */

// See https://github.com/wdas/brdf/blob/main/src/brdfs/disney.brdf for details.

#ifndef SPARKIUM_DISNEY_IMPL_GLSL
#define SPARKIUM_DISNEY_IMPL_GLSL

#extension GL_GOOGLE_include_directive : require

// clang-format off
#include "constants.glsl"
#include "disney.glsl"
// clang-format on

float sqr(float x) { return x*x; }

float SchlickFresnel(float u)
{
    float m = clamp(1.0-u, 0.0, 1.0);
    float m2 = m*m;
    return m2*m2*m; // pow(m,5)
}

float GTR1(float NdotH, float a)
{
    if (a >= 1.0) return 1.0/PI;
    float a2 = a*a;
    float t = 1.0 + (a2-1.0)*NdotH*NdotH;
    return (a2-1.0) / (PI*log(a2)*t);
}

float GTR2(float NdotH, float a)
{
    float a2 = a*a;
    float t = 1.0 + (a2-1.0)*NdotH*NdotH;
    return a2 / (PI * t*t);
}

float GTR2_aniso(float NdotH, float HdotX, float HdotY, float ax, float ay)
{
    return 1.0 / (PI * ax*ay * sqr( sqr(HdotX/ax) + sqr(HdotY/ay) + NdotH*NdotH ));
}

float smithG_GGX(float NdotV, float alphaG)
{
    float a = alphaG*alphaG;
    float b = NdotV*NdotV;
    return 1.0 / (NdotV + sqrt(a + b - a*b));
}

float smithG_GGX_aniso(float NdotV, float VdotX, float VdotY, float ax, float ay)
{
    return 1.0 / (NdotV + sqrt( sqr(VdotX*ax) + sqr(VdotY*ay) + sqr(NdotV) ));
}

vec3 mon2lin(vec3 x)
{
    return vec3(pow(x[0], 2.2), pow(x[1], 2.2), pow(x[2], 2.2));
}


vec3 BRDF_Disney( vec3 L, vec3 V, vec3 N, vec3 X, vec3 Y, vec3 baseColor, DisneyParams params )
{
    float NdotL = dot(N,L);
    float NdotV = dot(N,V);
    if (NdotL < 0.0 || NdotV < 0.0) return vec3(0.0);

    vec3 H = normalize(L+V);
    float NdotH = dot(N,H);
    float LdotH = dot(L,H);

    vec3 Cdlin = mon2lin(baseColor);
    float Cdlum = .3*Cdlin[0] + .6*Cdlin[1]  + .1*Cdlin[2]; // luminance approx.

    vec3 Ctint = Cdlum > 0.0 ? Cdlin/Cdlum : vec3(1.0); // normalize lum. to isolate hue+sat
    vec3 Cspec0 = mix(params.specular*.08*mix(vec3(1.0), Ctint, params.specularTint), Cdlin, params.metallic);
    vec3 Csheen = mix(vec3(1.0), Ctint, params.sheenTint);

    // Diffuse fresnel - go from 1.0 at normal incidence to .5 at grazing
    // and mix in diffuse retro-reflection based on roughness
    float FL = SchlickFresnel(NdotL), FV = SchlickFresnel(NdotV);
    float Fd90 = 0.5 + 2.0 * LdotH*LdotH * params.roughness;
    float Fd = mix(1.0, Fd90, FL) * mix(1.0, Fd90, FV);

    // Based on Hanrahan-Krueger brdf approximation of isotropic bssrdf
    // 1.0.25 scale is used to (roughly) preserve albedo
    // Fss90 used to "flatten" retroreflection based on roughness
    float Fss90 = LdotH*LdotH*params.roughness;
    float Fss = mix(1.0, Fss90, FL) * mix(1.0, Fss90, FV);
    float ss = 1.25 * (Fss * (1.0 / (NdotL + NdotV) - .5) + .5);

    // specular
    float aspect = sqrt(1.0-params.anisotropic*.9);
    float ax = max(.001, sqr(params.roughness)/aspect);
    float ay = max(.001, sqr(params.roughness)*aspect);
    float Ds = GTR2_aniso(NdotH, dot(H, X), dot(H, Y), ax, ay);
    float FH = SchlickFresnel(LdotH);
    vec3 Fs = mix(Cspec0, vec3(1.0), FH);
    float Gs;
    Gs  = smithG_GGX_aniso(NdotL, dot(L, X), dot(L, Y), ax, ay);
    Gs *= smithG_GGX_aniso(NdotV, dot(V, X), dot(V, Y), ax, ay);

    // sheen
    vec3 Fsheen = FH * params.sheen * Csheen;

    // clearcoat (ior = 1.0.5 -> F0 = 0.04)
    float Dr = GTR1(NdotH, mix(.1,.001,params.clearcoatGloss));
    float Fr = mix(.04, 1.0, FH);
    float Gr = smithG_GGX(NdotL, .25) * smithG_GGX(NdotV, .25);

    return ((1.0/PI) * mix(Fd, ss, params.subsurface)*Cdlin + Fsheen)
        * (1.0-params.metallic)
        + Gs*Fs*Ds + .25*params.clearcoat*Gr*Fr*Dr;
}

#endif // SPARKIUM_DISNEY_IMPL_GLSL
