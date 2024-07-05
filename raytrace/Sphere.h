#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "Util.h"
#include <cmath>

//using namespace std;
using Colour = cv::Vec3b;
typedef Eigen::Matrix<float, 3, 1> Vec3;

// Sphere Class
class Sphere
{
public:
    // Constructor that initializes a sphere with its center position and radius
    Sphere(Vec3 center, float radius) : mCenter_(center), mRadius_(radius) {}

    // Compute the intersection point of a ray and the sphere, returning the parameter 't'
    inline float intersectRay(Vec3 const& rayOrigin, Vec3 const& rayDirection) const {
        //Build the coefficients of the quadratic equation, A, B & C
        float A = rayDirection.dot(rayDirection);
        float B = 2.0f * rayDirection.dot(rayOrigin - mCenter_);
        float C = (rayOrigin - mCenter_).dot(rayOrigin - mCenter_) - mRadius_ * mRadius_;

        // Calculate the discriminant
        float discriminant = B * B - 4 * A * C;

        // If the discriminant is greater than zero, there is intersection
        if (discriminant > 0) {

            //find roots t1 and t2
            float t1 = (-B + sqrt(discriminant)) / (2 * A);
            float t2 = (-B - sqrt(discriminant)) / (2 * A);

            //determine which one is the real intersection point, the closest one min(t1, t2)
            float t = MIN(t1, t2);

            return t;
        }

        // returns NAN if no Intersection
        return NAN;
    }

    // Color the sphere at the intersection point using Blinn-Phong shading
    inline Colour colourSphere(Vec3 rayOrigin, Vec3 rayDirection, float t) { 

        // Calculate the intersection point using the ray parameter 't'
        Vec3 Intersection = rayDirection * t + rayOrigin;

        // Normal vector of the sphere at the intersection point
        Vec3 Normal = Intersection - mCenter_;
        Normal.normalize();

        //Compute difffuse term
        Vec3 LightVector = lightSource - Intersection;
        LightVector.normalize();
        float diffuseTerm = LightVector.dot(Normal);
        if (diffuseTerm < 0) diffuseTerm = 0;
        
        //Compute specular term
        Vec3 ViewVector = -rayDirection;
        Vec3 HalfVec = (LightVector + ViewVector).normalized();
        float specularTerm = HalfVec.dot(Normal); 
        if (specularTerm < 0) specularTerm = 0;

        // The reflection coefficients of sphere, k_a, k_d, k_s
        Vec3 k_a = Vec3(0.2f, 0.2f, 0.2f);
        Vec3 k_d = Vec3(0.3f, 0.8f, 1.0f);
        Vec3 k_s = Vec3(1.0f, 1.0f, 1.0f);

        // Specifying the ambient, diffuse and specular light colors 
        Vec3 ambLight(1.0f, 1.0f, 1.0f); // white
        Vec3 diffLight(1.0f, 1.0f, 1.0f); // white
        Vec3 specLight(1.0f, 1.0f, 1.0f); // white
        
        // Sphere Colour
        Vec3 objectColor; 
        Colour colour(0, 0, 0);

        // call the texture mapping function
        colour = shade_surface_point(Intersection);
        objectColor[0] = colour[0] / 255.0f;
        objectColor[1] = colour[1] / 255.0f;
        objectColor[2] = colour[2] / 255.0f;

        // Calculate the final color using the Phong shading model
        int shininess = 1000;
        colour[0] = Util::Clamp(((k_a[0] * ambLight[0] + k_d[0] * diffLight[0] * diffuseTerm + k_s[0] * specLight[0] * pow(specularTerm, shininess)) * objectColor[0]) * 255);
        colour[1] = Util::Clamp(((k_a[1] * ambLight[1] + k_d[1] * diffLight[1] * diffuseTerm + k_s[1] * specLight[1] * pow(specularTerm, shininess)) * objectColor[1]) * 255);
        colour[2] = Util::Clamp(((k_a[2] * ambLight[2] + k_d[2] * diffLight[2] * diffuseTerm + k_s[2] * specLight[2] * pow(specularTerm, shininess)) * objectColor[2]) * 255);

        return colour;
    }

    // Look up the texture from the world image
    Colour texture_lookup(float u, float v) {
        int i = round(v * world.cols - 0.5);
        int j = round(u * world.rows - 0.5);
        return world(i, j);
    }

    // get the Colour of texture at sphere point p 
    Colour shade_surface_point(Vec3 p) {
        Vec3 normal = p - mCenter_;
        float z = p[2];
        float y = p[1];
        float x = p[0];

        double u = (M_PI + atan2(y, x)) / (2 * M_PI);
        double v = (M_PI - acos(z / sqrt(x * x + y * y + z * z)));
        Colour diffuse_color = texture_lookup(u, v);
        return diffuse_color;
    }


private:
    // Coordinates of the sphere center (x, y, z)
    Vec3 mCenter_;
    // Radius of the sphere
    float mRadius_;
};