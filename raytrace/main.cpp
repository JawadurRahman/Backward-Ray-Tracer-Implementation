/* Assignment 1 of CSCI 4471/6671: Computer Graphics By Jawadur Rahman (2023). */

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "Image.h"
#include "Camera.h"
#include "ImagePlane.h"
#include "Sphere.h"
#include "Triangle.h"
#include <iostream>
using namespace std;

//Data types
typedef Eigen::Matrix<float, 3, 1> Vec3;
typedef Eigen::Matrix<float, 2, 1> Vec2;

//Color functions
using Colour = cv::Vec3b; // BGR Value: Type aliasing (Colour is just an alias for cv:Vec3b type)
Colour white() { return Colour(255, 255, 255); }
Colour black() { return Colour(0, 0, 0); }

int main(int, char**) {
	float epsilon = 0.00001;

	//Create an image object with 500 x 500 resolution.
	Image image = Image(500, 500);

	//Coordinates of image rectangle
	Vec3 llc = Vec3(-1, -1, -1);//lower left corner
	Vec3 urc = Vec3(1, 1, -1);//upper right corner
	ImagePlane imgPlane = ImagePlane(llc, urc, image.cols, image.rows); 

	// Specifying the camera position
	Vec3 cameraPoint = Vec3(0, 0, 0);
	Camera camera = Camera(cameraPoint);
	Vec3 Origin = camera.getPosition();

	// Defining the sphere: spherePos, sphereRadius 
	Vec3 SpherePos = Vec3(0, 0, -4); // Sphere Position
	float sphereRadius = 1.0f; // Sphere Radius
	Sphere sphere = Sphere(SpherePos, sphereRadius);

	// Floor Triangles
	Triangle triangle = Triangle(Vec3(-1, -1, -1), Vec3(10, -1, -10), Vec3(-10, -1, -10)); // floor triangle
	Triangle triangle2 = Triangle(Vec3(-1, -1, -1), Vec3(1, -1, -1), Vec3(10, -1, -10)); // floor triangle 2


    // Super-sampling parameters
    int h_samples = 5; // Number of horizontal samples per pixel
    int v_samples = 5; // Number of vertical samples per pixel
    int total_samples = h_samples * v_samples;

    for (int i = 0; i < image.rows; ++i) {
        for (int j = 0; j < image.cols; ++j) {
            Vec3 totalPixelColour(0, 0, 0);
            for (int dx = 0; dx < h_samples; ++dx) {
                for (int dy = 0; dy < v_samples; ++dy) {
                    // Build primary rays.
                    Vec3 pixelPos = imgPlane.generatePixelPos(i, j, dx, dy, h_samples, v_samples);
                    Vec3 rayDirec = camera.generateRay(pixelPos);
                    Vec3 unitRay = rayDirec.normalized();

                    // Purple to blue gradient Background
                    float a = 0.5 * (4 * unitRay.y());
                    Colour pixelColour = (a)*Colour(64, 6, 64) + (1.0 - a) * Colour(248, 206, 96);

                    Vec3 Intersection;
                    Colour floorColour;
                    // Render floor triangles if ray intersects with triangles
                    if (triangle.intersectRay(Origin, unitRay, epsilon, sphere, Intersection, floorColour)) {
                        pixelColour = floorColour;
                    }
                    if (triangle2.intersectRay(Origin, unitRay, epsilon, sphere, Intersection, floorColour)) {
                        pixelColour = floorColour;
                    }

                    // Check if ray intersects with sphere
                    float t = sphere.intersectRay(Origin, unitRay);
                    // if ray intersects with sphere then render it
                    if (t > 0) {
                        pixelColour = sphere.colourSphere(Origin, unitRay, t);
                    }
                    // The camera is inside the sphere
                    else if (!isnan(t)) {
                        pixelColour = black();
                    }

                    totalPixelColour += Vec3(pixelColour[0], pixelColour[1], pixelColour[2]);
                }
            }
            // Average the colors obtained from all samples
            totalPixelColour = totalPixelColour / total_samples;
            image(i, j) = Colour(totalPixelColour[0], totalPixelColour[1], totalPixelColour[2]);
        }
    }


	image.save("./result.png");
	image.display();

	return EXIT_SUCCESS;
}

