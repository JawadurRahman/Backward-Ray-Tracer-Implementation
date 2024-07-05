#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>

typedef Eigen::Matrix<float, 3, 1> Vec3;

// ImagePlane Class
class ImagePlane
{
public:
    // Constructor that initializes the image plane with its corners and dimensions
    ImagePlane(Vec3 lowerLeft, Vec3 upperRight, int xCols, int yRows)
        : lowerLeft_(lowerLeft), upperRight_(upperRight), xCols_(xCols), yRows_(yRows) {};

    // Method to generate the position of a pixel given its row (i) and column (j) indices
    Vec3 generatePixelPos(int i, int j) {
        int width = upperRight_(0) - lowerLeft_(0); // width of the image plane
        int height = upperRight_(1) - lowerLeft_(1); // height of the image plane
        float u = lowerLeft_(0) + (width) * (j + 0.5f) / xCols_;
        float v = lowerLeft_(1) + 2 - (height) * (i + 0.5f) / yRows_;
        return Vec3(u, v, -1);
    }

    Vec3 generatePixelPos(int i, int j, int dx, int dy, int h_samples, int v_samples) {
        float randomX = ((double)rand()) / RAND_MAX;
        float randomY = ((double)rand()) / RAND_MAX;
        int width = upperRight_(0) - lowerLeft_(0); // width of the image plane
        int height = upperRight_(1) - lowerLeft_(1); // height of the image plane
        float u = lowerLeft_(0) + (width) * (j + (randomX + dx)/h_samples) / xCols_;
        float v = lowerLeft_(1) + 2 - (height) * (i + (randomY + dy)/v_samples) / yRows_;
        return Vec3(u, v, -1);
    }

private:
    //lower-left corner of the image plane
    Vec3 lowerLeft_;
    //upper-right corner of the image plane
    Vec3 upperRight_;
    // number of columns
    int xCols_;
    // number of rows
    int yRows_;
};

