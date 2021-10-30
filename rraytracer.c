#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "rtga/rtga.h"

#define HEIGHT 400
#define WIDTH 600
#define ASPECT_RATIO ((float)WIDTH / HEIGHT)

#define VIEW_HEIGHT 2.0f
#define VIEW_WIDTH (ASPECT_RATIO * VIEW_HEIGHT)
#define FOCAL_LENGTH 1.0f

// This program is based off of Ray Tracing in One Weekend by Peter Shirley (https://raytracing.github.io/books/RayTracingInOneWeekend.html)

struct vec3 {
    float x, y, z;
};

struct vec3 vec_add(struct vec3 v1, struct vec3 v2) {
    return (struct vec3){v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

struct vec3 vec_sub(struct vec3 v1, struct vec3 v2) {
    return (struct vec3){v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}

struct vec3 vec_multf(struct vec3 v, float f) {
    return (struct vec3){v.x * f, v.y * f, v.z * f};
}

struct vec3 vec_div(struct vec3 v1, struct vec3 v2) {
    return (struct vec3){v1.x / v2.x, v1.y / v2.y, v1.z / v2.z};
}

struct vec3 vec_divf(struct vec3 v, float f) {
    return (struct vec3){v.x / f, v.y / f, v.z / f};
}

float vec_length(struct vec3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

struct vec3 vec_unit(struct vec3 v) {
    return vec_divf(v, vec_length(v));
}

struct ray {
    struct vec3 origin, direction;
};

TgaColor ray_color(struct ray r) {
    struct vec3 unit_direction = vec_unit(r.direction);
    float t = 0.5f * (unit_direction.y + 1.0f);
    struct vec3 c1 = {1.0f, 1.0f, 1.0f};
    struct vec3 c2 = {0.5f, 0.7f, 1.0f};
    c1 = vec_multf(c1, 1.0f - t);
    c2 = vec_multf(c2, t);
    c1 = vec_add(c1, c2);
    return (TgaColor)COLOR24(c1.x * 255, c1.y * 255, c1.z * 255);
}

int main(void) {
    TgaImage image;
    int success;

    struct vec3 origin = {0, 0, 0};
    struct vec3 horizontal = {VIEW_WIDTH, 0, 0};
    struct vec3 vertical = {0, VIEW_HEIGHT, 0};
    struct vec3 focal_length = {0, 0, FOCAL_LENGTH};
    struct vec3 low_left_corner = vec_sub(vec_sub(vec_sub(
        origin, vec_divf(horizontal, 2.0f)),
        vec_divf(vertical, 2.0f)),
        focal_length);

    // Set image struct to all 0s
    memset(&image, 0, sizeof(image));

    // Allocate image
    success = tga_alloc(UNCOMPRESSED_TRUE_COLOR_IMAGE, WIDTH, HEIGHT, 24, &image);
    if (success != TGA_SUCCESS) {
        printf("Error %i occured while allocating image.\n", success);
        return 1;
    }

    // Trace rays for each pixel in image
    for (int y = 0; y < HEIGHT; ++y) {
        printf("\r\e[2K");
        printf("\rScanlines remaining: %i", HEIGHT - y);
        fflush(stdout);
        for (int x = 0; x < WIDTH; ++x) {
            float u = (float)x / (WIDTH - 1.0f);
            float v = (float)y / (HEIGHT - 1.0f);
            struct ray r = {origin, vec_sub(vec_add(vec_add(
                low_left_corner, vec_multf(horizontal, u)),
                vec_multf(vertical, v)), origin)};
            TgaColor color = ray_color(r);
            tga_set_pixel(&image, x, y, color);
        }
    }
    printf("\r\e[2K");
    printf("\rRaytracer finished!\n");


    // Write image to file
    success = tga_write_file(&image, "raytracer.tga");
    if (success != TGA_SUCCESS) {
        printf("Error occured while writing image to a file.\n");
        return 2;
    }

    // Free image
    tga_free(&image);

    return 0;
}
