#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "rtga/rtga.h"

#define HEIGHT 400
#define WIDTH 600
#define ASPECT_RATIO ((float)WIDTH / HEIGHT)
#define SAMPLES 100
#define MAX_DEPTH 50

#define VIEW_HEIGHT 2.0f
#define VIEW_WIDTH (ASPECT_RATIO * VIEW_HEIGHT)
#define FOCAL_LENGTH 1.0f

// This program is based off of Ray Tracing in One Weekend by Peter Shirley (https://raytracing.github.io/books/RayTracingInOneWeekend.html)

float randf() {
    return rand() / (RAND_MAX + 1.0f);
}

float randf_range(float min, float max) {
    return min + (max - min) * randf();
}

float clamp(float x, float min, float max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

struct vec3 {
    float x, y, z;
};

struct sphere {
    struct vec3 center;
    float radius;
};

struct sphere spheres[] = {
    {{0, 0, -1}, 0.5f},
    {{0, -100.5f, -1}, 100},
};

struct vec3 vec_add(struct vec3 v1, struct vec3 v2) {
    return (struct vec3){v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}

struct vec3 vec_addf(struct vec3 v, float f) {
    return (struct vec3){v.x + f, v.y + f, v.z + f};
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

float vec_length_squ(struct vec3 v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

float vec_length(struct vec3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

struct vec3 vec_unit(struct vec3 v) {
    return vec_divf(v, vec_length(v));
}

float vec_dot(struct vec3 v1, struct vec3 v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

struct vec3 vec_clamp(struct vec3 v, float min, float max) {
    return (struct vec3){
        clamp(v.x, min, max),
        clamp(v.y, min, max),
        clamp(v.z, min, max)
    };
}

struct vec3 vec_randf() {
    return (struct vec3){randf(), randf(), randf()};
}

struct vec3 vec_randf_range(float min, float max) {
    return (struct vec3){
        randf_range(min, max),
        randf_range(min, max),
        randf_range(min, max)
    };
}

struct ray {
    struct vec3 origin, direction;
};

struct vec3 ray_at(struct ray r, float t) {
    return vec_add(r.origin, vec_multf(r.direction, t));
}

struct hit_record {
    struct vec3 p, normal;
    float t;
    bool front_face;
};

void record_face_normal(struct hit_record *rec, struct ray r, struct vec3 outward_normal) {
    rec->front_face = (vec_dot(r.direction, outward_normal) < 0);
    rec->normal = rec->front_face 
        ? outward_normal 
        : vec_sub((struct vec3){0, 0, 0}, outward_normal);
}

bool sphere_hit(struct sphere s, struct ray r, float t_min, float t_max, struct hit_record *rec) {
    struct vec3 oc = vec_sub(r.origin, s.center);
    float a = vec_length_squ(r.direction);
    float half_b = vec_dot(oc, r.direction);
    float c = vec_length_squ(oc) - s.radius * s.radius;
    float discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return false;

    float sqrtd = sqrt(discriminant);

    float root = (-half_b - sqrtd) / a;
    if (root < t_min || root > t_max) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || root > t_max) return false;
    }

    rec->t = root;
    rec->p = ray_at(r, rec->t);
    struct vec3 outward_normal = vec_divf(vec_sub(rec->p, s.center), s.radius);
    record_face_normal(rec, r, outward_normal);

    return true;
}

bool world_hit(struct ray r, float t_min, float t_max, struct hit_record *rec) {
    struct hit_record tmp_rec;
    float closest = t_max;
    bool has_hit = false;

    for (int i = 0; i < sizeof(spheres) / sizeof(struct sphere); ++i) {
        if (sphere_hit(spheres[i], r, t_min, closest, &tmp_rec)) {
            has_hit = true;
            closest = tmp_rec.t;
            memcpy(rec, &tmp_rec, sizeof(struct hit_record));
        }
    }

    return has_hit;
}

struct vec3 randf_in_unit_sphere() {
    while (true) {
        struct vec3 p = vec_randf_range(-1.0f, 1.0f);
        if (vec_length_squ(p) >= 1.0f) continue;
        return p;
    }
}

struct vec3 ray_color(struct ray r, int depth) {
    struct hit_record rec;
    if (depth <= 0) return (struct vec3){0, 0, 0};
    if (world_hit(r, 0, INFINITY, &rec)) {
        struct vec3 target = vec_add(vec_add(rec.p, rec.normal), randf_in_unit_sphere());
        return vec_multf(ray_color((struct ray){rec.p, vec_sub(target, rec.p)}, depth - 1), 0.5f);
    }
    struct vec3 unit_direction = vec_unit(r.direction);
    float t = 0.5f * (unit_direction.y + 1.0f);
    struct vec3 c1 = {1.0f, 1.0f, 1.0f};
    struct vec3 c2 = {0.5f, 0.7f, 1.0f};
    c1 = vec_multf(c1, 1.0f - t);
    c2 = vec_multf(c2, t);
    c1 = vec_add(c1, c2);
    return vec_multf(c1, 255);
}

void write_color(TgaImage *image, int x, int y, struct vec3 color, int samples) {
    float scale = 1.0f / samples;
    color = vec_multf(color, scale);
    color = vec_clamp(color, 0.0f, 255.0f);

    tga_set_pixel(image, x, y, (TgaColor)COLOR24(color.x, color.y, color.z));
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
            struct vec3 color = {0, 0, 0};
            for (int s = 0; s < SAMPLES; ++s) {
                float u = (x + randf()) / (WIDTH - 1.0f);
                float v = (y + randf()) / (HEIGHT - 1.0f);
                struct ray r = {origin, vec_sub(vec_add(vec_add(
                    low_left_corner, vec_multf(horizontal, u)),
                    vec_multf(vertical, v)), origin)};
                color = vec_add(color, ray_color(r, MAX_DEPTH));
            }
            write_color(&image, x, y, color, SAMPLES);
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
