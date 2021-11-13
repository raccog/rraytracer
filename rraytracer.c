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
#define SAMPLES 500
#define MAX_DEPTH 50

#define VIEW_HEIGHT 2.0f
#define VIEW_WIDTH (ASPECT_RATIO * VIEW_HEIGHT)
#define FOCAL_LENGTH 1.0f

#define GROUND {{0.8f, 0.8f, 0.0f}, METAL, 1, 0}
#define CENTER {{0.3f, 0.3f, 0.7f}, DIFFUSE, 0, 0}
#define LEFT {{0.8f, 0.8f, 0.8f}, DIELECTRIC, 0, 1.5f}
#define RIGHT {{0.8f, 0.2f, 0.2f}, METAL, 0.1f, 0}

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

float reflectance(float cosine, float ref_idx) {
    float r0 = (1 - ref_idx) / (1 + ref_idx);
    r0 *= r0;
    return r0 + (1 - r0) * pow(1 - cosine, 5);
}

struct vec3 {
    float x, y, z;
};

enum material_type {
    DIFFUSE,
    METAL,
    DIELECTRIC
};

struct material {
    struct vec3 color;
    enum material_type type;
    float fuzz;
    float refraction;
};

struct sphere {
    struct vec3 center;
    float radius;
    struct material mat;
};

struct sphere spheres[] = {
    {{0, 0, -1}, 0.5f, CENTER},
    {{0, -100.5f, -1}, 100, GROUND},
    {{-1, 0, -1}, -0.5f, LEFT},
    {{1, 0.5f, -1}, 0.5f, RIGHT}
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

struct vec3 vec_subf(struct vec3 v, float f) {
    return (struct vec3){v.x - f, v.y - f, v.z - f};
}

struct vec3 vec_mult(struct vec3 v1, struct vec3 v2) {
    return (struct vec3){v1.x * v2.x, v1.y * v2.y, v1.z * v2.z};
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

bool vec_near_zero(struct vec3 v) {
    float s = 1e-8;
    return (fabs(v.x) < s) && (fabs(v.y) < s) && (fabs(v.z) < s);
}

struct vec3 vec_reflect(struct vec3 v, struct vec3 n) {
    return vec_sub(v, vec_multf(n, 2 * vec_dot(v, n)));
}

struct vec3 vec_refract(struct vec3 v, struct vec3 n, float etai_over_etat) {
    float cos_theta = fmin(vec_dot(vec_sub((struct vec3){0,0,0}, v), n), 1.0f);
    struct vec3 r_out_perp = vec_multf(vec_add(v, vec_multf(n, cos_theta)), etai_over_etat);
    struct vec3 r_out_parallel = vec_multf(n, -sqrt(fabs(1.0f - vec_length_squ(r_out_perp))));
    return vec_add(r_out_perp, r_out_parallel);
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
    struct material mat;
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
    rec->mat = s.mat;

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

struct vec3 randf_unit_vec() {
    return vec_unit(randf_in_unit_sphere());
}

struct vec3 randf_in_hemisphere(struct vec3 normal) {
    struct vec3 in_unit_sphere = randf_in_unit_sphere();
    if (vec_dot(in_unit_sphere, normal) > 0.0f) {
        return in_unit_sphere;
    } else {
        return vec_sub((struct vec3){0, 0, 0}, in_unit_sphere);
    }
}

bool scatter_diffuse(struct ray input, struct hit_record rec, struct vec3 *attenuation, struct ray *scattered) {
    struct vec3 scatter_direction = vec_add(rec.normal , randf_unit_vec());
    if (vec_near_zero(scatter_direction)) {
        scatter_direction = rec.normal;
    }
    scattered->origin = rec.p;
    scattered->direction = scatter_direction;
    memcpy(attenuation, &rec.mat.color, sizeof(struct vec3));

    return true;
}

bool scatter_metal(struct ray input, struct hit_record rec, struct vec3 *attenuation, struct ray *scattered) {
    struct vec3 reflected = vec_reflect(vec_unit(input.direction), rec.normal);
    scattered->origin = rec.p;
    scattered->direction = vec_add(reflected, vec_multf(randf_in_unit_sphere(), rec.mat.fuzz));
    memcpy(attenuation, &rec.mat.color, sizeof(struct vec3));
    return (vec_dot(scattered->direction, rec.normal) > 0);
}

bool scatter_refract(struct ray input, struct hit_record rec, struct vec3 *attenuation, struct ray *scattered) {
    attenuation->x = 1;
    attenuation->y = 1;
    attenuation->z = 1;
    float refraction_ratio = rec.front_face ? (1.0f / rec.mat.refraction) : rec.mat.refraction;

    struct vec3 unit_direction = vec_unit(input.direction);
    float cos_theta = fmin(vec_dot(vec_sub((struct vec3){0,0,0}, unit_direction), rec.normal), 1.0f);
    float sin_theta = sqrt(1.0f - cos_theta * cos_theta);

    bool cannot_refract = refraction_ratio * sin_theta > 1.0f;
    struct vec3 direction;

    if (cannot_refract || reflectance(cos_theta, refraction_ratio) > randf()) {
        direction = vec_reflect(unit_direction, rec.normal);
    } else {
        direction = vec_refract(unit_direction, rec.normal, refraction_ratio);
    }

    scattered->origin = rec.p;
    scattered->direction = direction;

    return true;
}

struct vec3 ray_color(struct ray r, int depth) {
    struct hit_record rec;
    if (depth <= 0) return (struct vec3){0, 0, 0};
    if (world_hit(r, 0.001f, INFINITY, &rec)) {
        struct ray scattered;
        struct vec3 attenuation;
        bool success;
        if (rec.mat.type == DIFFUSE) {
            success = scatter_diffuse(r, rec, &attenuation, &scattered);
        } else if (rec.mat.type == METAL) {
            success = scatter_metal(r, rec, &attenuation, &scattered);
        } else if (rec.mat.type == DIELECTRIC) {
            success = scatter_refract(r, rec, &attenuation, &scattered);
        }
        if (success) {
            return vec_mult(attenuation, ray_color(scattered, depth - 1));
        }
        return (struct vec3){0, 0, 0};
    }
    struct vec3 unit_direction = vec_unit(r.direction);
    float t = 0.5f * (unit_direction.y + 1.0f);
    struct vec3 c1 = {1.0f, 1.0f, 1.0f};
    struct vec3 c2 = {0.5f, 0.7f, 1.0f};
    c1 = vec_multf(c1, 1.0f - t);
    c2 = vec_multf(c2, t);
    c1 = vec_add(c1, c2);
    return c1;
}

void write_color(TgaImage *image, int x, int y, struct vec3 color, int samples) {
    float scale = 1.0f / samples;
    color = vec_multf(color, scale);
    color.x = sqrt(color.x);
    color.y = sqrt(color.y);
    color.z = sqrt(color.z);
    color = vec_clamp(color, 0.0f, 1.0f);
    color = vec_multf(color, 255);

    tga_set_pixel(image, x, y, COLOR24(color.x, color.y, color.z));
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
