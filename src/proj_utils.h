/**
 * @file proj_utils.h
 * @brief Utility functions to generate and manage arrays of random 2D points.
 *
 * This header provides functions to create an array of random 2D points
 * (x, y), and optionally save them to a text file.
 */

#ifndef PROJ_UTILS_H
#define PROJ_UTILS_H

#include <stddef.h>

/**
 * @struct Point
 * @brief Represents a 2D point with integer coordinates.
 */
typedef struct {
  int x; /**< X coordinate */
  int y; /**< Y coordinate */
} Point;

/**
 * Retrieves or generates an array of Point structures.
 *
 * If `filename` is NULL, the function generates a set of points without saving
 * them. If `filename` is not NULL, the function attempts to retrieve points
 * from the specified file. If the file does not exist, the function generates
 * the points and saves them to the given filename.
 *
 * @param num_points The number of points to generate or retrieve.
 * @param filename Path to the file to retrieve or save points. NULL to
 *                 generate points without saving.
 * @return Pointer to an array of Point structures. The caller is responsible
 *         for freeing the memory.
 */

Point* get_points(size_t num_points, const char* filename);

/**
 * Create (allocate) an array of points copying from `points_arr`.
 * 
 * @param points_arr Pointer to the array of pointers.
 * @param num_points Array size.
 * @return A new array with the duplicated array of points. The caller is
 *         responsible for freeing the memory.
 */
Point* dup_points(Point* points_arr, size_t num_points);

/**
 * @brief Generates an array of unique random 2D points.
 *
 * Allocates an array of `num_points` random points on the heap.
 * Optionally saves the generated points to a text file if `filename` is not
 * NULL.
 *
 * @param num_points The number of points to generate.
 * @param filename The name of the file to save the points to (one point per
 *                 line as "x y"), or NULL to skip saving.
 * @return Pointer to the dynamically allocated array of points, or NULL on
 *         failure.
 * @note The caller is responsible for freeing the returned array using
 *       free_points().
 */
Point* generate_random_points(size_t num_points, const char* filename);

/**
 * @brief Frees an array of points allocated by generate_random_points().
 *
 * @param points Pointer to the array of points to free.
 */
void free_points(Point* points);

#endif // PROJ_UTILS_H