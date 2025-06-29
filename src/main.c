#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "proj_utils.h"

/**
 * @brief Calculates the Euclidean distance between two points.
 * 
 * Time complexity: O(1)
 * Space complexity: O(1)
 *
 * @param p1 Pointer to the first point.
 * @param p2 Pointer to the second point.
 * @return The Euclidean distance between p1 and p2.
 */
double points_dist(Point* p1, Point* p2);
 
/**
 * @brief Finds the closest pair of points using the brute-force approach.
 *
 * The brute-force approach consists in checking the distance of every point
 * with everything other point. It uses two nested loops.
 * 
 * Time complexity: O(n²)
 * Space complexity: O(1)
 * 
 * @param points_arr Array of points.
 * @param points_num Number of points in the array.
 * @return A pointer to an array of integers representing the 0-indexed indices
 *         of the closets pair of points. E.g., { 1, 5 }.
 */
int* bf_closest_pair(Point* points_arr, int points_num);

// === Implementation ===

int main(void) {
  Point* points_arr = generate_random_points(100000, NULL);
  printf("Generated 100000 points.\n");

  for (int i = 0; i < 100000; i++) {
    printf(
      "Point %i: (%i, %i)\n", i,
      (points_arr + i)->x,
      (points_arr + i)->y
    );
  }

  int* closets_pair = bf_closest_pair(points_arr, 100000);
  Point* p1 = points_arr + *closets_pair;
  Point* p2 = points_arr + *(closets_pair + 1);

  printf("\nThe closests points are:\n");
  printf("- Point %d: (%d, %d)\n", *closets_pair, p1->x, p1->y);
  printf("- Point %d: (%d, %d)\n", *(closets_pair + 1), p2->x, p2->y);
  printf("Distance: %f\n\n", points_dist(p1, p2));

  free(closets_pair);
  free_points(points_arr);

  return 0;
}

double points_dist(Point* p1, Point* p2) {
  double p1x = (double) p1->x;
  double p1y = (double) p1->y;
  double p2x = (double) p2->x;
  double p2y = (double) p2->y;

  return sqrt( // √[(x₂ - x₁)² + (y₂ - y₁)²]
    pow((double) (p2x - p1x), 2.0) +
    pow((double) (p2y - p1y), 2.0)
  );
}

int* bf_closest_pair(Point* points_arr, int points_num) {
  if (points_num < 2) return NULL;

  int* pairs = calloc(2, sizeof(int));
  pairs[0] = 0, pairs[1] = 1;
  
  // d_min is not inf, the min distance must be the dist of one pair.
  double d_min = points_dist(points_arr, points_arr + 1);
  double d_temp = d_min;

  for (int i = 0; i < points_num - 1; i++) {
    for (int j = i + 1; j < points_num; j++){
      d_temp = points_dist(points_arr + i, points_arr + j);
      if (d_temp < d_min) {
        d_min = d_temp;
        pairs[0] = i;
        pairs[1] = j;
      }
    }
  }

  return pairs;
}