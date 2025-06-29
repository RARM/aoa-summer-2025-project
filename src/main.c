#include <stdio.h>
#include <limits.h>
#include <math.h>

#include "proj_utils.h"

/**
 * @brief Calculates the Euclidean distance between two points.
 *
 * @param p1 Pointer to the first point.
 * @param p2 Pointer to the second point.
 * @return The Euclidean distance between p1 and p2.
 */
double points_dist(Point* p1, Point* p2);
 
/**
 * @brief Finds the closest pair of points using the brute-force approach.
 *
 * @param points_arr Array of points.
 * @param points_num Number of points in the array.
 * @return The closest pair of points as a Point structure.
 */
// Point bf_closest_pair(Point* points_arr, int points_num);

// === Implementation ===

int main(void) {
  printf("Generating 2 points.\n");
  Point* points_arr = generate_random_points(2, NULL);

  for (int i = 0; i < 2; i++) {
    printf(
      "Point %d: (%i, %i)\n", i,
      (points_arr + i)->x,
      (points_arr + i)->y
    );
  }

  printf("\nDistance: %.2f\n", points_dist(points_arr, (points_arr + 1)));

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