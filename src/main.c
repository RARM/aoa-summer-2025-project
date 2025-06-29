#include <stdio.h>
#include <limits.h>

#include "proj_utils.h"

int main(void) {
  printf("Generating 5 points.\n");
  Point* points_arr = generate_random_points(5, NULL);

  for (int i = 0; i < 5; i++) {
    printf(
      "Point %d: (%i, %i)\n", i,
      (points_arr + i)->x,
      (points_arr + i)->y
    );
  }

  free_points(points_arr);

  return 0;
}