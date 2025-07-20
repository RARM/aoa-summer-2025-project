#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "proj_utils.h"

// === Util Function Prototypes ===

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

typedef struct {
    Point p1;
    Point p2;
    double dist;
} PointPair;

/**
 * @brief Find closest pair of points using the divide-and-conquer approach.
 *
 * The divide-and-conquer approach consists in dividing the points into two
 * halves, finding the closest pair in each half, and then checking if there is
 * a closer pair that crosses the division line.
 * 
 * This methods alters the points_arr array of points. It sorts the points by
 * their x values in place. The values return apply to the sorted points_arr
 * array.
 *
 * Time complexity: O(n log n)
 * Space complexity: O(n)
 *
 * @param points_arr Array of points.
 * @param points_num Number of points in the array.
 * @return A pointer to an array of integers representing the 0-indexed indices
 *         of the closest pair of points. E.g., { 1, 5 }.
 */
int* dnc_closest_pair(Point* points_arr, int points_num);

/**
 * Helper function to recursively find the closets points.
 * 
 * @param Px
 * @param Py
 * @param points_num
 * @returns
 */
PointPair dnc_closest_pair_rec(Point* Px, Point* Py, int points_num);

// Heap Sort

/**
 * Return the index of the parent of the given index.
 *
 * It performs the operation by shifting the number to the right by one bit.
 *
 * @param i Index of the element in the heap.
 * @return Index of the parent of the element in i.
 */
inline int heap_parent(int i) { return (i - 1) >> 1; }

/**
 * Get the index of the left child of a node in a heap given the element index.
 *
 * It performs the operation by shifting the index one bit to the left.
 *
 * @param i Index of the element in the heap.
 * @return Index of the left child in the heap.
 */
inline int heap_left(int i) {return (i << 1) + 1; }

/**
 * Get the index of the right child of a node in a heap given the index.
 *
 * The operation is performed shifting the index by the left one bit and adding
 * one.
 *
 * @param i Index of the element in the heap.
 * @return Index of the right child.
 */
inline int heap_right(int i) { return (i << 1) + 2; }

/**
 * Maintain the max-heap property after adding an item to the array.
 *
 * This function assume that the left and right subtrees are max-heaph (they
 * maintain the property of a max-heap). It performs the operation in-place.
 *
 * Note that the function orders the points by the x-coordinates.
 *
 * @param heap Pointer to the root of the max-heap.
 * @param s Size of the max-heap.
 * @param i Index of the element to heapify.
 * @param use_y Use the y value of each point to maintain the max heap.
 */
void max_heapify_points_by_x(Point* heap, int s, int i, int use_y);

/**
 * Given an array, this function generates a max heap in place.
 *
 * Note that this functions creates a max heap using the array of points and
 * their x value as reference.
 *
 * @param heap The contiguous array to convert to a max heap.
 * @param s Size of the input array.
 * @param use_y Use the y values of the points to build the max heap instead.
 */
void build_max_heap_points_by_x(Point* arr, int s, int use_y);

/**
 * Sort an array of points by their x coordinate using the heap sort algorithm.
 *
 * Time compexity: O(n log n)
 * Space complexity: O(1)
 *
 * @param arr Array of points to sort.
 * @param s Size of the array.
 * @param use_y Use the y value of the points for sorting instead.
 */
void heap_sort_points_by_x(Point* arr, int s, int use_y);

// === Implementation ===

int main(void) {
  FILE *out_fp = fopen("build/bin/performance.csv", "w");
  if (!out_fp) {
    fprintf(stderr, "[Main] Error: Could not open output file for writing.\n");
    return EXIT_FAILURE;
  }

  fprintf(out_fp, "Points,Brute Force Time (ms),Brute Force Time (ns),Divide and Conquer Time (ms),Divide and Conquer Time (ns)\n");

  // Get the points (generate or retrieve from file).
  for (int i = 0; i < 10; i++) {
    int points_num = 10000 * (i + 1);
    printf("=== [Main] Analyzing %d points. ===\n\n", points_num);

    double total_bf_time_ms = 0.0;
    double total_bf_time_ns = 0.0;
    double total_dnc_time_ms = 0.0;
    double total_dnc_time_ns = 0.0;

    for (int j = 0; j < 10; j++) { // 10 variations of points_num.
      char filename[50];
      snprintf(filename, sizeof(filename), "build/bin/data/points-%d-var-%d.txt", points_num, j + 1);
      Point* points_arr = get_points(points_num, filename);
      if (!points_arr) {
        fprintf(stderr, "[Main] Error: Failed to generate points.\n");
        return EXIT_FAILURE;
      }

      struct timespec start, end;
      double elapsed_time_ms, elapsed_time_ns;

      // I am using clock_gettime instead of time() because is available in 
      // non-POSIX systems and is more precise.

      clock_gettime(CLOCK_MONOTONIC, &start);
      int *bf_pair = bf_closest_pair(points_arr, points_num);
      clock_gettime(CLOCK_MONOTONIC, &end);

      Point* bf_p1 = points_arr + bf_pair[0];
      Point* bf_p2 = points_arr + bf_pair[1];

      // Calculate elapsed time in milliseconds and nanoseconds.
      elapsed_time_ms = (end.tv_sec - start.tv_sec) * 1000.0; // seconds to milliseconds
      elapsed_time_ms += (end.tv_nsec - start.tv_nsec) / 1000000.0; // nanoseconds to milliseconds

      elapsed_time_ns = (end.tv_sec - start.tv_sec) * 1000000000.0; // seconds to nanoseconds
      elapsed_time_ns += (end.tv_nsec - start.tv_nsec); // nanoseconds

      total_bf_time_ms += elapsed_time_ms;
      total_bf_time_ns += elapsed_time_ns;

      printf("Method: Brute Force        | Time: %10.2f ms | %20.2f ns | Points found: (%d, %d) and (%d, %d)\n",
             elapsed_time_ms, elapsed_time_ns, bf_p1->x, bf_p1->y, bf_p2->x, bf_p2->y);

      free(bf_pair); // Free the allocated memory for the pair.

      clock_gettime(CLOCK_MONOTONIC, &start);
      int *dnc_pair = dnc_closest_pair(points_arr, points_num);
      clock_gettime(CLOCK_MONOTONIC, &end);

      Point* dnc_p1 = points_arr + dnc_pair[0];
      Point* dnc_p2 = points_arr + dnc_pair[1];

      // Calculate elapsed time in milliseconds and nanoseconds.
      elapsed_time_ms = (end.tv_sec - start.tv_sec) * 1000.0; // seconds to milliseconds
      elapsed_time_ms += (end.tv_nsec - start.tv_nsec) / 1000000.0; // nanoseconds to milliseconds

      elapsed_time_ns = (end.tv_sec - start.tv_sec) * 1000000000.0; // seconds to nanoseconds
      elapsed_time_ns += (end.tv_nsec - start.tv_nsec); // nanoseconds

      printf("Method: Divide and Conquer | Time: %10.2f ms | %20.2f ns | Points found: (%d, %d) and (%d, %d)\n",
             elapsed_time_ms, elapsed_time_ns, dnc_p1->x, dnc_p1->y, dnc_p2->x, dnc_p2->y);

      total_dnc_time_ms += elapsed_time_ms;
      total_dnc_time_ns += elapsed_time_ns;

      free(dnc_pair); // Free the allocated memory for the pair.
      
      free_points(points_arr); // Free the allocated memory for points.
      // printf("[Main] Successfully generated %d points.\n", points_num);
    }

    fprintf(out_fp, "%d,%10.2f,%20.2f,%10.2f,%20.2f\n",
            points_num,
            total_bf_time_ms / 10.0,
            total_bf_time_ns / 10.0,
            total_dnc_time_ms / 10.0,
            total_dnc_time_ns / 10.0
    );

    printf("\n\n");
  }
  fclose(out_fp);
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

int* dnc_closest_pair(Point* points_arr, int points_num) {
  Point* Px = dup_points(points_arr, points_num);
  Point* Py = dup_points(points_arr, points_num);
  if (!Px || !Py) {
    fprintf(stderr, "[Main] Error: Failed to allocate memory for sorted point arrays.\n");
    free(Px);
    free(Py);
    return NULL;
  }

  heap_sort_points_by_x(Px, points_num, 0); // 0 for sorting by x
  heap_sort_points_by_x(Py, points_num, 1); // 1 for sorting by y

  PointPair closest_pair = dnc_closest_pair_rec(Px, Py, points_num);

  // The recursive function gives us the points. Now we must find their
  // indices in the Px array to meet the return requirement. This O(n) search
  // only happens once, so it doesn't affect the overall O(n log n) complexity.
  int* final_indices = malloc(sizeof(int) * 2);
  if (!final_indices) return NULL;
  
  final_indices[0] = -1;
  final_indices[1] = -1;

  for (int i = 0; i < points_num; i++) {
  if (points_arr[i].x == closest_pair.p1.x && points_arr[i].y == closest_pair.p1.y) {
    final_indices[0] = i;
  }
  if (points_arr[i].x == closest_pair.p2.x && points_arr[i].y == closest_pair.p2.y) {
    final_indices[1] = i;
  }
}
  
  // Handle case where the two points are the same in the final result
  if (final_indices[0] == final_indices[1]) {
    fprintf(stderr, "[Main] Error: POINTS ARE THE SAME.\n");
    // for (int i = 0; i < points_num; i++) {
    //     if (i != final_indices[0] && Px[i].x == closest_pair.p1.x && Px[i].y == closest_pair.p1.y) {
    //         final_indices[1] = i;
    //         break;
    //     }
    // }
  }

  free(Px);
  free(Py);

  return final_indices;
}

PointPair dnc_closest_pair_rec(Point* Px, Point* Py, int points_num) {
  // BASE CASE: If there are less than 2 points, return an invalid pair.
  if (points_num < 2) {
    PointPair result;
    result.dist = INFINITY; // Return infinite distance for invalid pairs
    return result;
  }

  // If there are 3 points, check the distances between them.
  if (points_num <= 3) {
    int *pairs = bf_closest_pair(Px, points_num);
    PointPair result;
    result.p1 = Px[pairs[0]];
    result.p2 = Px[pairs[1]];
    result.dist = points_dist(&result.p1, &result.p2);
    free(pairs);
    return result;
  }

  // DIVIDE: Find the midpoint and create Qx, Rx, Qy, Ry.
  int mid = points_num / 2;
  Point mid_point = Px[mid - 1];

  // Create Qx and Rx arrays
  Point* Qx = (Point*)malloc(sizeof(Point) * mid);
  Point* Rx = (Point*)malloc(sizeof(Point) * (points_num - mid));
  for(int i = 0; i < mid; i++) {
      Qx[i] = Px[i];
  }
  for(int i = 0; i < (points_num - mid); i++) {
      Rx[i] = Px[mid + i];
  }

  // Create Qy and Ry (the y-sorted arrays for the left and right halves).
  Point* Qy = (Point*)malloc(sizeof(Point) * mid);
  Point* Ry = (Point*)malloc(sizeof(Point) * (points_num - mid));
  int q_idx = 0;
  int r_idx = 0;
  for (int i = 0; i < points_num; i++) {
    if (Py[i].x < mid_point.x || (Py[i].x == mid_point.x && Py[i].y <= mid_point.y)) {
        if (q_idx < mid) {
            Qy[q_idx++] = Py[i];
        }
    } else {
      if (r_idx < (points_num - mid)) {
            Ry[r_idx++] = Py[i];
      }
    }
  }

  // CONQUER: Recursively find the closest pair in the left and right halves.
  PointPair pair_q = dnc_closest_pair_rec(Qx, Qy, mid);
  PointPair pair_r = dnc_closest_pair_rec(Rx, Ry, points_num - mid);

  // Free the memory for the subarrays
  free(Qx);
  free(Rx);
  free(Qy);
  free(Ry);

  // COMBINE
  // Determine the minimum distance (delta) from the two halves.
  PointPair min_pair = (pair_q.dist < pair_r.dist) ? pair_q : pair_r;
  double delta = min_pair.dist;

  // Constructing S. O(n) to create the strip.
  // Check for a closer pair across the dividing line (in the "strip").
  Point* strip = (Point*)malloc(sizeof(Point) * points_num);
  int strip_len = 0;
  for (int i = 0; i < points_num; i++) {
    if (abs(Py[i].x - mid_point.x) < delta) {
      strip[strip_len++] = Py[i];
    }
  }

  // Check for closer pairs within the strip.
  // For each point, we only need to check the next 15 points in the y-sorted strip.
  for (int i = 0; i < strip_len; i++) {
    for (int j = i + 1; j < strip_len && (j - i) <= 15; j++) {
      double d = points_dist(&strip[i], &strip[j]);
      if (d < min_pair.dist) {
        min_pair.dist = d;
        min_pair.p1 = strip[i];
        min_pair.p2 = strip[j];
      }
    }
  }

  // We don't need the final comparison in the original code because this code
  // is already including the pairs q and r in the comparision. The inclusion
  // happens when calculating delta.

  free(strip);
  return min_pair;
}

void max_heapify_points_by_x(Point* heap, int s, int i, int use_y) {
  int l = heap_left(i);
  int r = heap_right(i);
  int largest;

  if (
    l < s &&
    (use_y ? (heap + l)->y > (heap + i)->y : (heap + l)->x > (heap + i)->x)
  ) {
    largest = l;
  } else {
    largest = i;
  }

  if (
    r < s &&
    (use_y ? (heap + r)->y > (heap + largest)->y : (heap + r)->x > (heap + largest)->x)
  ) {
    largest = r;
  }

  if (largest != i) {
    Point temp = *(heap + i);
    *(heap + i) = *(heap + largest);
    *(heap + largest) = temp;

    max_heapify_points_by_x(heap, s, largest, use_y);
  }
}

void build_max_heap_points_by_x(Point* arr, int s, int use_y) {
  for (int i = (s >> 1) - 1; i >= 0; i--) {
    max_heapify_points_by_x(arr, s, i, use_y);
  }
}

void heap_sort_points_by_x(Point* arr, int s, int use_y) {
  build_max_heap_points_by_x(arr, s, use_y);
  int heap_s = s;

  for (int i = s - 1; i > 0; i--) {
    Point temp = *arr;
    *arr = *(arr + i);
    *(arr + i) = temp;

    heap_s--;
    max_heapify_points_by_x(arr, heap_s, 0, use_y);
  }
}