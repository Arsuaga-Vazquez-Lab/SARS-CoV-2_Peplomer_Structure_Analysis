/* Calculates ACN and space writhe of a piecewise-linear path.
   This is much faster than doing it in python, since python is interpreted and C is compiled
   This formula, as well as all the variable names, is taken from:
   https://en.wikipedia.org/wiki/Writhe#Numerically_approximating_the_Gauss_integral_for_writhe_of_a_curve_in_space

*/
#define PI 3.14159265358979323846

#include <tuple>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>

std::vector<std::array<double, 3>> points_csv(std::string filename) {
  // returns contents of file as a list of 3D points
  std::ifstream file (filename);
  std::vector<std::array<double, 3>> points;
  std::string line, coord;
  while (std::getline(file, line)) {
    std::array<double, 3> new_point;
    std::stringstream s (line);
    for (int i = 0; i <3 ; i++) {
      getline(s, coord, ',');
      new_point[i] = std::stod(coord);
    }
    points.push_back(new_point);
  }
  file.close();
  return points;
}

inline double magnitude(std::array<double, 3> vec) {
  return sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
}

inline std::array<double, 3> unit(std::array<double, 3> vec) {
  // returns a vector (math vector, called array in C++) that points
  // in the same direction as the input, but with magnitude = 1
  if (magnitude(vec) == 0.) {
    return {0., 0., 0.};
  }
  return { vec[0] / magnitude(vec), vec[1] / magnitude(vec), vec[2] / magnitude(vec) };
}

inline std::array<double, 3> cross(std::array<double, 3> vec1,
                                   std::array<double, 3> vec2) {
  // returns the cross product of inputs vec1 & vec2
  double x = vec1[1] * vec2[2] - vec1[2] * vec2[1];
  double y = vec1[2] * vec2[0] - vec1[0] * vec2[2];
  double z = vec1[0] * vec2[1] - vec1[1] * vec2[0];
  return { x, y, z };
}

inline double dot(std::array<double, 3> vec1,
                  std::array<double, 3> vec2) {
  // Returns the dot product of inputs vec1 & vec2
  return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

std::array<double, 3> operator - (std::array<double, 3> vec1,
                                  std::array<double, 3> vec2) {
  // Returns vec1 minus vec2
  return { vec1[0]-vec2[0], vec1[1]-vec2[1], vec1[2]-vec2[2] };
}

inline double sign(double num) { return num > 0. ? 1. : num < 0. ? -1. : 0.; }

inline double vec_sum(std::vector<double> vec) {
  double sum = 0;
  for (double num: vec) {
    sum += num;
  }
  return sum;
}

inline double expected_crossing_value(std::array<double, 3> p1,
                                      std::array<double, 3> p2,
                                      std::array<double, 3> p3,
                                      std::array<double, 3> p4) {
  // Returns the chance that line segment p1 -> p2 will cross p3 -> p4
  // when viewed from a random 2D projection. If the crossing is negative,
  // returns the same thing but negative
  std::array<double, 3> n1 = unit(cross(p3 - p1, p4 - p1));
  std::array<double, 3> n2 = unit(cross(p4 - p1, p4 - p2));
  std::array<double, 3> n3 = unit(cross(p4 - p2, p3 - p2));
  std::array<double, 3> n4 = unit(cross(p3 - p2, p3 - p1));
  for (std::array<double, 3> unit_vector : {n1, n2, n3, n4}) {
    // If any of the cross products resulted in magnitude zero,
    // the final result should be zero, and we need to catch that now
    if (magnitude(unit_vector) == 0.) {
      return 0.;
    }
  }
  double omega_star = asin(dot(n1, n2)) +
                      asin(dot(n2, n3)) +
                      asin(dot(n3, n4)) +
                      asin(dot(n4, n1));
  omega_star /= 4 * PI;
  double omega_signed = omega_star * sign(dot(p3 - p1, cross(p4 - p3, p2 - p1)));
  return omega_signed;
}

void add_crossing_values(int i,
                         std::vector<std::array<double, 3>>* knot_1,
                         std::vector<std::array<double, 3>>* knot_2,
                         std::vector<double>* acn_subcalcs,
                         std::vector<double>* space_writhe_subcalcs) {
  // Compares segment i of knot_1 to all segments of knot_2
  // This is to break up the calculations. This can be run in one thread,
  // but it would be slow to compare every segment of knot_1 to every segment
  // of knot_2 in just one thread. Since we don't want multiple threads touching
  // the same final calcs, this modifies one element of a vector passed by
  // pointer. That's the vector of all the subcalcs, which will be added together later.
  for (int j = 1; j < knot_2->size(); j++) {
    if (i == j) { continue; }
    std::array<double, 3> p1 = knot_1->at(i - 1);
    std::array<double, 3> p2 = knot_1->at(i);
    std::array<double, 3> p3 = knot_2->at(j - 1);
    std::array<double, 3> p4 = knot_2->at(j);
    double crossing_value = expected_crossing_value(p1, p2, p3, p4);
    space_writhe_subcalcs->at(i) += crossing_value;
    acn_subcalcs->at(i) += fabs(crossing_value);
  }
}

std::tuple<double, double> acn_and_writhe(std::vector<std::array<double, 3>> knot) {
  // Given a knot as a vector of points, returns the tuple (ACN, space writhe)
  // Since the calculation is a bunch of smaller things to be added together,
  // it splits that up into knot.size() different groups of calcs (subcalcs)
  // which each have their own thread. They modify values in an vector of subcalcs,
  // which is then added up to create the final calculations
  std::vector<double> acn_subcalcs (knot.size());
  std::vector<double> space_writhe_subcalcs (knot.size());
  std::vector<std::thread> threads = {};
  for (int i = 1; i < knot.size(); i++) {
    std::thread new_thread(add_crossing_values, i, &knot, &knot, &acn_subcalcs, &space_writhe_subcalcs);
    threads.push_back(std::move(new_thread));
  }
  for (std::thread& current_thread: threads) {
    current_thread.join();
  }
  double acn = vec_sum(acn_subcalcs);
  double space_writhe = vec_sum(space_writhe_subcalcs);
  return std::make_tuple(acn, space_writhe);
}

std::tuple<double, double> ACN_and_linking_number(
    std::vector<std::array<double, 3>> knot_A,
    std::vector<std::array<double, 3>> knot_B) {
  // This is very similar to the acn and space writhe calculations (see above)
  // except that it's iterating over two different knots. That is, it's
  // comparing each segment of one knot to a segment in the other, instead of
  // comparing every segment in a knot to every other segment in the same knot
  std::vector<double> acn_subcalcs (knot_A.size());
  std::vector<double> space_writhe_subcalcs (knot_A.size());
  std::vector<std::thread> threads = {};
  for (int i = 1; i < knot_A.size(); i++) {
    std::thread new_thread(add_crossing_values, i, &knot_A, &knot_B, &acn_subcalcs, &space_writhe_subcalcs);
    threads.push_back(std::move(new_thread));
  }
  for (std::thread& current_thread: threads) {
    current_thread.join();
  }
  // In this context, ACN is the number of times one knot crosses the other
  double acn = vec_sum(acn_subcalcs);
  double linking_number = vec_sum(space_writhe_subcalcs);
  return std::make_tuple(acn, linking_number);
}

int main() {
  std::vector<std::array<double, 3>> points_A = points_csv("points_A.csv");
  std::vector<std::array<double, 3>> points_B = points_csv("points_B.csv");
  double lnk_num = std::get<1>(ACN_and_linking_number(points_A, points_B));
  std::cout << "Linking number: " << lnk_num << std::endl;
  std::tuple<double, double> tuple_acn_and_space_writhe = acn_and_writhe(points_A);
  double acn          = std::get<0>(tuple_acn_and_space_writhe);
  double space_writhe = std::get<1>(tuple_acn_and_space_writhe);
  std::cout << "Hardware concurrency: " << std::thread::hardware_concurrency() << std::endl;
  std::cout << "ACN: " << acn << "\nSpace Writhe: " << space_writhe << std::endl;
  return 0;  // executed successfully
}
