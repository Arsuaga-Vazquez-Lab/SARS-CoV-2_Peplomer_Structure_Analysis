/* Calculates ACN and space writhe. This should be much faster than doing it
   in python, since python is interpreted and C is compiled */
#define PI 3.14159265358979323846

#include <tuple>
#include <math.h>
#include <vector>
#include <iostream>

inline std::array<double, 3> unit(std::array<double, 3> vec) {
  // returns a vector (math vector, called array in C++) that points
  // in the same direction as the input, but with magnitude = 1
  double magnitude = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  if (magnitude == 0.) {
    throw std::runtime_error("Cannot find unit vector for the zero vector");
  }
  return { vec[0] / magnitude, vec[1] / magnitude, vec[2] / magnitude };
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
  // returns the dot product of inputs vec1 & vec2
  return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

std::array<double, 3> operator - (std::array<double, 3> vec1,
                                  std::array<double, 3> vec2) {
  // returns vec1 minus vec2
  return { vec1[0] - vec2[0], vec1[1] - vec2[1], vec1[2] - vec2[2] };
}

inline double sign(double num) {
  if (num > 0.) { return 1.; }
  else if (num < 0.) { return -1.; }
  else { return 0.; }
}

std::tuple<double, double> acn_and_writhe(std::vector<std::array<double, 3>> knot) {
  // given a knot as and array of points, returns a tuple (ACN, space writhe)
  double acn = 0;
  double space_writhe = 0;
  for (int i = 1; i < knot.size(); i++) {
    for (int j = 1; j < i - 1; i++) {
      std::array<double, 3> p1 = knot[i - 1];
      std::array<double, 3> p2 = knot[i];
      std::array<double, 3> p3 = knot[j - 1];
      std::array<double, 3> p4 = knot[j];
      // C definitely aint gonna let me just subtract vectors from each other that easily
      std::array<double, 3> n1 = unit(cross(p3 - p1, p4 - p1));
      std::array<double, 3> n2 = unit(cross(p4 - p1, p4 - p2));
      std::array<double, 3> n3 = unit(cross(p4 - p2, p3 - p2));
      std::array<double, 3> n4 = unit(cross(p3 - p2, p3 - p1));
      double omega_star = asin(dot(n1, n2)) +
                          asin(dot(n2, n3)) +
                          asin(dot(n3, n4)) +
                          asin(dot(n4, n1));
      omega_star /= 2 * PI;
      double omega_signed = omega_star * sign(dot(p3 - p1, cross(p4 - p3, p2 - p1)));
      acn += omega_star;
      space_writhe += omega_signed;
    }
  }
  return std::make_tuple(acn, space_writhe);
}
