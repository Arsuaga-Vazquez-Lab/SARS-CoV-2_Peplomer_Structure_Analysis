/* Calculates ACN and space writhe. This should be much faster than doing it
   in python, since python is interpreted and C is compiled */
#include <tuple>
#include <math>

inline std::vector<3, double> unit(std::vector<3, double> vec) {
  // returns a vector in the same direction as he input, but with magnitude one
  double magnitude = math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  if (magnitude == 0.) {
    throw Exception;
  }
  return [vec[0] / magnitude, vec[1] / magnitude, vec[2] / magnitude];
}

inline std::vector<3, double> cross(std::vector<3, double> vec1, std::vector<3, double> vec2) {
  // returns the cross product of inputs vec1 & vec2
  double x = vec1[1] * vec2[2] - vec1[2] * vec2[1];
  double y = vec1[2] * vec2[0] - vec1[0] * vec2[2];
  double z = vec1[0] * vec2[1] - vec1[1] * vec2[0];
  return [x, y, z];
}

inline double dot(std::vector<3, double> vec1, std::vector<3, double> vec2) {
  // returns the dot product of inputs vec1 & vec2
  return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

inline double sign(double num) {
  if (num > 0.) { return 1.; }
  else if (num < 0.) { return -1.; }
  else { return 0.; }
}

tuple<double, double> acn_and_writhe(std::array<std::vector<3, double>> knot) {
  // given a knot as and array of points, returns a tuple (ACN, space writhe)
  double acn = 0;
  double space_writhe = 0;
  for (int i = 1; i < knot.length(); i++) {
    for (int j = 1; j < i - 1; i++) {
      std::vector<3, double> p1 = knot[i - 1];
      std::vector<3, double> p2 = knot[i];
      std::vector<3, double> p3 = knot[j - 1];
      std::vector<3, double> p4 = knot[j];
      // C definitely aint gonna let me just subtract vectors from each other that easily
      std::vector<3, double> n1 = unit(cross(p3 - p1, p4 - p1));
      std::vector<3, double> n2 = unit(cross(p4 - p1, p4 - p2));
      std::vector<3, double> n3 = unit(cross(p4 - p2, p3 - p2));
      std::vector<3, double> n4 = unit(cross(p3 - p2, p3 - p1));
      double omega_star = math.asin(dot(n1, n2) +
                          math.asin(dot(n2, n3) +
                          math.asin(dot(n3, n4) +
                          math.asin(dot(n4, n1);
      omega_star /= 2 * math.pi;
      double omega_signed = omega_star * sign(dot(p3 - p1, cross(p4 - p3, p2 - p1)));
      acn += omega_star;
      space_writhe += omega_signed;
  return (acn, space_writhe);
}
