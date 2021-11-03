#include <Rcpp.h>
using namespace Rcpp;

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

const double SQRT_2 = 1.4142135623730950;
const double TWO_PI = 6.283185307179586;

class Point {
public:
    Point()
        : x_(0.0), y_(0.0) {}

    Point(double x, double y)
        : x_(x), y_(y) {}

    double get_x() const { return x_; }
    double get_y() const { return y_; }

    double squaredDistance(const Point& other) const {
        // squared distance between this point and other
        double xDist = x_ - other.get_x();
        double yDist = y_ - other.get_y();
        return xDist*xDist + yDist*yDist;
    }

private:
    double x_, y_;
};

// FOR DEBUGGING
// std::ostream& operator<<(std::ostream& out, const Point& point) {
//     return out << std::fixed << std::setprecision(3) << "Point(" << point.get_x() << ", " << point.get_y() << ")";
// }


Point randomPoint(double x_min, double x_max, double y_min, double y_max) {
    auto x = R::runif(x_min, x_max);
    auto y = R::runif(y_min, y_max);
    return Point(x, y);
}

int randomIndex(int from, int to) {
    // random integer in range [from, to)
    return sample(to - from, 1)[0] + from - 1;
}

double randomAngle() {
    return R::runif(0, TWO_PI);
}

std::array<double, 2> randomDirection() {
    double angle = randomAngle();
    return std::array<double, 2>{std::sin(angle), std::cos(angle)};
}

double randomDistance(double from, double to) {
    return R::runif(from, to);
}


class Grid {
public:
    Grid(int nx, int ny)
        : grid_(nx*ny, 0),
          x_(nx),
          y_(ny){}

    int get(int x, int y) const {
        return grid_[array_index(x, y)];
    }

    void set(int x, int y, int value) {
        grid_[array_index(x, y)] = value;
    }

    int getMaxX() const { return x_; }
    int getMaxY() const { return y_; }

private:
    std::vector<int> grid_;
    int x_;
    int y_;

    int array_index(int x, int y) const {
        return y_ * x + y;
    }
};

bool isValid(const Point& candidate,
             double size_x,
             double size_y,
             double cellSize,
             double radius,
             const std::vector<Point>& pointsList,
             const Grid& grid) {
    if (candidate.get_x() >= 0 && candidate.get_x() < size_x &&
        candidate.get_y() >= 0 && candidate.get_y() < size_y) {
        int cellX = (int) (candidate.get_x() / cellSize);
        int cellY = (int) (candidate.get_y() / cellSize);
        int searchStartX = std::max(0, cellX - 2);
        int searchStartY = std::max(0, cellY - 2);
        int searchEndX = std::min(cellX + 2, grid.getMaxX() - 1);
        int searchEndY = std::min(cellY + 2, grid.getMaxY() - 1);

        for (int x = searchStartX; x <= searchEndX; ++x) {
            for (int y = searchStartY; y <= searchEndY; ++y) {
                int pointIndex = grid.get(x, y) - 1;
                if (pointIndex != -1) {
                    double sqrDst = candidate.squaredDistance(pointsList[pointIndex]);
                    if (sqrDst < radius*radius) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
    return false;
}

// Implementation pretty much mimics Sebastian Lague's Procedural Object Placement
// youtube vid: https://www.youtube.com/watch?v=7WcmyxyFO7o
std::vector<Point> generatePoints(double radius,
                                  double size_x,
                                  double size_y,
                                  int numSamplesBeforeRejection = 30) {
    double cellSize = radius / SQRT_2;

    int gridSizeX = (int)std::ceil(size_x / cellSize);
    int gridSizeY = (int)std::ceil(size_y / cellSize);
    Grid grid(gridSizeX, gridSizeY);

    std::vector<Point> points;
    std::vector<Point> spawnPoints;

    // Initialise spawnPoints with a random point
    spawnPoints.push_back(randomPoint(0, size_x, 0, size_y));

    while (!spawnPoints.empty()) {
        // Pick a random spawn point from the available list
        int spawnIndex = randomIndex(0, spawnPoints.size());
        Point spawnCentre = spawnPoints[spawnIndex];

        bool candidateAccepted = false;

        for (int i = 0; i < numSamplesBeforeRejection; ++i) {
            // Make a new candidate at a random distance (between r and 2r) away from
            // the spawn point, in a random direction...
            auto direction = randomDirection();
            auto distance = randomDistance(radius, 2*radius);
            Point candidate(spawnCentre.get_x() + direction[0] * distance,
                            spawnCentre.get_y() + direction[1] * distance);

            // ...if it doesn't collide with anything nearby in the grid...
            if (isValid(candidate, size_x, size_y, cellSize, radius, points, grid)) {

                // ...then accept it as a new point
                candidateAccepted = true;
                points.push_back(candidate);
                spawnPoints.push_back(candidate);

                grid.set((int)(candidate.get_x() / cellSize),
                         (int)(candidate.get_y() / cellSize),
                         (int)points.size());
                break;
            }
        }

        // If no accepted candidates after x tries, this spawn point is no good, so delete it
        if (!candidateAccepted) {
            spawnPoints.erase(spawnPoints.begin() + spawnIndex);
        }
    }
    return points;
}


// [[Rcpp::export]]
NumericMatrix poisson_disc_sampling(double radius,
                                    double size_x,
                                    double size_y,
                                    int num_samples_before_rejection = 30) {
    auto points = generatePoints(radius, size_x, size_y, num_samples_before_rejection);

    // The hardest part is extracting the data into an R matrix
    std::vector<double> xs_and_ys;
    std::transform(begin(points), end(points), std::back_inserter(xs_and_ys),
                   [](const Point& pt) { return pt.get_x(); });

    std::transform(begin(points), end(points), std::back_inserter(xs_and_ys),
                   [](const Point& pt) { return pt.get_y(); });

    NumericMatrix result(points.size(), 2, xs_and_ys.begin());

    return result;
}
