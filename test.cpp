#include <math.h>
#include <iostream>

struct Coordinate
{
    double x;
    double y;
    double z;

    void assign(double givenX, double givenY, double givenZ)
    {
        x = givenX;
        y = givenY;
        z = givenZ;
    }

    void divideBy(double h){
        x /= h;
        y /= h;
        z /= h;
    }

    void clear(){
        x = 0;
        y = 0;
        z = 0;
    }

    void add(double a, double b, double c){
        x += a;
        y += b;
        z += c;
    }

    void print(){
        std::cout << ("x: ");
        std::cout << (x);
        std::cout << (", ");
        std::cout << ("y: ");
        std::cout << (y);
        std::cout << (", ");
        std::cout << ("z: ");
        std::cout << (z) << "\n";
    }
};

Coordinate knownCoordinates[3];

Coordinate unknownCoordinate;


// precompute constant data and store for each triplet of coordinates to do trilateration with.
struct TripletConstants{

    uint8_t coordinateIndexes[3];

    Coordinate e1, e2, e3;

    double h, i, j;

    double t; 

    void assign(const uint8_t left, const uint8_t center, const uint8_t right){
        coordinateIndexes[0] = left;
        coordinateIndexes[1] = center;
        coordinateIndexes[2] = right;
    }

    void precompute(){

        Coordinate * p1 = &knownCoordinates[coordinateIndexes[0]];
        Coordinate * p2 = &knownCoordinates[coordinateIndexes[1]];
        Coordinate * p3 = &knownCoordinates[coordinateIndexes[2]];

        // e1 = p2 - p1
        // e1 is the vector from p1 to p2
        e1.assign(
            p2->x - p1->x, 
            p2->y - p1->y, 
            p2->z - p1->z);

        // h = ||p2 - p1||
        h = sqrt( 
            e1.x * e1.x + 
            e1.y * e1.y + 
            e1.z * e1.z );

        // convert e1 into unit vector
        e1.divideBy(h);

        // i = e1 dotProduct (p3 - p1)
        i = e1.x * (p3->x - p1->x) + 
            e1.y * (p3->y - p1->y) + 
            e1.z * (p3->z - p1->z);


        // e2 = p3 - p1 - e1(e1 dotProduct (p3 - p1))
        // substitution
        // e2 = p3 - p1 - e1(i)                       
        e2.assign(
            p3->x - p1->x - (i * e1.x), 
            p3->y - p1->y - (i * e1.y),
            p3->z - p1->z - (i * e1.z));

        // convert e2 to unit vector
        t = sqrt(e2.x * e2.x + 
                 e2.y * e2.y + 
                 e2.z * e2.z );
        e2.divideBy(t);

        // don't bother checking if three fixed points are too close to being on the same line.
        // If t <= epsilon:
        //     Error: the three fixed points are too close to being on the same line.

        // j = e2 dotProduct (p3 - p1)
        j = e2.x * (p3->x - p1->x) 
          + e2.y * (p3->y - p1->y) 
          + e2.z * (p3->z - p1->z);

        // skip checking to see if points are too close to being on the same line
        // If j <= epsilon and j >= -epsilon:
        //     Error: the three fixed points are too close to being on the same line.
        
        // e3 = e1 crossProduct e2
        e3.assign(
            e1.y * e2.z - e1.z * e2.y,
            e1.z * e2.x - e1.x * e2.z,
            e1.x * e2.y - e2.x * e1.y);
    }
    
};

TripletConstants sensorTriplets[1];


void performTrilaterationIteration(){

    double d1 = 13.341664;
    double d2 = 14.3527;
    double d3 = 14.282857;

    // precomputed constants
    double h = sensorTriplets[0].h;
    double i = sensorTriplets[0].i;
    double j = sensorTriplets[0].j;
    Coordinate * e1 = &sensorTriplets[0].e1;
    Coordinate * e2 = &sensorTriplets[0].e2;
    Coordinate * e3 = &sensorTriplets[0].e3;

    Coordinate * p1 = &knownCoordinates[sensorTriplets[0].coordinateIndexes[0]];
    
    // computations
    double u = (d1*d1 - d2*d2 + h*h) / (2*h);
    double v = (d1*d1 - d3*d3 + i*(i - 2*u) + j*j) / (2*j);
    double w = sqrt(d1*d1 - u*u - v*v);

    // store the computation in unknownCoordinate
    unknownCoordinate.add(
        p1->x + u*e1->x + v*e2->x + w*e3->x,
        p1->y + u*e1->y + v*e2->y + w*e3->y,
        p1->z + u*e1->z + v*e2->z + w*e3->z
    );

}



int main(){

    knownCoordinates[0].assign(3.0, -4.0, 5.0);

    knownCoordinates[1].assign(1.0, 2.0, -3.0);

    knownCoordinates[2].assign(-6.0, 6.0, 6.0);

    sensorTriplets[0].assign(0, 1, 2);

    sensorTriplets[0].precompute();

    unknownCoordinate.clear();

    performTrilaterationIteration();

    // value should be (8,8,8)
    unknownCoordinate.print();

    return 0;
}