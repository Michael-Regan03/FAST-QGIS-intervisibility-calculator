// open a window of DSM data centred on the observer
// convert heights to slope angles relative to observer eye level
// for each target — draw a Bresenham line to it
// find max slope along that line — that's the blocking horizon
// compare target slope to horizon — record depth

#ifndef INTERVISIBILITY_H
#define INTERVISIBILITY_H



#include "dsm.h"
#include "types.h"
#include <vector>

class Intervisibilty {

public:
    Intervisibilty(const DSM& dsm, const types::ObserverPoints& points, double observer_height);
    std::vector<int>& calculate();

private:
    const DSM& dsm;
    const types::ObserverPoints& points;
    std::vector<int> result{};
    double observer_height{};
    

    bool checkLOS(int col, int row, int dx, int dy, double z_abs);

};

#endif