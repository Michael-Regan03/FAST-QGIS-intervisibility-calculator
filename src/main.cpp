
#include "dsm.h"
#include "types.h"
#include "interVisibility.h"
#include <iostream>



int main() {
    std::cout << "Loading DSM...\n";
    DSM dsm("assets/test.asc");
    std::cout << "DSM loaded\n";
    std::cout << "nrows: " << dsm.get_nrows() << " ncols: " << dsm.get_ncols() << "\n";

    std::vector<std::pair<double,double>> coords = {
        {715200.0, 733500.0},
        {715300.0, 733600.0},
        {715400.0, 733700.0},
    };

    types::ObserverPoints targets;
    for(auto& [x, y] : coords) {
        std::cout << "Converting " << x << ", " << y << "\n";
        auto [px, py] = dsm.geo_to_pix(x, y);
        std::cout << "Pixel: " << px << ", " << py << "\n";
        targets.push_back({x, y, px, py});
    }

    std::cout << "Starting calculation...\n";
    Intervisibilty iv(dsm, targets, 1.6);
    auto& result = iv.calculate();
    std::cout << "Calculation done\n";

    dsm.write_asc("assets/output.asc", result);
    std::cout << "Done\n";
    return 0;
}