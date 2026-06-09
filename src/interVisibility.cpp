#include "interVisibility.h"
#include "types.h"
#include <vector>
#include <cmath>

Intervisibilty::Intervisibilty(const DSM& dsm, const types::ObserverPoints& points, double observer_height)
    : dsm(dsm), points(points), observer_height(observer_height) {
    
};



std::vector<int>& Intervisibilty::calculate(){
    

    result.resize(dsm.get_nrows() * dsm.get_ncols(), 0);


    if(points.empty()){
        return result;
    }



    for(auto row =0;  row < dsm.get_nrows(); ++row){
        for(auto col=0; col < dsm.get_ncols(); ++col ){
             // setup — done once per observer pixel
            float terrain = dsm.windowSample(col, row, 0, 0);
            if(std::isnan(terrain)) 
                continue;
            double z_abs = terrain + observer_height;
            for(auto target: points){
                int dx = target.pix_x - col;
                int dy = target.pix_y - row;
                
                if (checkLOS(col, row, dx, dy, z_abs)){
                    result[row * dsm.get_ncols() + col]++;
                }
            }
        }
        
    }

    return result;
};




bool Intervisibilty::checkLOS(int col, int row, int dx, int dy, double z_abs){
    
    if(dx == 0 && dy == 0) return false;

    // skip if target is outside analysis radius
    if(std::abs(dx) > dsm.get_radius_pix() || 
       std::abs(dy) > dsm.get_radius_pix()) 
        return false;

    float target_terrain = dsm.windowSample(col, row, dx, dy);
    if(std::isnan(target_terrain)) return false;

    types::PixelCoord target_offset = {std::abs(dx), std::abs(dy)};
    double target_dist = dsm.index_mx_dist(target_offset);
    double target_slope = (target_terrain - z_abs) / target_dist;

    bool steep = std::abs(dy) > std::abs(dx);
    int steps = steep ? std::abs(dy) : std::abs(dx);
    double max_slope = -std::numeric_limits<double>::infinity();


    for(int step = 0; step < steps; ++step) {
        
        types::PixelCoord offset = dsm.transformOffset(dx, dy, step);
        float terrain = dsm.windowSample(col, row, offset.x, offset.y);
        if(std::isnan(terrain)) continue;

        const double dist = dsm.index_mx_dist(offset);



        // calculate slope, update horizon
        double slope = (terrain - z_abs) / dist;


        if(slope > max_slope) max_slope = slope;
        
    }

    if(target_slope >= max_slope) return true;
    return false;


}