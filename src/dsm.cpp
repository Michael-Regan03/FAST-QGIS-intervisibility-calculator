#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <cmath>
#include "dsm.h"


static constexpr int DEFAULT_RADIUS_PIX = 200;


DSM::DSM(const std::string& filePath){
        load_asc(filePath);
        precompute(DEFAULT_RADIUS_PIX, false, 0);

};

int DSM::get_nrows()       const { return nrows; };
int DSM::get_ncols()       const { return ncols; };
double DSM::get_cellsize() const { return cellsize; };
double DSM::get_nodata()   const { return nodata; };
double DSM::get_xll()      const { return xll; };
double DSM::get_yll()      const { return yll; };
const std::vector<float>& DSM::get_los_errors() const { return los_errors; };

const std::vector<float>& DSM::get_data() const { return data; };





void DSM::load_asc(const std::string& path){
        
        std::ifstream file(path);
        
        if (!file.is_open())
            throw std::runtime_error("ERROR!!! cannot open file: " + path);
        

        bool has_xll   = false;
        bool is_center = false;        
    
        std::string line{};
        int headers_found {};

        
        
        while(headers_found < 6 && std::getline(file,line)) {

            std::istringstream ss{line};
            std::string key;
            ss >> key;


            if (key == "ncols" ) { ss >> ncols; headers_found++; }
            else if (key == "nrows")        { ss >> nrows;    headers_found++; }
            else if (key == "xllcorner")    { ss >> xll;      headers_found++; has_xll = true;  is_center = false; }
            else if (key == "xllcenter")    { ss >> xll;      headers_found++; has_xll = true;  is_center = true;  }
            else if (key == "yllcorner")    { ss >> yll;      headers_found++; }
            else if (key == "yllcenter")    { ss >> yll;      headers_found++; is_center = true; }
            else if (key == "cellsize")     { ss >> cellsize; headers_found++; }
            else if (key == "nodata_value") { ss >> nodata;   headers_found++; }
    
        }

        if (ncols <= 0 || nrows <= 0 || cellsize <= 0.0)
            throw std::runtime_error("ERROR!!! DSM is invalid or incomplete header in " + path);
        

        if (!has_xll)
            throw std::runtime_error("ERROR!!! Missing xll coordinate in " + path);

        // xllcenter means the coord is the  of the corner cell,
        // convert to corner so all downstream code is consistent
        if (is_center) {
            xll -= cellsize * 0.5;
            yll -= cellsize * 0.5;
        }

        
        data.reserve(static_cast<size_t>(nrows) * ncols);
        float val{};
        size_t expected {static_cast<size_t>(nrows) * ncols};

        
         while (data.size() < expected && file >> val) {
            // Normalise nodata to NaN so comparisons are clean later
            data.push_back(val == static_cast<float>(nodata)
                           ? std::numeric_limits<float>::quiet_NaN()
                           : val);
        }

        if (data.size() != expected) {
            throw std::runtime_error(
                "DSM: expected " + std::to_string(expected) +
                " values, got "  + std::to_string(data.size()) +
                " in " + path);
        }


    };

void DSM::precompute(int radius_pix, bool curvature, float refraction){
        int window_size{2 * radius_pix + 1};
        RadiusPix = radius_pix;
        windowSize = window_size;
        mx_dist.resize(window_size * window_size);
        mx_angles.resize(window_size*window_size);
        for(auto r{0}; r < window_size; ++r){
            for(auto c{0}; c < window_size; ++c){
                int index {r *window_size + c};
                
                double a {static_cast<double>(c) - radius_pix};
                //double b {abs(r - radius_pix)};
                double b {static_cast<double>(radius_pix) - r};
                
                
                double angle {std::atan2(b,a) * (180/ M_PI)};
                if(angle < 0) {angle += 360;}
                       
                mx_dist[index] = std::sqrt(a*a + b*b); //pythagoras theorem
                mx_angles[index] = angle;
            }
        }


        
        mx_curvature.resize(window_size*window_size, 0.0f);
        if(curvature){
            for(int i = 0; i < window_size * window_size; ++i){
                double real_distance = mx_dist[i] * cellsize;
                double correction = (real_distance * real_distance )/earth_diameter * (1-refraction);
                mx_curvature[i] = correction;
            }
        }

        bresenhamLoop(radius_pix);

    };


void DSM::bresenhamLoop(int radius_pix){
        if (radius_pix <= 0)
            throw std::runtime_error("ERROR!!! radius_pix must be greater than zero");
        
        los_indices.resize((radius_pix + 1 ) * radius_pix);
        los_errors.resize((radius_pix + 1 ) * radius_pix);
        int dx = radius_pix;
        for(int dy{0}; dy <= radius_pix; ++dy){
            
            int D{};
            int x{};
            int y{};

            for(int step = 0; step < radius_pix; ++step){
                // bresenham step
                x += 1;

                D += dy;

                if(2 * D >= dx){
                    y += 1;
                    D -= dx;
                }

                double error_fraction = static_cast<double>(D)/dx;
                int index = dy * radius_pix + step;
                los_errors[index] = error_fraction;
                los_indices[index] = {x, y};
            }
      


        }
        los_mask.resize((radius_pix + 1) * radius_pix, false);
        for(int col{0}; col < radius_pix; ++col){
            int best {0};
            double minimum {1.0};
            for(int row = 0; row <= radius_pix; ++row){
                int index = row * radius_pix + col;
                if(los_errors[index] < minimum){
                    minimum = los_errors[index];
                    best = row;
                }
            }
            int index = best * radius_pix + col;

            los_mask[index] = true;
        }
    };



    
    




std::pair<int,int> DSM::geo_to_pix(double x, double y) const{
        int col {(x - xll) / cellsize};
        int row {(yll + nrows * cellsize - y) / cellsize};


        if (col < 0 || col >= ncols || row < 0 || row >= nrows)
            throw std::runtime_error(
                "ERROR!!! Coordinate outside raster extent: (" +
                std::to_string(x) + ", " + std::to_string(y) + ")");

        return {col, row};
    };



float DSM::sample(int x, int y) const {
        if (x >= ncols || x < 0)
            throw std::runtime_error("ERROR!!! Out of bounds error x: " 
                + std::to_string(x) + " nrows: " 
                + std::to_string(nrows));

        if (y >= nrows || y < 0)
            throw std::runtime_error("ERROR!!! Out of bounds error y: " 
                + std::to_string(y) + " ncols: " 
                + std::to_string(ncols));


        return data[y * ncols + x];
    };


float DSM::windowSample(int obs_x, int obs_y, int offset_x, int offset_y) const {
    int col = obs_x + offset_x;
    int row = obs_y + offset_y;

    if (col < 0 || col >= ncols || row < 0 || row >= nrows)
        return std::numeric_limits<float>::quiet_NaN();

    return data[row * ncols + col];
}


// octant 1: x+, y+, not steep  →  { x,  y}   (stored directly)
// octant 2: x+, y+, steep      →  { y,  x}   (swap)
// octant 3: x-, y+, not steep  →  {-x,  y}   (negate x)
// octant 4: x-, y+, steep      →  { y, -x}   (swap, negate x)
// octant 5: x-, y-, not steep  →  {-x, -y}   (negate both)
// octant 6: x-, y-, steep      →  {-y, -x}   (swap, negate both)
// octant 7: x+, y-, not steep  →  { x, -y}   (negate y)
// octant 8: x+, y-, steep      →  {-y,  x}   (swap, negate y)

types::PixelCoord DSM::transformOffset(int dx, int dy, int step) const {
    bool flip_x = dx < 0;
    bool flip_y = dy < 0;
    bool steep  = std::abs(dy) > std::abs(dx);

    // which row in the bresenham table
    int table_row = steep ? std::abs(dx) : std::abs(dy);

    // get the stored first-octant offset
    int index = table_row * RadiusPix + step;
    types::PixelCoord offset = los_indices[index];

    // swap if steep
    if(steep) std::swap(offset.x, offset.y);

    // apply sign flips
    if(flip_x) offset.x = -offset.x;
    if(flip_y) offset.y = -offset.y;

    return offset;
}

const float DSM::indexAngle(int dx, int dy) const { 
    bool flip_x = dx < 0;
    bool flip_y = dy < 0;
    bool steep  = std::abs(dy) > std::abs(dx);

    int tx = steep ? std::abs(dy) : std::abs(dx);
    int ty = steep ? std::abs(dx) : std::abs(dy);

    if (flip_x) tx = -tx;
    if (flip_y) ty = -ty;

    return mx_angles[(RadiusPix + ty) * windowSize + (RadiusPix + tx)];
}


void DSM::write_asc(const std::string& path, 
                    const std::vector<int>& result) const {
    std::ofstream file(path);
    
    if (!file.is_open())
        throw std::runtime_error("Cannot open output file: " + path);

    // write header
    file << "ncols "         << ncols     << "\n";
    file << "nrows "         << nrows     << "\n";
    file << "xllcorner "     << xll       << "\n";
    file << "yllcorner "     << yll       << "\n";
    file << "cellsize "      << cellsize  << "\n";
    file << "nodata_value "  << -9999     << "\n";

    // write data row by row
    for(int row = 0; row < nrows; ++row) {
        for(int col = 0; col < ncols; ++col) {
            file << result[row * ncols + col];
            if(col < ncols - 1) file << " ";
        }
        file << "\n";
    }
}


const double DSM::index_mx_dist(types::PixelCoord offset) const{

    return mx_dist[(RadiusPix + offset.y) * windowSize + (RadiusPix + offset.x)];
}

int DSM::get_radius_pix() const { return RadiusPix; }
