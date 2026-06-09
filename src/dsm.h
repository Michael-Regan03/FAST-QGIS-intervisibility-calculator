#include <vector>
#include <string>
#include "types.h"

#ifndef DSM_H
#define DSM_H

constexpr int  earth_diameter {12756000};


class DSM {
public:
    DSM(const std::string& filePath);
    
    int get_nrows() const;
    int get_ncols() const;
    double get_cellsize() const;
    double get_nodata() const;
    double get_xll() const;
    double get_yll() const;
    const std::vector<float>& get_data() const; 
    const std::vector<float>& get_los_errors() const;
    float windowSample(int obs_x, int obs_y, int offset_x, int offset_y) const ;
    const float indexAngle(int dx, int dy) const;
    const double index_mx_dist(types::PixelCoord offset) const;
    types::PixelCoord transformOffset(int dx, int dy, int step) const;
    void write_asc(const std::string& path, const std::vector<int>& result) const ;
    std::pair<int,int> geo_to_pix(double x, double y) const;
    int get_radius_pix() const;


private:
    std::vector<float> data {};
    int nrows {};
    int ncols {};
    double xll {};
    double yll {};
    double cellsize {}; 
    double nodata {};
    int windowSize{}; 
    int RadiusPix{};

    std::vector<float> mx_dist {};
    std::vector<float> mx_curvature {};
    std::vector<float> mx_angles {};
    // Bresenham tables
    std::vector<types::PixelCoord>  los_indices {};
    std::vector<float> los_errors {};
    std::vector<bool>  los_mask {};
    
    void load_asc(const std::string& path);
    void precompute(int radius_pix, bool curvature, float refraction);
    void bresenhamLoop(int radius_pix);
    float sample(int x, int y) const;
    
};


#endif