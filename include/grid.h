// PolyMorph
// Copyright (c) 2024-2025
// Nicolas Müller, nicolmueller@ethz.ch
// Roman Vetter, vetterro@ethz.ch
// ETH Zurich

#ifndef GRID_H
#define GRID_H

#include <vector>
#include <fstream>
#include <cmath>
#include <sstream>
#include <iostream>

#include "geometry.h"

template<typename T>
struct Grid {
    std::vector<std::vector<T>> data; // 2D grid of type T
    Grid (size_t Nx, size_t Ny) : data(Nx, std::vector<T>(Ny)) {}
    Grid (size_t Nx, size_t Ny, T value) : data(Nx, std::vector<T>(Ny, value)) {}
    Grid () {}
    T& operator()(int i, int j) { return data[i][j]; }
    T& operator()(Index idx) { return data[idx.i][idx.j]; }
    T operator()(int i, int j) const { return data[i][j]; }
    T operator()(Index idx) const { return data[idx.i][idx.j]; }
    size_t sizeX() const { return data.size(); }
    size_t sizeY() const { return data.empty() ? 0 : data[0].size(); }
    size_t sizeZ() const; // only defined for T=vector
    void parallel_copy_from(const Grid<T>& other); // parallelized assignment operator
    std::string to_vtk(std::string name); // convert to VTK format (DataArray)
    void rescale(size_t Nx, size_t Ny, int offset_x, int offset_y, T fill_value);
};

template<typename T>
size_t Grid<T>::sizeZ() const {
    return 0;
}

template<>
size_t Grid<std::vector<double>>::sizeZ() const {
    return data[0][0].size();
}

template<>
size_t Grid<std::vector<Point>>::sizeZ() const {
    return data[0][0].size();
}

// parallelized assignment operator for large grids of same size
template<typename T>
void Grid<T>::parallel_copy_from(const Grid<T>& other) {
    if (this != &other) {
        #pragma omp parallel for
        for (size_t i = 0; i < data.size(); i++) {
            data[i] = other.data[i]; // note: we are copying vectors here
        }
    }
}

template<typename T> // for scalar grids
std::string Grid<T>::to_vtk(std::string name) {
    std::stringstream xml;
    xml << "<DataArray type=\"Float64\" Name=\"" << name 
    << "\" format=\"ascii\">" << std::endl;
    for (int i = 0; i < sizeX(); i++) {
        for (int j = 0; j < sizeY(); j++) {
            xml << data[i][j] << " ";
        }
        xml << std::endl;
    }
    xml << "</DataArray>" << std::endl;
    return xml.str();
}

template<> // for std-vector grids
std::string Grid<std::vector<double>>::to_vtk(std::string name) {
    std::stringstream xml;
    xml << "<DataArray type=\"Float64\" Name=\"" << name 
        << "\" NumberOfComponents=\"" << sizeZ() 
        << "\" format=\"ascii\">" << std::endl;
    for (int i = 0; i < sizeX(); i++) {
        for (int j = 0; j < sizeY(); j++) {
            for (int k = 0; k < sizeZ(); k++) {
                xml << data[i][j][k] << " ";
            }
        }
        xml << std::endl;
    }
    xml << "</DataArray>" << std::endl;
    
    return xml.str();
}

template<> // for Point (2D vector) grids
std::string Grid<Point>::to_vtk(std::string name) {
    std::stringstream xml;
    xml << "<DataArray type=\"Float64\" Name=\"" << name 
        << "\" NumberOfComponents=\"2\" format=\"ascii\">" << std::endl;
    for (int i = 0; i < sizeX(); i++) {
        for (int j = 0; j < sizeY(); j++) {
            xml << data[i][j].x << " " << data[i][j].y << " ";
        }
        xml << std::endl;
    }
    xml << "</DataArray>" << std::endl;
    
    return xml.str();
}

template<> // for gradient (vector of 2D vectors) grids
std::string Grid<std::vector<Point>>::to_vtk(std::string name) {
    std::stringstream xml;
    xml << "<DataArray type=\"Float64\" Name=\"" << name 
        << "_x\" NumberOfComponents=\"" << sizeZ() 
        << "\" format=\"ascii\">" << std::endl;
    for (int i = 0; i < sizeX(); i++) {
        for (int j = 0; j < sizeY(); j++) {
            for (int k = 0; k < sizeZ(); k++) {
                xml << data[i][j][k].x << " ";
            }
        }
        xml << std::endl;
    }
    xml << "</DataArray>" << std::endl;
    xml << "<DataArray type=\"Float64\" Name=\"" << name 
        << "_y\" NumberOfComponents=\"" << sizeZ() 
        << "\" format=\"ascii\">" << std::endl;
    for (int i = 0; i < sizeX(); i++) {
        for (int j = 0; j < sizeY(); j++) {
            for (int k = 0; k < sizeZ(); k++) {
                xml << data[i][j][k].y << " ";
            }
        }
        xml << std::endl;
    }
    xml << "</DataArray>" << std::endl;
    return xml.str();
}

#endif