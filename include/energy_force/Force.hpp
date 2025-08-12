/**
 * @file Force.hpp
 * @author Y Ying (yying7@jh.edu)
 * @brief The Force class defines corresponding force terms
 *        type. Note that all the forces are gsl
 *        matrices with a size of (3,1).
 * @date 2023-01-05
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <sstream>
#include <ostream>
#include <omp.h>
#include <iomanip>
#include <algorithm>

#include "linalg/Linear_algebra.hpp"
#include "Parameters.hpp"


/**
 * @brief This class includes 7 components in calculating force
 *        exerted on vertrices: force due to membrane bending,
 *        area constraint, volume constraint, membrane thickness
 *        constraint, lipid tilting, regularization, and spline.
 *
 */
class Force
{
public:
    Matrix forceCurvature;      ///< Force due to membrane bending.
    Matrix forceArea;           ///< Force due to area constraint.
    Matrix forceVolume;         ///< Force due to volume constraint; set 0 if membrane is flat.
    Matrix forceThickness;      ///< Force due to membrane thickness change.
    Matrix forceTilt;           ///< Force due to lipid tilting.
    Matrix forceRegularization; ///< Regularization force; this does not change the
                                ///< geometry of actual limit surface. Instead, it
                                ///< restrains the shape of each mesh triangle to be
                                ///< as close to equilateral triangle as possible to avoid
                                ///< deformation of the triangular mesh.
    Matrix forceHarmonicBond;   ///< Force due to harmonic bond between the membrane
                                ///< and the attached spline, which can be a single particle
                                ///< or a lattice. The calculations are defined separately in
                                ///< spline_points.cpp due to it requiring the coordinate
                                ///< information of the spline points.
    Matrix forceTotal;          ///< Total force by summing up the 7 force terms.

    /**
     * @brief Construct a new Force object with all force components set to zero.
     */
    Force();

    /**
     * @brief Construct a new Force object with all force components
     * set to the same as the input force, i.e. create a deepcopy
     * of the input force.
     *
     * @param force
     */
    Force(const Force &force);

    /**
     * @brief Set forceTotal to the sum of 7 force components.
     * @return gsl_matrix of calculated total force.
     */
    Matrix &calculate_total_force();

    /**
     * @brief Get the magnitude (norm) of total force
     * @return double magnitude of the total force
     * @note This method assumes total force is up to date. If not, please use
     * calculate_total_force() to update the total force.
     *
     */
    double get_total_force_magnitude();

    /**
     * @brief Returns a string representation of the force values for various properties.
     *
     * @param force A `Force` object to extract the force values from.
     *
     * @return A string containing the force values for each property, separated by newlines.
     */
    std::string to_string_full();

    /**
     * @brief Returns a string representation of the forces acting on a vertex for each property.
     *
     * @param i The index of the vertex to extract the forces from.
     *
     * @return A string containing the force values for each property at the specified vertex, separated by newlines.
     */
    std::string to_string_full_at_vertex(const int i);

    /**
     * @brief Overrides the operator<< in ostream. Outputs
     * force in the format of (forceTotal[0], forceTotal[1], forceTotal[2]).
     * Used in command line output.
     *
     * @note Only outputs forceTotal. Please use Force::to_full_string() or
     * Force::to_full_string_at_vertex() for all components.
     */
    friend std::ostream &operator<<(std::ostream &stream, const Force &force);
};

/**
 * @brief This code defines a binary operator + for the Force class. The operator takes two 
 * const references to Force objects (lhs and rhs) as its arguments.
 * 
 * @param lhs 
 * @param rhs 
 * @return Force 
 */
Force operator+(const Force& lhs, const Force& rhs);

/**
 * @brief This code defines a compound assignment operator += for the Force class. 
 *        The operator takes a reference to a Force object (lhs) and a const reference 
 *        to another Force object (rhs) as its arguments. It modifies the lhs object by 
 *        adding the components of the rhs object to it.
 * 
 * @param lhs 
 * @param rhs 
 * @return Force& 
 */
Force& operator+=(Force& lhs, const Force& rhs);

