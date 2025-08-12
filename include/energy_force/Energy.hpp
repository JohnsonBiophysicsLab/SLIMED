/**
 * @file Energy.hpp
 * @author yying7@jh.edu
 * @brief This file defines the Energy class with 7 energy terms.
 * @date 2023-01-11
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#pragma once

#include <math.h>
#include <cmath>
#include <vector>
#include <string>

/**
 * @brief This class includes 7 components in calculating energy
 *        on vertrices including membrane bending energy,
 *        area constraint energy, volume constraint energy,
 *        membrane thickness constraint energy, lipid tilting energy,
 *        regularization energy, and harmonic bonding energy.
 * 
 */
class Energy{
public:

    double energyCurvature = 0.0; ///< Variable to store curvature energy
    double energyArea = 0.0; ///< Variable to store area constraint energy
    double energyVolume = 0.0; ///< Variable to store volume constraint energy
    double energyThickness = 0.0; ///< Variable to store thickness constraint energy
    double energyTilt = 0.0; ///< Variable to store tilt angle constraint energy
    double energyRegularization = 0.0; ///< Variable to store regularization energy
    double energyHarmonicBond = 0.0; ///< Variable to store harmonic bond energy
    double energyTotal = 0.0; ///< Variable to store total energy

    
    /**
    * @brief Generate a new Energy object with all energy terms
    * set to zero 
    */
    Energy();

    /**
     * @brief Construct a new Energy object by copying all member values
     * to the new Energy instance, i.e. deepcopy of input energy.
     * 
     * @param energy 
     */
    Energy(const Energy& energy);

    /**
     * @brief Set energyTotal to the sum of energy components.
     * @return Calculated total energy.
     */
    double calculateTotalEnergy();

    /**
     * @brief Overrides the operator << in ostream. Outputs
     * force in the format of (forceTotal[0], forceTotal[1], forceTotal[2]).
     * Used in command line output. 
     */
    friend std::ostream& operator<<(std::ostream& stream, const Energy& energy);
};

/**
 * @brief Add the energy components of two Energy objects together and update
 * the first Energy object with the result. 
 * 
 * @param e1 
 * @param e2 
 * @return Energy 
 */
Energy& operator+=(Energy& e1, const Energy& e2);