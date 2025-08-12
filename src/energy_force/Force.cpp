#include "energy_force/Force.hpp"

/**
 * @brief Construct a new Force object with all force components set to zero.
 */
Force::Force() : forceCurvature(mat_calloc(3, 1)),
                 forceArea(mat_calloc(3, 1)),
                 forceVolume(mat_calloc(3, 1)),
                 forceThickness(mat_calloc(3, 1)),
                 forceTilt(mat_calloc(3, 1)),
                 forceRegularization(mat_calloc(3, 1)),
                 forceHarmonicBond(mat_calloc(3, 1)),
                 forceTotal(mat_calloc(3, 1))
{
}

/**
 * @brief Construct a new Force object with all force components
 * set to the same as the input force, i.e. create a deepcopy
 * of the input force.
 *
 * @param force
 */
Force::Force(const Force &force) : forceCurvature(force.forceCurvature),
                                   forceArea(force.forceArea),
                                   forceVolume(force.forceVolume),
                                   forceThickness(force.forceThickness),
                                   forceTilt(force.forceTilt),
                                   forceRegularization(force.forceRegularization),
                                   forceHarmonicBond(force.forceHarmonicBond),
                                   forceTotal(force.forceTotal)
{
}

/**
 * @brief Set forceTotal to the sum of 7 force components.
 * @return gsl_matrix of calculated total force.
 */
Matrix &Force::calculate_total_force()
{
    // zero the total force
    forceTotal.set_all(0);

    // add each term to the total force
    forceTotal += forceCurvature;
    forceTotal += forceArea;
    forceTotal += forceVolume;
    forceTotal += forceThickness;
    forceTotal += forceTilt;
    forceTotal += forceRegularization;
    forceTotal += forceHarmonicBond;

    return forceTotal;
}

/**
 * @brief Get the magnitude (norm) of total force
 * @return double magnitude of the total force
 * @note This method assumes total force is up to date. If not, please use
 * calculate_total_force() to update the total force.
 *
 */
double Force::get_total_force_magnitude(){
    return forceTotal.calculate_norm();
}

/**
 * @brief returns a string representation of the force values for various properties.
 *
 * @param force A `Force` object to extract the force values from.
 *
 * @return A string containing the force values for each property, separated by newlines.
 */
std::string Force::to_string_full()
{
    std::stringstream ss;
    ss << "force_curve = "
       << forceCurvature(0,0) << " , "
       << forceCurvature(1,0) << " , "
       << forceCurvature(2,0) << " , " << std::endl
       << "force_area = "
       << forceArea(0,0) << " , "
       << forceArea(1,0) << " , "
       << forceArea(2,0) << " , " << std::endl
       << "force_vol = "
       << forceVolume(0,0) << " , "
       << forceVolume(1,0) << " , "
       << forceVolume(2,0) << " , " << std::endl
       << "force_thic = "
       << forceThickness(0,0) << " , "
       << forceThickness(1,0) << " , "
       << forceThickness(2,0) << " , " << std::endl
       << "force_tilt = "
       << forceTilt(0,0) << " , "
       << forceTilt(1,0) << " , "
       << forceTilt(2,0) << " , " << std::endl
       << "force_reg = "
       << forceRegularization(0,0) << " , "
       << forceRegularization(1,0) << " , "
       << forceRegularization(2,0) << " , " << std::endl
       << "force_total = "
       << forceTotal(0,0) << " , "
       << forceTotal(1,0) << " , "
       << forceTotal(2,0) << " , " << std::endl;
    return ss.str();
}

/**
 * @brief Returns a string representation of the forces acting on a vertex for each property.
 *
 * @param i The index of the vertex to extract the forces from.
 *
 * @return A string containing the force values for each property at the specified vertex, separated by newlines.
 */
std::string Force::to_string_full_at_vertex(const int i)
{
    std::stringstream ss;
    ss << "force_curve @ " << i << " = "
       << forceCurvature(0,0) << " , "
       << forceCurvature(1,0) << " , "
       << forceCurvature(2,0) << " , " << std::endl
       << "force_area @ " << i << " = "
       << forceArea(0,0) << " , "
       << forceArea(1,0) << " , "
       << forceArea(2,0) << " , " << std::endl
       << "force_vol @ " << i << " = "
       << forceVolume(0,0) << " , "
       << forceVolume(1,0) << " , "
       << forceVolume(2,0) << " , " << std::endl
       << "force_thic @ " << i << " = "
       << forceThickness(0,0) << " , "
       << forceThickness(1,0) << " , "
       << forceThickness(2,0) << " , " << std::endl
       << "force_tilt @ " << i << " = "
       << forceTilt(0,0) << " , "
       << forceTilt(1,0) << " , "
       << forceTilt(2,0) << " , " << std::endl
       << "force_reg @ " << i << " = "
       << forceRegularization(0,0) << " , "
       << forceRegularization(1,0) << " , "
       << forceRegularization(2,0) << " , " << std::endl
       << "force_total @ " << i << " = "
       << forceTotal(0,0) << " , "
       << forceTotal(1,0) << " , "
       << forceTotal(2,0) << " , " << std::endl;
    return ss.str();
}

/**
 * @brief Overrides the operator<< in ostream. Outputs
 * force in the format of (forceTotal(0,0), forceTotal(1,0), forceTotal(2,0)).
 * Used in command line output.
 *
 * @note Only outputs forceTotal. Please use Force::to_full_string() or
 * Force::to_full_string_at_vertex() for all components.
 */
std::ostream &operator<<(std::ostream &stream, const Force &force)
{
    stream << force.forceTotal(0, 0) << ", "
           << force.forceTotal(1, 0) << ", "
           << force.forceTotal(2, 0);
    return stream;
}

/**
 * @brief This code defines a binary operator + for the Force class. The operator takes two 
 * const references to Force objects (lhs and rhs) as its arguments.
 * 
 * @param lhs 
 * @param rhs 
 * @return Force 
 */
Force operator+(const Force& lhs, const Force& rhs){
    Force sum(lhs);
    sum.forceArea += rhs.forceArea;
    sum.forceCurvature += rhs.forceCurvature;
    sum.forceHarmonicBond += rhs.forceHarmonicBond;
    sum.forceRegularization += rhs.forceRegularization;
    sum.forceThickness += rhs.forceThickness;
    sum.forceTilt += rhs.forceTilt;
    sum.forceVolume += rhs.forceVolume;
    sum.forceTotal += rhs.forceTotal;
    return sum;
}

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
Force& operator+=(Force& lhs, const Force& rhs){
    lhs.forceArea += rhs.forceArea;
    lhs.forceCurvature += rhs.forceCurvature;
    lhs.forceHarmonicBond += rhs.forceHarmonicBond;
    lhs.forceRegularization += rhs.forceRegularization;
    lhs.forceThickness += rhs.forceThickness;
    lhs.forceTilt += rhs.forceTilt;
    lhs.forceVolume += rhs.forceVolume;
    lhs.forceTotal += rhs.forceTotal;
    return lhs;
}