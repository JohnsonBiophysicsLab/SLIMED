#include "test_optimization_algorithm.hpp"

// Test for reseting NCG direction vector to (0, 0, 0) vectors
/*
TEST(NCGOptimizerTest, ResetNCGDirectionTest)
{
    Param param = Param();
    Mesh mesh = Mesh(param);
    Record record(param.maxIterations);
    Model model = Model(mesh, record);

    model.ncgDirection0[0].forceTotal.set(0, 0, 0.369);
    model.ncgDirection0[0].forceTotal.set(1, 0, 0.8666033);

    // Reset
    model.reset_ncg_direction();

    // Assert
    EXPECT_DOUBLE_EQ(model.ncgDirection0[0].forceTotal.get(0, 0), 0);
    EXPECT_DOUBLE_EQ(model.ncgDirection0[0].forceTotal.get(1, 0), 0);

    // Add additional assertions if needed
}
*/