#include "mesh/Gauss_quadrature.hpp"

void get_gauss_quadrature_weight_VWU(const int &gaussQuadratureN,
                                     Matrix &vwu,
                                     Matrix &weight)
{                                 // To setup the Gaussian points for integral calculation. 'n' here is the Gausssian order.
    std::vector<double> vwuCoeff; // weight
    switch (gaussQuadratureN)
    {
    case 1:
    {
        // VMU
        std::vector<std::vector<double>> vwuTemp = {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}};
        vwu = Matrix(vwuTemp);
        // weight
        vwuCoeff = {1.0};
    }
    break;
    case 2:
    {

        std::vector<std::vector<double>> vwuTemp = {{1.0 / 6.0, 1.0 / 6.0, 4.0 / 6.0},
                                                    {1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0},
                                                    {4.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0}};

        vwu = Matrix(vwuTemp);

        vwuCoeff = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
    }
    break;
    case 3: // third order Guasssian may not converge. Weird! -Yiben
    {
        std::vector<std::vector<double>> vwuTemp = {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
                                                    {1.0 / 5.0, 1.0 / 5.0, 3.0 / 5.0},
                                                    {1.0 / 5.0, 3.0 / 5.0, 1.0 / 5.0},
                                                    {3.0 / 5.0, 1.0 / 5.0, 1.0 / 5.0}};
        vwuCoeff = {-0.56250000000000, 0.52083333333333, 0.52083333333333, 0.52083333333333};
    }
    break;
    case 4:
    {
        std::vector<std::vector<double>> vwuTemp = {{0.44594849091597, 0.44594849091597, 0.10810301816807},
                                                    {0.44594849091597, 0.10810301816807, 0.44594849091597},
                                                    {0.10810301816807, 0.44594849091597, 0.44594849091597},
                                                    {0.09157621350977, 0.09157621350977, 0.81684757298046},
                                                    {0.09157621350977, 0.81684757298046, 0.09157621350977},
                                                    {0.81684757298046, 0.09157621350977, 0.09157621350977}};
        vwu = Matrix(vwuTemp);
        vwuCoeff = {0.22338158967801, 0.22338158967801, 0.22338158967801, 0.10995174365532, 0.10995174365532, 0.10995174365532};
    }
    break;
    case 5:
    {
        std::vector<std::vector<double>> vwuTemp = {{0.33333333333333, 0.33333333333333, 0.33333333333333},
                                                    {0.47014206410511, 0.47014206410511, 0.05971587178977},
                                                    {0.47014206410511, 0.05971587178977, 0.47014206410511},
                                                    {0.05971587178977, 0.47014206410511, 0.47014206410511},
                                                    {0.10128650732346, 0.10128650732346, 0.79742698535309},
                                                    {0.10128650732346, 0.79742698535309, 0.10128650732346},
                                                    {0.79742698535309, 0.10128650732346, 0.10128650732346}};
        vwu = Matrix(vwuTemp);
        vwuCoeff = {0.22500000000000, 0.13239415278851, 0.13239415278851, 0.13239415278851, 0.12593918054483, 0.12593918054483, 0.12593918054483};
    }
    break;
    case 6:
    {
        std::vector<std::vector<double>> vwuTemp = {{0.24928674517091, 0.24928674517091, 0.50142650965818},
                                                    {0.24928674517091, 0.50142650965818, 0.24928674517091},
                                                    {0.50142650965818, 0.24928674517091, 0.24928674517091},
                                                    {0.06308901449150, 0.06308901449150, 0.87382197101700},
                                                    {0.06308901449150, 0.87382197101700, 0.06308901449150},
                                                    {0.87382197101700, 0.06308901449150, 0.06308901449150},
                                                    {0.31035245103378, 0.63650249912140, 0.05314504984482},
                                                    {0.63650249912140, 0.05314504984482, 0.31035245103378},
                                                    {0.05314504984482, 0.31035245103378, 0.63650249912140},
                                                    {0.63650249912140, 0.31035245103378, 0.05314504984482},
                                                    {0.31035245103378, 0.05314504984482, 0.63650249912140},
                                                    {0.05314504984482, 0.63650249912140, 0.31035245103378}};
        vwu = Matrix(vwuTemp);
        vwuCoeff = {0.11678627572638, 0.11678627572638, 0.11678627572638, 0.05084490637021, 0.05084490637021, 0.05084490637021, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837, 0.08285107561837};
    }
    break;
    default: // default to simplest middle point VMU
    {        // VMU
        std::vector<std::vector<double>> vwuTemp = {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}};
        vwu = Matrix(vwuTemp);
        // weight
        vwuCoeff = {1.0};
    }
    break;
    }

    // set weight
    weight = Matrix(vwuCoeff.size(), 1);
    for (int i = 0; i < vwuCoeff.size(); i++)
    {
        weight.set(i, 0, vwuCoeff[i]);
    }
}

void get_shapefunction(const Matrix &vwu, Matrix &shapefunction)
{
    // get vwu
    double v = vwu.get(0, 0);
    double w = vwu.get(0, 1);
    double u = vwu.get(0, 2);
    /*
     *   shape_functions(0,:), shape functions;
     *   shape_functions(1,:), differential to v;
     *   shape_functions(2,:), differential to w;
     *   shape_functions(3,:), double differential to v;
     *   shape_functions(4,:), double differential to w;
     *   shape_functions(5,:), differential to v and w;
     *   shape_functions(6,:), differential to w and v;
     */
    shapefunction.set(0, 0, 1.0 / 12.0 * (pow(u, 4.0) + 2.0 * pow(u, 3.0) * v));
    shapefunction.set(1, 0, 1.0 / 12.0 * (-2.0 * pow(u, 3.0) - 6.0 * pow(u, 2.0) * v));
    shapefunction.set(2, 0, 1.0 / 12.0 * (-4.0 * pow(u, 3.0) - 6.0 * pow(u, 2.0) * v));
    shapefunction.set(3, 0, u * v);
    shapefunction.set(4, 0, pow(u, 2.0) + u * v);
    shapefunction.set(5, 0, 1.0 / 2.0 * (pow(u, 2.0) + 2.0 * u * v));
    shapefunction.set(6, 0, 1.0 / 2.0 * (pow(u, 2.0) + 2.0 * u * v));

    shapefunction.set(0, 1, 1.0 / 12.0 * (pow(u, 4.0) + 2.0 * pow(u, 3.0) * w));
    shapefunction.set(1, 1, 1.0 / 12.0 * (-4.0 * pow(u, 3.0) - 6.0 * pow(u, 2.0) * w));
    shapefunction.set(2, 1, 1.0 / 12.0 * (-2.0 * pow(u, 3.0) - 6.0 * pow(u, 2.0) * w));
    shapefunction.set(3, 1, pow(u, 2.0) + u * w);
    shapefunction.set(4, 1, u * w);
    shapefunction.set(5, 1, 1.0 / 2.0 * (pow(u, 2.0) + 2.0 * u * w));
    shapefunction.set(6, 1, 1.0 / 2.0 * (pow(u, 2.0) + 2.0 * u * w));

    shapefunction.set(0, 2, 1.0 / 12.0 * (pow(u, 4.0) + 2.0 * pow(u, 3.0) * w + 6.0 * pow(u, 3.0) * v + 6.0 * pow(u, 2.0) * v * w + 12.0 * pow(u, 2.0) * pow(v, 2.0) + 6.0 * u * pow(v, 2.0) * w + 6.0 * u * pow(v, 3.0) + 2.0 * pow(v, 3.0) * w + pow(v, 4.0)));
    shapefunction.set(1, 2, 1.0 / 12.0 * (2.0 * pow(u, 3.0) + 6.0 * pow(u, 2.0) * v - 6.0 * u * pow(v, 2.0) - 2.0 * pow(v, 3.0)));
    shapefunction.set(2, 2, 1.0 / 12.0 * (-2.0 * pow(u, 3.0) - 6.0 * pow(u, 2.0) * w - 12.0 * pow(u, 2) * v - 12.0 * u * v * w - 18.0 * u * pow(v, 2.0) - 6.0 * pow(v, 2.0) * w - 4.0 * pow(v, 3.0)));
    shapefunction.set(3, 2, -2.0 * u * v);
    shapefunction.set(4, 2, u * w + v * w + u * v + pow(v, 2.0));
    shapefunction.set(5, 2, 1.0 / 2.0 * (-pow(u, 2.0) - 2.0 * u * v + pow(v, 2.0)));
    shapefunction.set(6, 2, 1.0 / 2.0 * (-pow(u, 2.0) - 2.0 * u * v + pow(v, 2.0)));

    shapefunction.set(0, 3, 1.0 / 12.0 * (6.0 * pow(u, 4.0) + 24.0 * pow(u, 3.0) * w + 24.0 * pow(u, 2.0) * pow(w, 2.0) + 8.0 * u * pow(w, 3.0) + pow(w, 4.0) + 24.0 * pow(u, 3.0) * v + 60.0 * pow(u, 2.0) * v * w + 36.0 * u * v * pow(w, 2.0) + 6.0 * v * pow(w, 3.0) + 24.0 * pow(u, 2.0) * pow(v, 2.0) + 36.0 * u * pow(v, 2.0) * w + 12.0 * pow(v, 2.0) * pow(w, 2.0) + 8.0 * u * pow(v, 3.0) + 6.0 * pow(v, 3.0) * w + pow(v, 4.0)));
    shapefunction.set(1, 3, 1.0 / 12.0 * (-12.0 * pow(u, 2.0) * w - 12.0 * u * pow(w, 2.0) - 2.0 * pow(w, 3.0) - 24.0 * pow(u, 2.0) * v - 48.0 * u * v * w - 12.0 * v * pow(w, 2.0) - 24.0 * u * pow(v, 2.0) - 18.0 * pow(v, 2.0) * w - 4.0 * pow(v, 3.0)));
    shapefunction.set(2, 3, 1.0 / 12.0 * (-24.0 * pow(u, 2.0) * w - 24.0 * u * pow(w, 2.0) - 4.0 * pow(w, 3.0) - 12.0 * pow(u, 2.0) * v - 48.0 * u * v * w - 18.0 * v * pow(w, 2.0) - 12.0 * u * pow(v, 2.0) - 12.0 * pow(v, 2.0) * w - 2.0 * pow(v, 3.0)));
    shapefunction.set(3, 3, -2.0 * u * w - 2.0 * pow(u, 2.0) + v * w + pow(v, 2.0));
    shapefunction.set(4, 3, -2.0 * pow(u, 2.0) + pow(w, 2.0) - 2.0 * u * v + v * w);
    shapefunction.set(5, 3, 1.0 / 2.0 * (-2.0 * pow(u, 2.0) + pow(w, 2.0) + 4.0 * v * w + pow(v, 2.0)));
    shapefunction.set(6, 3, 1.0 / 2.0 * (-2.0 * pow(u, 2.0) + pow(w, 2.0) + 4.0 * v * w + pow(v, 2.0)));

    shapefunction.set(0, 4, 1.0 / 12.0 * (pow(u, 4.0) + 6.0 * pow(u, 3.0) * w + 12.0 * pow(u, 2.0) * pow(w, 2.0) + 6.0 * u * pow(w, 3.0) + pow(w, 4.0) + 2.0 * pow(u, 3.0) * v + 6.0 * pow(u, 2.0) * v * w + 6.0 * u * v * pow(w, 2.0) + 2.0 * v * pow(w, 3.0)));
    shapefunction.set(1, 4, 1.0 / 12.0 * (-2.0 * pow(u, 3.0) - 12.0 * pow(u, 2.0) * w - 18.0 * u * pow(w, 2.0) - 4.0 * pow(w, 3.0) - 6.0 * pow(u, 2.0) * v - 12.0 * u * v * w - 6.0 * v * pow(w, 2.0)));
    shapefunction.set(2, 4, 1.0 / 12.0 * (2.0 * pow(u, 3.0) + 6.0 * pow(u, 2.0) * w - 6.0 * u * pow(w, 2.0) - 2.0 * pow(w, 3.0)));
    shapefunction.set(3, 4, u * w + pow(w, 2.0) + u * v + v * w);
    shapefunction.set(4, 4, -2.0 * u * w);
    shapefunction.set(5, 4, 1.0 / 2.0 * (-pow(u, 2.0) - 2.0 * u * w + pow(w, 2.0)));
    shapefunction.set(6, 4, 1.0 / 2.0 * (-pow(u, 2.0) - 2.0 * u * w + pow(w, 2.0)));

    shapefunction.set(0, 5, 1.0 / 12.0 * (2.0 * u * pow(v, 3.0) + pow(v, 4.0)));
    shapefunction.set(1, 5, 1.0 / 12.0 * (6.0 * u * pow(v, 2.0) + 2.0 * pow(v, 3.0)));
    shapefunction.set(2, 5, -1.0 / 6.0 * pow(v, 3.0));
    shapefunction.set(3, 5, u * v);
    shapefunction.set(4, 5, 0.0);
    shapefunction.set(5, 5, -1.0 / 2.0 * pow(v, 2.0));
    shapefunction.set(6, 5, -1.0 / 2.0 * pow(v, 2.0));

    shapefunction.set(0, 6, 1.0 / 12.0 * (pow(u, 4.0) + 6.0 * pow(u, 3.0) * w + 12.0 * pow(u, 2.0) * pow(w, 2.0) + 6.0 * u * pow(w, 3.0) + pow(w, 4.0) + 8.0 * pow(u, 3.0) * v + 36.0 * pow(u, 2.0) * v * w + 36.0 * u * v * pow(w, 2.0) + 8.0 * v * pow(w, 3.0) + 24.0 * pow(u, 2.0) * pow(v, 2.0) + 60.0 * u * pow(v, 2.0) * w + 24.0 * pow(v, 2.0) * pow(w, 2.0) + 24.0 * u * pow(v, 3.0) + 24.0 * pow(v, 3.0) * w + 6.0 * pow(v, 4.0)));
    shapefunction.set(1, 6, 1.0 / 12.0 * (4.0 * pow(u, 3.0) + 18.0 * pow(u, 2.0) * w + 12.0 * u * pow(w, 2.0) + 2.0 * pow(w, 3.0) + 24.0 * pow(u, 2.0) * v + 48.0 * u * v * w + 12.0 * v * pow(w, 2.0) + 24.0 * u * pow(v, 2.0) + 12.0 * pow(v, 2.0) * w));
    shapefunction.set(2, 6, 1.0 / 12.0 * (2.0 * pow(u, 3.0) + 6.0 * pow(u, 2.0) * w - 6.0 * u * pow(w, 2.0) - 2.0 * pow(w, 3.0) + 12.0 * pow(u, 2.0) * v - 12.0 * v * pow(w, 2.0) + 12.0 * u * pow(v, 2.0) - 12.0 * pow(v, 2.0) * w));
    shapefunction.set(3, 6, pow(u, 2.0) + u * w - 2.0 * v * w - 2.0 * pow(v, 2.0));
    shapefunction.set(4, 6, -2.0 * u * w - 2.0 * u * v - 2.0 * v * w - 2.0 * pow(v, 2.0));
    shapefunction.set(5, 6, 1.0 / 2.0 * (pow(u, 2.0) - 2.0 * u * w - pow(w, 2.0) - 4.0 * v * w - 2.0 * pow(v, 2.0)));
    shapefunction.set(6, 6, 1.0 / 2.0 * (pow(u, 2.0) - 2.0 * u * w - pow(w, 2.0) - 4.0 * v * w - 2.0 * pow(v, 2.0)));

    shapefunction.set(0, 7, 1.0 / 12.0 * (pow(u, 4.0) + 8.0 * pow(u, 3.0) * w + 24.0 * pow(u, 2.0) * pow(w, 2.0) + 24.0 * u * pow(w, 3.0) + 6.0 * pow(w, 4.0) + 6.0 * pow(u, 3.0) * v + 36.0 * pow(u, 2.0) * v * w + 60.0 * u * v * pow(w, 2.0) + 24.0 * v * pow(w, 3.0) + 12.0 * pow(u, 2.0) * pow(v, 2.0) + 36.0 * u * pow(v, 2.0) * w + 24.0 * pow(v, 2.0) * pow(w, 2.0) + 6.0 * u * pow(v, 3.0) + 8.0 * pow(v, 3.0) * w + pow(v, 4.0)));
    shapefunction.set(1, 7, 1.0 / 12.0 * (2.0 * pow(u, 3.0) + 12.0 * pow(u, 2.0) * w + 12.0 * u * pow(w, 2.0) + 6.0 * pow(u, 2.0) * v - 12.0 * v * pow(w, 2.0) - 6.0 * u * pow(v, 2.0) - 12.0 * pow(v, 2.0) * w - 2.0 * pow(v, 3.0)));
    shapefunction.set(2, 7, 1.0 / 12.0 * (4.0 * pow(u, 3.0) + 24.0 * pow(u, 2.0) * w + 24.0 * u * pow(w, 2.0) + 18.0 * pow(u, 2.0) * v + 48.0 * u * v * w + 12.0 * v * pow(w, 2.0) + 12.0 * u * pow(v, 2.0) + 12.0 * pow(v, 2.0) * w + 2.0 * pow(v, 3.0)));
    shapefunction.set(3, 7, -2.0 * u * w - 2.0 * pow(w, 2.0) - 2.0 * u * v - 2.0 * v * w);
    shapefunction.set(4, 7, pow(u, 2.0) - 2.0 * pow(w, 2.0) + u * v - 2.0 * v * w);
    shapefunction.set(5, 7, 1.0 / 2.0 * (pow(u, 2.0) - 2.0 * pow(w, 2.0) - 2.0 * u * v - 4.0 * v * w - pow(v, 2.0)));
    shapefunction.set(6, 7, 1.0 / 2.0 * (pow(u, 2.0) - 2.0 * pow(w, 2.0) - 2.0 * u * v - 4.0 * v * w - pow(v, 2.0)));

    shapefunction.set(0, 8, 1.0 / 12.0 * (2.0 * u * pow(w, 3.0) + pow(w, 4.0)));
    shapefunction.set(1, 8, -1.0 / 6.0 * pow(w, 3.0));
    shapefunction.set(2, 8, 1.0 / 12.0 * (6.0 * u * pow(w, 2.0) + 2.0 * pow(w, 3.0)));
    shapefunction.set(3, 8, 0.0);
    shapefunction.set(4, 8, u * w);
    shapefunction.set(5, 8, -1.0 / 2.0 * pow(w, 2.0));
    shapefunction.set(6, 8, -1.0 / 2.0 * pow(w, 2.0));

    shapefunction.set(0, 9, 1.0 / 12.0 * (2.0 * pow(v, 3.0) * w + pow(v, 4.0)));
    shapefunction.set(1, 9, 1.0 / 12.0 * (6.0 * pow(v, 2.0) * w + 4.0 * pow(v, 3.0)));
    shapefunction.set(2, 9, 1.0 / 6.0 * pow(v, 3.0));
    shapefunction.set(3, 9, v * w + pow(v, 2.0));
    shapefunction.set(4, 9, 0.0);
    shapefunction.set(5, 9, 1.0 / 2.0 * pow(v, 2.0));
    shapefunction.set(6, 9, 1.0 / 2.0 * pow(v, 2.0));

    shapefunction.set(0, 10, 1.0 / 12.0 * (2.0 * u * pow(w, 3.0) + pow(w, 4.0) + 6.0 * u * v * pow(w, 2.0) + 6.0 * v * pow(w, 3.0) + 6.0 * u * pow(v, 2.0) * w + 12.0 * pow(v, 2.0) * pow(w, 2.0) + 2.0 * u * pow(v, 3.0) + 6.0 * pow(v, 3.0) * w + pow(v, 4.0)));
    shapefunction.set(1, 10, 1.0 / 12.0 * (4.0 * pow(w, 3.0) + 18.0 * v * pow(w, 2.0) + 6.0 * u * pow(w, 2.0) + 12.0 * pow(v, 2.0) * w + 12.0 * u * v * w + 2.0 * pow(v, 3.0) + 6.0 * u * pow(v, 2.0)));
    shapefunction.set(2, 10, 1.0 / 12.0 * (2.0 * pow(w, 3.0) + 6.0 * u * pow(w, 2.0) + 12.0 * v * pow(w, 2.0) + 12.0 * u * v * w + 18.0 * pow(v, 2.0) * w + 6.0 * u * pow(v, 2.0) + 4.0 * pow(v, 3.0)));
    shapefunction.set(3, 10, pow(w, 2.0) + v * w + u * w + u * v);
    shapefunction.set(4, 10, u * w + v * w + u * v + pow(v, 2.0));
    shapefunction.set(5, 10, 1.0 / 2.0 * (pow(w, 2.0) + 4.0 * v * w + 2.0 * u * w + pow(v, 2.0) + 2.0 * u * v));
    shapefunction.set(6, 10, 1.0 / 2.0 * (pow(w, 2.0) + 4.0 * v * w + 2.0 * u * w + pow(v, 2.0) + 2.0 * u * v));

    shapefunction.set(0, 11, 1.0 / 12.0 * (pow(w, 4.0) + 2.0 * v * pow(w, 3.0)));
    shapefunction.set(1, 11, 1.0 / 6.0 * pow(w, 3.0));
    shapefunction.set(2, 11, 1.0 / 12.0 * (4.0 * pow(w, 3.0) + 6.0 * v * pow(w, 2.0)));
    shapefunction.set(3, 11, 0.0);
    shapefunction.set(4, 11, pow(w, 2.0) + v * w);
    shapefunction.set(5, 11, 1.0 / 2.0 * pow(w, 2.0));
    shapefunction.set(6, 11, 1.0 / 2.0 * pow(w, 2.0));
}

void get_shapefunction_vector(const Matrix &vwuMat, std::vector<Matrix> &sfVec)
{
    sfVec.resize(vwuMat.nrow(), Matrix(7, 12));
    // iterate over vwuMat
    for (int i = 0; i < vwuMat.nrow(); i++)
    {
        // get shapefunction
        get_shapefunction(vwuMat.get_row(i), sfVec[i]);
    }
}

void get_subdivision_matrices(Matrix &mat,
                              Matrix &subMat1,
                              Matrix &subMat2,
                              Matrix &subMat3,
                              Matrix &subMat4)
{
    const int N = 6;
    const double w = 3.0 / 8.0 / N; // w=1/N*(5/8-(3/8+1/4*cos(2*pi/N))^2);
    const int N1 = 5;
    const double w1 = 3.0 / 8.0 / N1; // w1=1/N1*(5/8-(3/8+1/4*cos(2*pi/N1))^2);
    const double a = 3.0 / 8.0;
    const double b = 1.0 / 8.0;
    std::vector<std::vector<double>> Mtmp{{a, b, a, b, 0, 0, 0, 0, 0, 0, 0},
                                {b, a, a, 0, 0, b, 0, 0, 0, 0, 0},
                                {w1, w1, 1.0 - N1 * w1, w1, 0, w1, w1, 0, 0, 0, 0},
                                {b, 0, a, a, 0, 0, b, 0, 0, 0, 0},
                                {0, a, b, 0, b, a, 0, 0, 0, 0, 0},
                                {0, b, a, 0, 0, a, b, 0, 0, 0, 0},
                                {0, 0, a, b, 0, b, a, 0, 0, 0, 0},
                                {0, 0, b, a, 0, 0, a, b, 0, 0, 0},
                                {0, b, 0, 0, a, a, 0, 0, b, 0, 0},
                                {0, w, w, 0, w, 1.0 - N * w, w, 0, w, w, 0},
                                {0, 0, b, 0, 0, a, a, 0, 0, b, 0},
                                {0, 0, w, w, 0, w, 1.0 - N * w, w, 0, w, w},
                                {0, 0, 0, b, 0, 0, a, a, 0, 0, b},
                                {0, 0, 0, 0, b, a, 0, 0, a, b, 0},
                                {0, 0, 0, 0, 0, a, b, 0, b, a, 0},
                                {0, 0, 0, 0, 0, b, a, 0, 0, a, b},
                                {0, 0, 0, 0, 0, 0, a, b, 0, b, a}};

    mat = Matrix(Mtmp);

    std::vector<std::vector<double>> SM1(12, std::vector<double>(17, 0.0));
    std::vector<int> element1{2, 3, 5, 6, 7, 9, 10, 11, 12, 14, 15, 16};
    for (int i = 0; i < 12; i++)
    {
        SM1[i][element1[i]] = 1.0;
    }
    subMat1 = Matrix(SM1);

    std::vector<std::vector<double>> SM2(12, std::vector<double>(17, 0.0));
    std::vector<int> element2{4, 1, 9, 5, 2, 14, 10, 6, 3, 15, 11, 7};
    for (int i = 0; i < 12; i++)
    {
        SM2[i][element2[i]] = 1.0;
    }
    subMat2 = Matrix(SM2);

    std::vector<std::vector<double>> SM3(12, std::vector<double>(17, 0.0));
    std::vector<int> element3{1, 2, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15};
    for (int i = 0; i < 12; i++)
    {
        SM3[i][element3[i]] = 1.0;
    }
    subMat3 = Matrix(SM3);

    std::vector<std::vector<double>> SM4(11, std::vector<double>(17, 0.0));
    std::vector<int> element4{0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11};
    for (int i = 0; i < 11; i++)
    {
        SM4[i][element4[i]] = 1.0;
    }
    subMat4 = Matrix(SM4);
}