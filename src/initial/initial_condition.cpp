/**
 * @file initial_condition.cpp
 * @brief Initial condition factory implementation
 */

#include "euler1d/initial/initial_condition.hpp"

namespace euler1d {

InitialConditionVariant create_initial_condition(const InitialConditionConfig& config) {
    switch (config.type) {
        case InitialConditionType::PiecewiseConstant: {
            PiecewiseConstantIC ic;
            ic.regions = config.regions;
            return ic;
        }
        case InitialConditionType::ShockEntropyInteraction: {
            ShockEntropyInteractionIC ic;
            ic.discontinuity_position = config.discontinuity_position;
            ic.left_state = config.left_state;
            ic.right_state = config.right_state;
            return ic;
        }
    }
    // Should never reach here
    return PiecewiseConstantIC{};
}

}  // namespace euler1d
