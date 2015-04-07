from pyomo.core import *

import baa99_util

model = baa99_util.define_stochastic_model()

model.x1 = Var(bounds=(0,217))
model.x2 = Var(bounds=(0,217))
model.v1 = Var(within=NonNegativeReals)
model.v2 = Var(within=NonNegativeReals)
model.u1 = Var(within=NonNegativeReals)
model.u2 = Var(within=NonNegativeReals)
model.w11 = Var(within=NonNegativeReals)
model.w12 = Var(within=NonNegativeReals)
model.w22 = Var(within=NonNegativeReals)

model.FirstStageCost = \
    Expression(initialize=(4*model.x1 + 2*model.x2))
model.SecondStageCost = \
    Expression(initialize=(-8*model.w11 - 4*model.w12 - 4*model.w22 +\
                           0.2*model.v1 + 0.2*model.v2 + 10*model.u1 + 10*model.u2))

model.obj = Objective(expr=model.FirstStageCost + model.SecondStageCost)

model.d1 = Constraint(expr=model.w11 + model.u1 == model.d1_rhs)
model.PySP_StochasticRHS[model.d1_rhs] = model.d1

model.d2 = Constraint(expr=model.w12 + model.w22 + model.u2 == model.d2_rhs)
model.PySP_StochasticRHS[model.d2_rhs] = model.d2

model.s1 = Constraint(expr=-model.x1 + model.w11 + model.w12 + model.v1 == 0)

model.s2 = Constraint(expr=-model.x2 + model.w22 + model.v2 == 0)

num_scenarios = 100
sample_data = None
def pysp_scenario_tree_model_callback():
    global sample_data

    scenario_tree_model = \
        baa99_util.generate_scenario_tree_model(num_scenarios)

    sample_data = baa99_util.sample_into_scenario_tree_model(model, scenario_tree_model)

    return scenario_tree_model

def pysp_instance_creation_callback(scenario_name, node_names):
    assert sample_data is not None

    #
    # Clone a new instance and update the stochastic
    # parameters from the sampled scenario
    #

    instance = model.clone()

    for param_cuid, param_value in sample_data[scenario_name]:
        param_data = param_cuid.find_component(instance)
        param_data.value = param_value

    return instance
