import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tc_python import *

"""
This example shows how to create a property (step) diagram using TC-Python.
The alloy system Fe-Ni is used as an example. 
"""

df_lst = []

with TCPython() as start:
    start.set_cache_folder(os.path.basename(__file__) + "_cache")

    # temp_celsius = 900
    temp_step = 100
    for temp_celsius in np.arange(500, 1200 + temp_step, temp_step):
        calculation = (
            start.select_database_and_elements("TCCU5", ["Cu", "Ag"]).without_default_phases().select_phase("LIQUID").select_phase(
                "FCC_A1").get_system().
            with_property_diagram_calculation().
            with_axis(CalculationAxis(ThermodynamicQuantity.mass_fraction_of_a_component("Ag")).
                      set_min(0).
                      set_max(1).
                      with_axis_type(Linear().set_min_nr_of_steps(1000))).enable_step_separate_phases().
            set_condition(ThermodynamicQuantity.temperature(), temp_celsius + 273.15)
        )

        property_diagram = calculation.calculate()
        # property_diagram.set_phase_name_style(PhaseNameStyle.ALL)
        groups_liquid = property_diagram.get_values_grouped_by_quantity_of(
            ThermodynamicQuantity.mass_fraction_of_a_component("Ag"),
            ThermodynamicQuantity.gibbs_energy_of_a_phase("LIQUID", use_ser=True)
        )
        groups_fcc = property_diagram.get_values_grouped_by_quantity_of(
            ThermodynamicQuantity.mass_fraction_of_a_component("Ag"),
            ThermodynamicQuantity.gibbs_energy_of_a_phase("FCC_A1", use_ser=True)
        )

        for group_liquid, group_fcc in zip(groups_liquid.values(), groups_fcc.values()):
            liquid_x, liquid_y = np.array(group_liquid.x), np.array(group_liquid.y)
            fcc_x, fcc_y = np.array(group_fcc.x), np.array(group_fcc.y)

            idx_to_keep_liquid = [idx for idx, val in enumerate(liquid_y) if val != 0 and not np.isnan(val)]
            idx_to_keep_fcc = [idx for idx, val in enumerate(fcc_y) if val != 0 and not np.isnan(val)]

            liquid_x, liquid_y = liquid_x[idx_to_keep_liquid], liquid_y[idx_to_keep_liquid]
            fcc_x, fcc_y = fcc_x[idx_to_keep_fcc], fcc_y[idx_to_keep_fcc]

            df_liquid = pd.DataFrame.from_dict(
                {
                    "x": liquid_x,
                    "y": liquid_y,
                    "phase": "LIQUID",
                    "temp": temp_celsius
                }
            )
            df_lst.append(df_liquid)
            df_fcc = pd.DataFrame.from_dict(
                {
                    "x": fcc_x,
                    "y": fcc_y,
                    "phase": "FCC",
                    "temp": temp_celsius
                }
            )
            df_lst.append(df_fcc)

df_to_store = pd.concat(df_lst)

df_to_store.to_csv(r"C:\Users\rpw19\PycharmProjects\matsci_animation\data\gibbs_energy\binary_Cu_Ag.csv", index=False)

# plt.xlabel("Cr [wt fraction]")
# plt.ylabel("Gibbs energy")
# plt.legend(loc="upper left")
# # plt.title("Fe-10Ni")
# plt.tight_layout()
# plt.show()
