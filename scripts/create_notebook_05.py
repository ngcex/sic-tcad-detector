"""Create Phase 5 parametric studies notebook programmatically."""

import os

import nbformat

# Resolve paths relative to this script's location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
NOTEBOOK_PATH = os.path.join(PROJECT_DIR, "notebooks", "05_parametric_studies.ipynb")

nb = nbformat.v4.new_notebook()

# Cell 1: Title (markdown)
nb.cells.append(
    nbformat.v4.new_markdown_cell(
        "# Phase 5: Parametric Studies and Publication-Quality Results\n"
        "\n"
        "This notebook performs a multi-dimensional CCE vs dose-rate sweep across "
        "epitaxial layer thickness, bulk doping concentration, and bias voltage "
        "for the FLASH SiC paper.\n"
        "\n"
        "**Parameter space:** 4 epi thicknesses x 4 doping levels x 3 bias "
        "voltages = 48 conditions, each swept over 6 dose rates (20-230 Gy/s).\n"
        "\n"
        "**Scientific goal:** Determine whether any combination of detector "
        "design parameters leads to measurable CCE degradation under FLASH "
        "dose rates in 4H-SiC."
    )
)

# Cell 2: Imports (code)
nb.cells.append(
    nbformat.v4.new_code_cell(
        "import sys, os\n"
        "sys.path.insert(0, os.path.join(os.getcwd(), '..'))\n"
        "import numpy as np\n"
        "import matplotlib\n"
        "matplotlib.use('Agg')  # non-interactive backend for notebook creation\n"
        "import matplotlib.pyplot as plt\n"
        "from src.flash_recombination import (parametric_cce_sweep, save_parametric_results,\n"
        "                                      load_parametric_results, cce_vs_dose_rate)\n"
        "from src.plotting import (plot_parametric_epi, plot_parametric_doping, plot_parametric_bias,\n"
        "                          plot_cce_vs_dose_rate, save_figure)\n"
        "\n"
        "print('All imports successful')"
    )
)

# Cell 3: Configuration header (markdown)
nb.cells.append(
    nbformat.v4.new_markdown_cell(
        "## Configuration\n"
        "\n"
        "Modify these parameters to customize the parametric sweep. Set "
        "`RECOMPUTE = True` to run the full sweep (expect ~1-2 hours for 48 "
        "conditions). With `RECOMPUTE = False`, previously cached results are "
        "loaded from disk."
    )
)

# Cell 4: Configuration (code)
nb.cells.append(
    nbformat.v4.new_code_cell(
        "DOSE_RATES = np.array([20, 50, 100, 150, 200, 230], dtype=float)\n"
        "EPI_THICKNESSES = [5e-4, 10e-4, 15e-4, 20e-4]\n"
        "N_D_BULK_VALUES = [5e13, 1e14, 2e14, 5e14]\n"
        "BIAS_VOLTAGES = [-10.0, -30.0, -50.0]\n"
        "E_MEV = 62\n"
        "RECOMPUTE = False\n"
        'RESULTS_FILE = "../figures/parametric_results.json"'
    )
)

# Cell 5: Parametric sweep header (markdown)
nb.cells.append(
    nbformat.v4.new_markdown_cell(
        "## Parametric Sweep\n"
        "\n"
        "The sweep explores all combinations of the parameters above:\n"
        "- **Epi thickness:** 5, 10, 15, 20 um\n"
        "- **Bulk doping:** 5e13, 1e14, 2e14, 5e14 cm$^{-3}$\n"
        "- **Bias voltage:** -10, -30, -50 V\n"
        "\n"
        "Total: 4 x 4 x 3 = 48 parameter combinations, each swept over 6 dose rates."
    )
)

# Cell 6: Computation (code)
nb.cells.append(
    nbformat.v4.new_code_cell(
        "if RECOMPUTE:\n"
        "    results = parametric_cce_sweep(DOSE_RATES, EPI_THICKNESSES, N_D_BULK_VALUES, BIAS_VOLTAGES, E_MEV)\n"
        "    save_parametric_results(results, RESULTS_FILE)\n"
        '    print(f"Sweep complete: {len(results)} conditions saved to {RESULTS_FILE}")\n'
        "else:\n"
        "    results = load_parametric_results(RESULTS_FILE)\n"
        '    print(f"Loaded {len(results)} conditions from cache")'
    )
)

# Cell 7: Epi thickness header (markdown)
nb.cells.append(
    nbformat.v4.new_markdown_cell(
        "## Results: Epitaxial Thickness Dependence\n"
        "\n"
        "CCE vs dose rate for varying epitaxial layer thickness at fixed "
        "reference doping ($N_D$ = 8.50e13 cm$^{-3}$). One panel per bias voltage."
    )
)

# Cell 8: Epi thickness plot (code)
nb.cells.append(
    nbformat.v4.new_code_cell(
        "fig = plot_parametric_epi(results, EPI_THICKNESSES, N_D_bulk_ref=8.50e13)\n"
        'save_figure(fig, "flash_parametric_epi")\n'
        "plt.show()"
    )
)

# Cell 9: Doping header (markdown)
nb.cells.append(
    nbformat.v4.new_markdown_cell(
        "## Results: Doping Concentration Dependence\n"
        "\n"
        "CCE vs dose rate for varying bulk doping concentration at fixed "
        "reference epi thickness (10 um). One panel per bias voltage."
    )
)

# Cell 10: Doping plot (code)
nb.cells.append(
    nbformat.v4.new_code_cell(
        "fig = plot_parametric_doping(results, N_D_BULK_VALUES, epi_ref_cm=10e-4)\n"
        'save_figure(fig, "flash_parametric_doping")\n'
        "plt.show()"
    )
)

# Cell 11: Bias header (markdown)
nb.cells.append(
    nbformat.v4.new_markdown_cell(
        "## Results: Bias Voltage Dependence\n"
        "\n"
        "CCE vs dose rate for varying bias voltage at fixed reference epi "
        "thickness (10 um) and doping ($N_D$ = 8.50e13 cm$^{-3}$)."
    )
)

# Cell 12: Bias plot (code)
nb.cells.append(
    nbformat.v4.new_code_cell(
        "fig = plot_parametric_bias(results, BIAS_VOLTAGES, epi_ref_cm=10e-4, N_D_bulk_ref=8.50e13)\n"
        'save_figure(fig, "flash_parametric_bias")\n'
        "plt.show()"
    )
)

# Cell 13: Summary header (markdown)
nb.cells.append(
    nbformat.v4.new_markdown_cell(
        "## Summary and Key Findings\n"
        "\n"
        "This parametric study explores how detector design parameters affect "
        "charge collection under FLASH dose rates:\n"
        "\n"
        "1. **Number of parameter combinations:** 48 (4 epi x 4 doping x 3 bias), "
        "each swept over 6 dose rates from 20 to 230 Gy/s.\n"
        "\n"
        "2. **CCE degradation:** The analysis below checks whether any condition "
        "produces CCE below 0.99, indicating meaningful Auger-driven degradation.\n"
        "\n"
        "3. **Implications for SiC dosimeters under FLASH irradiation:** If CCE "
        "remains near unity across all conditions, 4H-SiC detectors are robust "
        "dose-rate-independent dosimeters regardless of design parameter choice "
        "within the explored range.\n"
        "\n"
        "---\n"
        "*Phase 5 of the Petringa 4H-SiC TCAD simulation project*"
    )
)

# Cell 14: Summary statistics (code)
nb.cells.append(
    nbformat.v4.new_code_cell(
        "total = len(results)\n"
        "failed = sum(1 for v in results.values() if v is None)\n"
        'print(f"Parameter combinations: {total}")\n'
        'print(f"Failed simulations: {failed}")\n'
        'print(f"Successful: {total - failed}")\n'
        "# Check if any CCE < 0.99 (indicating meaningful degradation)\n"
        "for key, val in results.items():\n"
        '    if val is not None and np.any(np.array(val["cce_values"]) < 0.99):\n'
        "        epi, nd, vb = key\n"
        "        print(f\"  CCE degradation at epi={epi*1e4:.0f}um, N_D={nd:.1e}, V={vb}V: min CCE={min(val['cce_values']):.4f}\")"
    )
)

# Write the notebook
nbformat.write(nb, NOTEBOOK_PATH)
print(f"Notebook created with {len(nb.cells)} cells")
print(f"Code cells: {len([c for c in nb.cells if c.cell_type == 'code'])}")
print(f"Markdown cells: {len([c for c in nb.cells if c.cell_type == 'markdown'])}")
