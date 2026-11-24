"""Experiment units and the sample-name mapping between the original MAGs and the re-assembly units."""
UNITS = ["R1_A_D1", "R2_B_D2", "R2A_C_D1", "SF2_A_D1", "R1_B_H2O", "SF2_B_H2O", "R2_restored_D_shore", "SF2_sediment"]
SF2_SEDIMENT = ["SF2_A_D1", "SF2_A_D2", "SF2_B_D1", "SF2_B_D2", "SF2_C_D1", "SF2_C_D2"]


def orig_to_unit_sample(s):
    """Original-MAG sample name -> re-assembly sample name (only R2_restored_DShore differs)."""
    return "R2_restored_D_shore" if s == "R2_restored_DShore" else s


def unit_samples(unit):
    """The (re-assembly-named) samples whose original MAGs a unit should recover; a site co-assembly gets its 6 cores."""
    if unit.endswith("_sediment"):
        site = unit.split("_sediment")[0]
        return [f"{site}_{r}_{l}" for r in "ABC" for l in ("D1", "D2")]
    return [unit]
