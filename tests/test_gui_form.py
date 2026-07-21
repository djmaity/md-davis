"""Tests for GUI-independent form serialization."""

from types import SimpleNamespace

from md_davis.gui.collation_form import build_collation_data


class _Control:
    def __init__(self, value=""):
        self.value = value

    def GetPath(self):
        return self.value

    def GetValue(self):
        return self.value


def _form(**values):
    fields = {
        "name_textbox": "simulation",
        "label_textbox": "",
        "text_label_textbox": "",
        "output_picker": "output.h5",
        "structure_picker": "structure.pdb",
        "trajectory_picker": "trajectory.xtc",
        "rmsd_picker": "",
        "rg_picker": "",
        "secstr_picker": "",
        "rmsf_picker": "",
        "sasa_picker": "",
        "electrostatics_picker": "",
        "calculate_dihedral": False,
        "chunk_intbox": 1000,
    }
    fields.update(values)
    return SimpleNamespace(**{name: _Control(value) for name, value in fields.items()})


def test_build_collation_data_omits_empty_groups():
    data = build_collation_data(_form())

    assert "timeseries" not in data
    assert "residue_property" not in data
    assert "dihedral" not in data


def test_build_collation_data_uses_collator_keys():
    data = build_collation_data(
        _form(
            rmsd_picker="rmsd.xvg",
            rg_picker="rg.xvg",
            electrostatics_picker="potential",
            calculate_dihedral=True,
            chunk_intbox=250,
        )
    )

    assert data["timeseries"] == {"rmsd": "rmsd.xvg", "rg": "rg.xvg"}
    assert data["residue_property"] == {"surface_potential": "potential"}
    assert data["dihedral"] == {"chunk": 250}
