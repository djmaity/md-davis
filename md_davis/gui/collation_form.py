"""Shared behavior for the GUI panels that build collation requests."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import wx

from md_davis.core.collation import create_hdf


def build_collation_data(form: Any) -> dict[str, Any]:
    """Build a collator request from controls exposed by a GUI form."""
    scalar_fields = {
        "name": form.name_textbox.GetValue(),
        "label": form.label_textbox.GetValue(),
        "text_label": form.text_label_textbox.GetValue(),
        "output": form.output_picker.GetPath(),
        "structure": form.structure_picker.GetPath(),
        "trajectory": form.trajectory_picker.GetPath(),
    }
    data = {key: value for key, value in scalar_fields.items() if value}

    grouped_fields = {
        "timeseries": {
            "rmsd": form.rmsd_picker.GetPath(),
            "rg": form.rg_picker.GetPath(),
        },
        "residue_property": {
            "secondary_structure": form.secstr_picker.GetPath(),
            "rmsf": form.rmsf_picker.GetPath(),
            "sasa": form.sasa_picker.GetPath(),
            "surface_potential": form.electrostatics_picker.GetPath(),
        },
    }
    for group, values in grouped_fields.items():
        selected = {key: value for key, value in values.items() if value}
        if selected:
            data[group] = selected

    if form.calculate_dihedral.GetValue():
        data["dihedral"] = {"chunk": form.chunk_intbox.GetValue()}
    return data


class CollationActionsMixin:
    """Actions shared by collation-oriented GUI panels."""

    def on_check_dihedral(self, event: wx.CommandEvent) -> None:
        visible = self.calculate_dihedral.GetValue()
        self.chunk_text.Show(visible)
        self.chunk_intbox.Show(visible)
        self.Layout()
        self.SetupScrolling(scrollToTop=False)

    def on_button_press(self, event: wx.CommandEvent) -> None:
        data = build_collation_data(self)
        if create_hdf(data):
            print(f"Collated all data into {self.output_picker.GetPath()}")
