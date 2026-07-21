from __future__ import annotations

from pathlib import Path

import wx

from md_davis.gui.collate_panel import CollatePanel
from md_davis.gui.electrodynamics_panel import ElectrodynamicsPanel
from md_davis.gui.landscape_panel import LandscapePanel
from md_davis.gui.sequence_panel import SequencePanel


class MainPanel(wx.Panel):
    def __init__(self, parent: wx.Window) -> None:
        super().__init__(parent)

        main_sizer = wx.BoxSizer(orient=wx.VERTICAL)
        notebook = wx.Notebook(self)
        notebook.AddPage(LandscapePanel(notebook), "Landscape")
        notebook.AddPage(ElectrodynamicsPanel(notebook), "Electrodynamics")
        notebook.AddPage(CollatePanel(notebook), "Collate")
        notebook.AddPage(SequencePanel(notebook), "Sequence")
        main_sizer.Add(notebook, 1, wx.ALL | wx.EXPAND, 5)
        self.SetSizer(main_sizer)
        self.Layout()


class MainFrame(wx.Frame):
    def __init__(self) -> None:
        super().__init__(None, title="MD DaVis")
        self.SetIcon(wx.Icon(str(Path(__file__).with_name("MD_DaVis.ico"))))
        self.panel = MainPanel(self)
        self.SetMinSize(wx.Size(width=600, height=400))
        self.Center()


def main() -> None:
    app = wx.App()
    frame = MainFrame()
    frame.Show()
    app.MainLoop()


if __name__ == "__main__":
    main()
